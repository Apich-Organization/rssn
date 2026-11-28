#include "ScalarLoadPack.h"
#include "Planning.h"
#include "llvm/ADT/MapVector.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/Analysis/AliasAnalysis.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/OptimizationRemarkEmitter.h"
#include "llvm/Analysis/TargetTransformInfo.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/Module.h"
#include <optional>

using namespace llvm;

namespace rssn {
namespace {

constexpr char PassName[] = "rssn-scalar-load-pack";
constexpr unsigned MaxWindow = 32;
constexpr unsigned MinFields = 3;
constexpr unsigned MaxFields = 4;
constexpr unsigned MaxAliases = 16;
constexpr unsigned MaxPacks = 24;

struct Field {
  LoadInst *Load;
  int64_t Offset;
  unsigned Position;
};

struct Pack {
  SmallVector<Field, MaxFields> Fields;
  unsigned Bytes;
  unsigned Begin;
  unsigned End;
  bool WholeRecord;
};

bool isFieldLoad(const LoadInst &L) {
  if (!L.isSimple() || L.getPointerAddressSpace() != 0 ||
      (!L.getType()->isIntegerTy(8) && !L.getType()->isIntegerTy(16) &&
       !L.getType()->isIntegerTy(32)) || L.use_empty())
    return false;
  SmallVector<std::pair<unsigned, MDNode *>, 4> Metadata;
  L.getAllMetadataOtherThanDebugLoc(Metadata);
  return all_of(Metadata, [](const auto &MD) {
    return MD.first == LLVMContext::MD_tbaa ||
           MD.first == LLVMContext::MD_alias_scope ||
           MD.first == LLVMContext::MD_noalias ||
           MD.first == LLVMContext::MD_noundef;
  });
}

bool isFourByteRecord(const Value *Base, const DataLayout &DL) {
  // A single dynamic index into a four-byte element is distinct from
  // the low half of an eight-byte entry with an independently used dword.
  const auto *GEP = dyn_cast<GEPOperator>(Base);
  return GEP && GEP->getNumIndices() == 1 &&
         !isa<Constant>(GEP->idx_begin()->get()) &&
         DL.getTypeAllocSize(GEP->getSourceElementType()) == TypeSize::getFixed(4);
}

bool hasDependentRecordReads(ArrayRef<Pack> Packs, ArrayRef<Instruction *> Code) {
  SmallPtrSet<const LoadInst *, 32> Fields;
  for (const Pack &P : Packs)
    if (P.WholeRecord)
      for (const Field &F : P.Fields)
        Fields.insert(F.Load);

  SmallPtrSet<const Value *, 32> FromField;
  for (const Instruction *I : Code) {
    if (const auto *Load = dyn_cast<LoadInst>(I); Load && Fields.contains(Load)) {
      if (FromField.contains(Load->getPointerOperand()))
        return true;
      FromField.insert(Load);
    } else if (!isa<PHINode, CallBase>(I) && !I->mayReadOrWriteMemory() &&
               any_of(I->operand_values(),
                      [&](const Value *V) { return FromField.contains(V); })) {
      FromField.insert(I);
    }
  }
  return false;
}

std::optional<Pack> identifyPack(ArrayRef<Field> Fields,
                                 ArrayRef<Instruction *> Code,
                                 const DominatorTree &DT, MemoryModel &Memory,
                                 const TargetTransformInfo &TTI, unsigned Bytes,
                                 bool WholeRecord) {
  if (Fields.size() < MinFields)
    return std::nullopt;
  const Field &Low = Fields.front();
  int64_t Expected = Low.Offset;
  for (const Field &F : Fields) {
    const int64_t Width = F.Load->getType()->getIntegerBitWidth() / 8;
    if (F.Offset != Expected ||
        AddOverflow(Expected, Width, Expected))
      return std::nullopt;
  }
  if (Expected - Low.Offset != Bytes)
    return std::nullopt;
  const auto [First, Last] = std::minmax_element(
      Fields.begin(), Fields.end(),
      [](const Field &A, const Field &B) { return A.Position < B.Position; });
  if (Last->Position - First->Position > MaxWindow ||
      !DT.dominates(Low.Load->getPointerOperand(), First->Load))
    return std::nullopt;
  unsigned Fast = 0;
  if (Low.Load->getAlign() < Align(Bytes) &&
      (!TTI.allowsMisalignedMemoryAccesses(Low.Load->getContext(), 8 * Bytes,
                                           0, Low.Load->getAlign(), &Fast) || !Fast))
    return std::nullopt;
  for (unsigned K = First->Position; K <= Last->Position; ++K) {
    const Instruction &I = *Code[K];
    if (isa<CallBase>(I) || !isGuaranteedToTransferExecutionToSuccessor(&I))
      return std::nullopt;
  }
  for (const Field &F : Fields)
    if (!Memory.unchanged(*F.Load, Code, First->Position, F.Position))
      return std::nullopt;
  return Pack{SmallVector<Field, MaxFields>(Fields), Bytes, First->Position,
              Last->Position, WholeRecord};
}

void packLoads(const Pack &Plan) {
  const auto First = std::min_element(Plan.Fields.begin(), Plan.Fields.end(),
      [](const Field &A, const Field &B) { return A.Position < B.Position; });
  LoadInst *Low = Plan.Fields.front().Load;
  IRBuilder<> B(First->Load);
  auto *Bytes = FixedVectorType::get(B.getInt8Ty(), Plan.Bytes);
  // An integer wide load would make poison in one field observable in an
  // unrelated field. Keep byte lanes independent until the field is extracted.
  auto *Wide = B.CreateAlignedLoad(Bytes, Low->getPointerOperand(), Low->getAlign(),
                                   "scalar.pack");
  if (all_of(Plan.Fields, [](const Field &F) {
        return F.Load->hasMetadata(LLVMContext::MD_noundef);
      }))
    Wide->setMetadata(LLVMContext::MD_noundef, Low->getMetadata(LLVMContext::MD_noundef));

  for (const Field &F : Plan.Fields) {
    B.SetInsertPoint(F.Load);
    B.SetCurrentDebugLocation(F.Load->getDebugLoc());
    const auto Offset = static_cast<unsigned>(F.Offset - Plan.Fields.front().Offset);
    const unsigned Count = F.Load->getType()->getIntegerBitWidth() / 8;
    Value *Extract;
    if (Count == 1) {
      Extract = B.CreateExtractElement(Wide, Offset, "scalar.field");
    } else {
      SmallVector<int, 4> Mask;
      for (unsigned I = 0; I < Count; ++I)
        Mask.push_back(Offset + I);
      Value *Slice = B.CreateShuffleVector(Wide, Mask, "scalar.bytes");
      Extract = B.CreateBitCast(Slice, F.Load->getType(), "scalar.field");
    }
    F.Load->replaceAllUsesWith(Extract);
  }
  for (const Field &F : Plan.Fields)
    F.Load->eraseFromParent();
}

} // namespace

PreservedAnalyses ScalarLoadPackPass::run(Function &F, FunctionAnalysisManager &AM) {
  const Triple &T = F.getParent()->getTargetTriple();
  if (T.getArch() != Triple::x86_64 ||
      (!T.isOSWindows() && (!T.isOSLinux() || T.isX32())) ||
      !F.getParent()->getDataLayout().isLittleEndian() || F.isDeclaration() ||
      !F.getName().contains(Filter) || F.hasOptNone() || F.hasOptSize() ||
      F.hasMinSize() || !Liveness::supports(F))
    return PreservedAnalyses::all();
  const auto &LI = AM.getResult<LoopAnalysis>(F);
  if (LI.empty())
    return PreservedAnalyses::all();
  const auto &TTI = AM.getResult<TargetIRAnalysis>(F);
  const auto &DT = AM.getResult<DominatorTreeAnalysis>(F);
  auto &ORE = AM.getResult<OptimizationRemarkEmitterAnalysis>(F);
  MemoryModel Memory(AM.getResult<AAManager>(F));
  std::optional<Liveness> Live;

  const auto &DL = F.getParent()->getDataLayout();
  SmallVector<Pack, MaxPacks> Work;
  for (BasicBlock &B : F) {
    if (!LI.getLoopFor(&B) || B.isEHPad() || !DT.isReachableFromEntry(&B))
      continue;
    SmallVector<Instruction *, 64> Code;
    MapVector<const Value *, SmallVector<Field, 4>> ByBase;
    for (auto [Position, I] : enumerate(B)) {
      Code.push_back(&I);
      if (auto *L = dyn_cast<LoadInst>(&I); L && isFieldLoad(*L))
        if (auto A = accessAddress(*L))
          ByBase[A->first].push_back({L, A->second, static_cast<unsigned>(Position)});
    }
    SmallVector<Pack, 8> Hits;
    for (auto &[Base, Fields] : ByBase) {
      if (Fields.size() < MinFields || Fields.size() > MaxAliases)
        continue;
      const bool FourByteRecord = isFourByteRecord(Base, DL);
      sort(Fields, [](const Field &A, const Field &B) { return A.Offset < B.Offset; });
      for (unsigned First = 0; First + MinFields <= Fields.size();) {
        const bool WholeRecord = FourByteRecord && Fields[First].Offset == 0;
        const unsigned Bytes = T.isOSWindows() && !WholeRecord ? 8 : 4;
        bool Packed = false;
        const unsigned Remaining = Fields.size() - First;
        for (unsigned Count = MinFields; Count <= std::min(Remaining, MaxFields); ++Count) {
          if (auto P = identifyPack(ArrayRef(Fields).slice(First, Count), Code, DT, Memory,
                                    TTI, Bytes, WholeRecord)) {
            Hits.push_back(std::move(*P));
            First += Count;
            Packed = true;
            break;
          }
        }
        if (!Packed)
          ++First;
      }
    }
    if (Hits.size() < 2)
      continue;
    // Only functions with legal packs reach the function-wide dataflow solve.
    if (!Live) {
      Live.emplace(TTI, F);
      if (!Live->compute(F))
        return PreservedAnalyses::all();
    }
    const BlockPressure BP(B, *Live, false, Code);
    unsigned Peak = 0;
    unsigned Reads = 0;
    unsigned RecordReads = 0;
    uint64_t RecordGP = 0;
    for (const Pack &P : Hits) {
      for (unsigned K = P.Begin; K <= P.End; ++K)
        Peak = std::max(Peak, static_cast<unsigned>(BP.Before[K].GP));
      Reads += P.Fields.size();
      if (P.WholeRecord) {
        RecordReads += P.Fields.size();
        for (const Field &F : P.Fields)
          RecordGP += static_cast<unsigned>(BP.Before[F.Position].GP);
      }
    }
    const unsigned Budget = static_cast<unsigned>(Live->Limit.GP);
    // A four-byte load replaces two memory operations with field extraction.
    // When a field feeds the address of a later entry, extraction can lengthen
    // this serial chain. Require pressure across the field reads to justify it.
    const bool DependentReads = RecordReads && hasDependentRecordReads(Hits, Code);
    const bool RecordProfitable = !DependentReads || RecordGP > 2 * Budget * RecordReads;
    ORE.emit([&] {
      using NV = ore::NV;
      auto Remark = OptimizationRemarkAnalysis(PassName, "ScalarLoadPack", B.getTerminator())
             << NV("Function", F.getName()) << " " << NV("Action", Transform ? "can pack" : "found")
             << " " << NV("Loads", Reads) << " scalar loads into " << NV("Packs", Hits.size())
             << " packs";
      const unsigned Bytes = Hits.front().Bytes;
      if (all_of(Hits, [Bytes](const Pack &P) { return P.Bytes == Bytes; }))
        Remark << " of " << NV("BytesPerPack", Bytes) << " bytes";
      else
        Remark << " of mixed widths";
      Remark << ", GP pressure " << NV("Peak", Peak)
             << " / " << NV("Budget", Live->Limit.GP);
      if (RecordReads)
        Remark << ", four-byte records have mean GP " << NV("RecordGP", RecordGP / RecordReads)
               << " (" << (DependentReads ? "dependent reads, " : "independent reads, ")
               << (RecordProfitable ? "eligible" : "below rewrite threshold") << ")";
      return Remark;
    });
    if (!Transform || Peak < Budget)
      continue;
    if (!RecordProfitable)
      llvm::erase_if(Hits, [](const Pack &P) { return P.WholeRecord; });
    if (Hits.size() < 2 || Work.size() + Hits.size() > MaxPacks)
      continue;
    append_range(Work, Hits);
  }
  for (const Pack &P : Work)
    packLoads(P);
  if (Work.empty())
    return PreservedAnalyses::all();
  PreservedAnalyses PA;
  PA.preserveSet<CFGAnalyses>();
  return PA;
}

} // namespace rssn
