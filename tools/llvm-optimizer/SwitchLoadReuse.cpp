#include "SwitchLoadReuse.h"
#include "llvm/ADT/STLExtras.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/Analysis/OptimizationRemarkEmitter.h"
#include "llvm/Analysis/TargetTransformInfo.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/Module.h"

using namespace llvm;

namespace rssn {
namespace {

constexpr char PassName[] = "rssn-switch-load-reuse";
constexpr unsigned MaxArms = 16;
constexpr unsigned MaxPrefix = 12;

bool isIntegerRead(const LoadInst &L) {
  if (!L.isSimple() || L.getPointerAddressSpace() != 0 || L.use_empty() ||
      !isa<IntegerType>(L.getType()))
    return false;
  switch (L.getType()->getIntegerBitWidth()) {
  case 8:
  case 16:
  case 32:
  case 64:
    break;
  default:
    return false;
  }
  SmallVector<std::pair<unsigned, MDNode *>, 4> Metadata;
  L.getAllMetadataOtherThanDebugLoc(Metadata);
  return all_of(Metadata, [](const auto &Entry) {
    switch (Entry.first) {
    case LLVMContext::MD_tbaa:
    case LLVMContext::MD_alias_scope:
    case LLVMContext::MD_noalias:
    case LLVMContext::MD_noundef:
      return true;
    default:
      return false;
    }
  });
}

struct Reuse {
  SmallVector<LoadInst *, 16> Wide;
  SmallVector<LoadInst *, 8> Narrow;
  Align Alignment;
};

Reuse findReuse(SwitchInst &Switch, const DominatorTree &DT,
                const TargetTransformInfo &TTI) {
  if (Switch.getNumSuccessors() > MaxArms)
    return {};
  SmallVector<SmallVector<LoadInst *, 4>, 8> Reads;
  SmallPtrSet<BasicBlock *, 8> Seen;
  for (BasicBlock *Arm : successors(Switch.getParent())) {
    if (!Seen.insert(Arm).second)
      continue;
    if (Arm->getUniquePredecessor() != Switch.getParent() ||
        Arm == Switch.getParent() || Arm->isEHPad())
      return {};
    auto &Prefix = Reads.emplace_back();
    unsigned Remaining = MaxPrefix;
    for (Instruction &I : *Arm) {
      if (isa<DbgInfoIntrinsic>(I))
        continue;
      if (I.isTerminator() || Remaining-- == 0 || I.isEHPad() ||
          isa<CallBase>(I) || I.mayHaveSideEffects() ||
          !isGuaranteedToTransferExecutionToSuccessor(&I))
        break;
      if (auto *L = dyn_cast<LoadInst>(&I)) {
        if (!L->isSimple())
          break;
        if (isIntegerRead(*L) && DT.dominates(L->getPointerOperand(), &Switch))
          Prefix.push_back(L);
      }
    }
    if (Prefix.empty())
      return {};
  }
  if (Reads.size() < 3)
    return {};

  for (LoadInst *First : Reads.front()) {
    Type *WideTy = First->getType();
    if (WideTy->isIntegerTy(8) || !TTI.isTypeLegal(WideTy))
      continue;
    const Value *Pointer = First->getPointerOperand();
    Reuse Plan{{}, {}, First->getAlign()};
    for (const auto &Arm : Reads) {
      auto Wide = find_if(Arm, [&](LoadInst *L) {
        return L->getPointerOperand() == Pointer && L->getType() == WideTy;
      });
      if (Wide == Arm.end())
        break;
      Plan.Wide.push_back(*Wide);
      Plan.Alignment = std::min(Plan.Alignment, (*Wide)->getAlign());
      for (LoadInst *L : Arm)
        if (L->getPointerOperand() == Pointer &&
            L->getType()->getIntegerBitWidth() < WideTy->getIntegerBitWidth() &&
            TTI.isTruncateFree(WideTy, L->getType()))
          Plan.Narrow.push_back(L);
    }
    if (Plan.Wide.size() == Reads.size() && !Plan.Narrow.empty())
      return Plan;
  }
  return {};
}

void reuseLoads(SwitchInst &Switch, const Reuse &Plan) {
  IRBuilder<> B(&Switch);
  auto *First = Plan.Wide.front();
  auto *WideTy = cast<IntegerType>(First->getType());
  auto *Bytes = B.getByteNTy(WideTy->getBitWidth());
  // A byte type carries poison per bit without forcing the payload into
  // a SIMD register. Project before converting either width to integer.
  auto *Read = B.CreateAlignedLoad(Bytes, First->getPointerOperand(),
                                   Plan.Alignment, "switch.bytes");
  auto *Wide = B.CreateBitCast(Read, WideTy, "switch.word");

  for (LoadInst *L : Plan.Wide) {
    L->replaceAllUsesWith(Wide);
    L->eraseFromParent();
  }
  for (LoadInst *L : Plan.Narrow) {
    B.SetInsertPoint(L);
    auto *Ty = cast<IntegerType>(L->getType());
    auto *PartTy = B.getByteNTy(Ty->getBitWidth());
    auto *PartsTy =
        FixedVectorType::get(PartTy, WideTy->getBitWidth() / Ty->getBitWidth());
    auto *Parts = B.CreateBitCast(Read, PartsTy, "switch.parts");
    auto *Low = B.CreateExtractElement(Parts, B.getInt32(0), "switch.low");
    auto *Part = B.CreateBitCast(Low, Ty, "switch.narrow");
    L->replaceAllUsesWith(Part);
    L->eraseFromParent();
  }
}

} // namespace

PreservedAnalyses SwitchLoadReusePass::run(Function &F,
                                           FunctionAnalysisManager &AM) {
  const Module &M = *F.getParent();
  if (M.getTargetTriple().getArch() != Triple::x86_64 ||
      !M.getDataLayout().isLittleEndian() || !F.getName().contains(Filter) ||
      F.isDeclaration() || F.hasOptNone())
    return PreservedAnalyses::all();

  bool Changed = false;
  for (BasicBlock &B : F) {
    auto *Switch = dyn_cast<SwitchInst>(B.getTerminator());
    if (!Switch || Switch->getNumSuccessors() < 3)
      continue;
    const auto &DT = AM.getResult<DominatorTreeAnalysis>(F);
    if (!DT.isReachableFromEntry(&B))
      continue;
    const auto &TTI = AM.getResult<TargetIRAnalysis>(F);
    const Reuse Plan = findReuse(*Switch, DT, TTI);
    if (Plan.Narrow.empty())
      continue;
    auto &ORE = AM.getResult<OptimizationRemarkEmitterAnalysis>(F);
    ORE.emit([&] {
      using NV = ore::NV;
      return OptimizationRemarkAnalysis(PassName, "SwitchReuse", Switch)
             << NV("Function", F.getName()) << ", common read in "
             << NV("Arms", Plan.Wide.size()) << " switch arms covers "
             << NV("CoveredReads", Plan.Narrow.size()) << " narrower reads"
             << (Transform ? ", reused" : ", eligible");
    });
    if (Transform) {
      reuseLoads(*Switch, Plan);
      Changed = true;
    }
  }

  if (!Changed)
    return PreservedAnalyses::all();
  PreservedAnalyses PA;
  PA.preserveSet<CFGAnalyses>();
  return PA;
}

} // namespace rssn
