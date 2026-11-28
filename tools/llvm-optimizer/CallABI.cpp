#include "CallABI.h"
#include "Pressure.h"
#include "llvm/ADT/DenseSet.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/Analysis/BlockFrequencyInfo.h"
#include "llvm/Analysis/OptimizationRemarkEmitter.h"
#include "llvm/Analysis/TargetTransformInfo.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/InlineAsm.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/Module.h"
#include <limits>
#include <optional>

using namespace llvm;

namespace rssn {
namespace {

constexpr char RemarkName[] = "rssn-call-abi";
constexpr unsigned SysVGP = 6;
constexpr unsigned RegCallVolatileGP = 7;
constexpr unsigned RegCallGP = 11;
constexpr int PreservedGPBudget = 5;
constexpr int FrameSavePair = 2;
constexpr int StackArgumentAccess = 2;
constexpr int CrossCallSpillPair = 4;

enum class Bank { Unknown, GP, SSE };

Bank bank(Type *T) {
  if ((T->isPointerTy() && T->getPointerAddressSpace() == 0) ||
      T->isIntegerTy(1) || T->isIntegerTy(8) || T->isIntegerTy(16) ||
      T->isIntegerTy(32) || T->isIntegerTy(64))
    return Bank::GP;
  if (T->isFloatTy() || T->isDoubleTy())
    return Bank::SSE;
  return Bank::Unknown;
}

bool originalConvention(CallingConv::ID CC) {
  return CC == CallingConv::C || CC == CallingConv::Fast ||
         CC == CallingConv::X86_64_SysV;
}

bool canUseSSE(const Function &F) {
  const auto TargetFeatures = F.getFnAttribute("target-features");
  const StringRef Features = TargetFeatures.isStringAttribute()
                                 ? TargetFeatures.getValueAsString()
                                 : StringRef();
  SmallVector<StringRef, 16> Tokens;
  Features.split(Tokens, ',', -1, false);
  return !is_contained(Tokens, "-sse") && !is_contained(Tokens, "-sse2") &&
         !F.hasFnAttribute(Attribute::NoImplicitFloat);
}

struct Layout {
  unsigned GP = 0;
  unsigned SSE = 0;
  unsigned UsedStackGP = 0;
  bool Supported = true;
};

bool signatureAttrs(AttributeSet Attrs, unsigned I, Type *Return) {
  for (Attribute::AttrKind K :
       {Attribute::ByVal, Attribute::ByRef, Attribute::InAlloca,
        Attribute::Preallocated, Attribute::InReg, Attribute::Nest,
        Attribute::SwiftSelf, Attribute::SwiftError, Attribute::SwiftAsync,
        Attribute::Returned, Attribute::ImmArg})
    if (Attrs.hasAttribute(K))
      return false;
  return !Attrs.hasAttribute(Attribute::StructRet) ||
         (I == 0 && Return->isVoidTy());
}

Layout layout(const Function &F) {
  Layout L;
  if (F.isVarArg() || F.hasFnAttribute(Attribute::AllocSize) ||
      (!F.getReturnType()->isVoidTy() &&
       bank(F.getReturnType()) == Bank::Unknown))
    L.Supported = false;
  for (const Argument &A : F.args()) {
    if (!signatureAttrs(F.getAttributes().getParamAttrs(A.getArgNo()),
                        A.getArgNo(), F.getReturnType()))
      L.Supported = false;
    switch (bank(A.getType())) {
    case Bank::Unknown:
      L.Supported = false;
      break;
    case Bank::GP:
      if (L.GP++ >= SysVGP && !A.use_empty())
        ++L.UsedStackGP;
      break;
    case Bank::SSE:
      ++L.SSE;
      break;
    }
  }
  return L;
}

bool ordinaryLeaf(const Function &F) {
  unsigned Instructions = 0;
  auto mayLowerOutOfLine = [](Type *Ty) {
    Type *Scalar = Ty->getScalarType();
    return (Scalar->isIntegerTy() && Scalar->getIntegerBitWidth() > 64) ||
           (Scalar->isFloatingPointTy() && !Scalar->isFloatTy() &&
            !Scalar->isDoubleTy());
  };
  for (const BasicBlock &B : F)
    for (const Instruction &I : B) {
      if (++Instructions > 192 || mayLowerOutOfLine(I.getType()) ||
          any_of(I.operand_values(), [&](const Value *V) {
            return mayLowerOutOfLine(V->getType());
          }))
        return false;
      if (isa<AllocaInst, VAArgInst, AtomicRMWInst, AtomicCmpXchgInst>(I))
        return false;
      switch (I.getOpcode()) {
      case Instruction::FRem:
        return false; // May become a libcall without an IR call instruction.
      case Instruction::SDiv:
      case Instruction::UDiv:
      case Instruction::SRem:
      case Instruction::URem:
        if (!I.getType()->isIntegerTy() ||
            I.getType()->getIntegerBitWidth() > 64)
          return false;
        break;
      default:
        break;
      }
      if (auto *II = dyn_cast<IntrinsicInst>(&I)) {
        const Intrinsic::ID ID = II->getIntrinsicID();
        if (isa<DbgInfoIntrinsic>(II) || II->isLifetimeStartOrEnd() ||
            ID == Intrinsic::assume ||
            ((ID == Intrinsic::fshl || ID == Intrinsic::fshr) &&
             (II->getType()->isIntegerTy(32) ||
              II->getType()->isIntegerTy(64))))
          continue;
        return false;
      }
      if (isa<CallBase>(I))
        return false;
    }
  return true;
}

bool likelySiblingCall(const CallInst &C) {
  if (!C.isTailCall())
    return false;
  const Instruction *Next = C.getNextNode();
  while (isa_and_nonnull<DbgInfoIntrinsic>(Next))
    Next = Next->getNextNode();
  const auto *Ret = dyn_cast_or_null<ReturnInst>(Next);
  return Ret && (Ret->getReturnValue() == &C || C.getType()->isVoidTy());
}

bool collectSites(Function &F, SmallVectorImpl<CallInst *> &Sites,
                  ArrayRef<StringRef> Assembly) {
  if (any_of(Assembly,
             [&](StringRef Asm) { return Asm.contains(F.getName()); }))
    return false;
  for (Use &U : F.uses()) {
    auto *C = dyn_cast<CallInst>(U.getUser());
    if (!C || !C->isCallee(&U) || C->getFunctionType() != F.getFunctionType() ||
        C->getCallingConv() != F.getCallingConv() || C->isMustTailCall() ||
        likelySiblingCall(*C) || !C->doesNotThrow() || C->hasOperandBundles() ||
        C->isConvergent() || C->hasFnAttr(Attribute::AllocSize) ||
        C->hasFnAttr(Attribute::ReturnsTwice) ||
        C->getFunction()->hasOptNone() || C->getFunction()->hasOptSize() ||
        C->getFunction()->hasMinSize() || !canUseSSE(*C->getFunction()))
      return false;
    for (unsigned I = 0; I < C->arg_size(); ++I) {
      const auto Attrs = C->getAttributes().getParamAttrs(I);
      if (!signatureAttrs(Attrs, I, C->getType()) ||
          Attrs.hasAttribute(Attribute::StructRet) !=
              F.getArg(I)->hasStructRetAttr())
        return false;
    }
    Sites.push_back(C);
    if (Sites.size() > 32)
      return false;
  }
  return !Sites.empty();
}

struct CallerFacts {
  DenseMap<const Instruction *, Pressure> Through;
  bool Complete = false;
};

class Model {
  FunctionAnalysisManager &FAM;
  DenseMap<const Function *, CallerFacts> Callers;

public:
  explicit Model(FunctionAnalysisManager &FAM) : FAM(FAM) {}

  std::optional<Pressure> through(CallInst &C) {
    Function &F = *C.getFunction();
    auto [It, New] = Callers.try_emplace(&F);
    CallerFacts &Facts = It->second;
    if (New && Liveness::supports(F)) {
      Liveness Live(FAM.getResult<TargetIRAnalysis>(F), F);
      Facts.Complete = Live.compute(F);
      if (Facts.Complete)
        for (BasicBlock &B : F) {
          if (!any_of(B, isCallPressurePoint))
            continue;
          BlockPressure P(B, Live);
          for (unsigned I : P.Calls.set_bits())
            Facts.Through[P.Code[I]] = P.Through[I];
        }
    }
    const auto P = Facts.Through.find(&C);
    if (!Facts.Complete || P == Facts.Through.end())
      return std::nullopt;
    return P->second;
  }

  bool simdHeadroom(Function &F) {
    if (!Liveness::supports(F))
      return false;
    Liveness Live(FAM.getResult<TargetIRAnalysis>(F), F);
    return Live.compute(F) && summarize(F, Live).Peak.SIMD <= 6;
  }

  double tripsPerEntry(CallInst &C) {
    const auto &BFI = FAM.getResult<BlockFrequencyAnalysis>(*C.getFunction());
    const uint64_t Count = BFI.getBlockFreq(C.getParent()).getFrequency();
    const uint64_t Entry = BFI.getEntryFreq().getFrequency();
    return Entry ? std::min(1024.0, static_cast<double>(Count) / Entry) : 0.0;
  }
};

void reportIndirectScalarFields(Function &F, FunctionAnalysisManager &FAM,
                                Model &Costs) {
  if (!F.hasLocalLinkage() || F.isVarArg() || !Liveness::supports(F))
    return;
  auto &ORE = FAM.getResult<OptimizationRemarkEmitterAnalysis>(F);
  if (!ORE.allowExtraAnalysis(RemarkName))
    return;
  BitVector Sources(F.arg_size());
  for (const Argument &A : F.args()) {
    if (!A.getType()->isPointerTy() ||
        A.getType()->getPointerAddressSpace() != 0 || !A.onlyReadsMemory())
      continue;
    auto Attrs = F.getAttributes().getParamAttrs(A.getArgNo());
    if (capturesNothing(Attrs.getCaptureInfo()))
      Sources.set(A.getArgNo());
  }
  if (Sources.none())
    return;

  const auto &DL = F.getDataLayout();
  const auto &DT = FAM.getResult<DominatorTreeAnalysis>(F);
  SmallDenseSet<std::pair<unsigned, int64_t>, 16> Fields;
  BitVector Used(F.arg_size());
  unsigned Reads = 0;
  for (BasicBlock &B : F) {
    if (!DT.isReachableFromEntry(&B))
      continue;
    for (Instruction &I : B) {
      const auto *Read = dyn_cast<LoadInst>(&I);
      if (!Read || !Read->isSimple() || Read->use_empty() ||
          Read->getType()->isPointerTy() ||
          bank(Read->getType()) == Bank::Unknown)
        continue;
      const auto *Pointer = dyn_cast<LoadInst>(Read->getPointerOperand());
      if (!Pointer || !Pointer->isSimple() ||
          Pointer->getType()->getPointerAddressSpace() != 0)
        continue;
      int64_t Offset;
      const auto *A = dyn_cast<Argument>(
          GetPointerBaseWithConstantOffset(Pointer->getPointerOperand(), Offset,
                                           DL, /*AllowNonInbounds=*/false));
      if (!A || !Sources.test(A->getArgNo()))
        continue;
      Fields.insert({A->getArgNo(), Offset});
      Used.set(A->getArgNo());
      ++Reads;
    }
  }
  if (Fields.size() < 2)
    return;

  unsigned Sites = 0;
  unsigned CallerStack = 0;
  double Hot = 0;
  for (const Use &U : F.uses()) {
    auto *Call = dyn_cast<CallInst>(U.getUser());
    if (!Call || !Call->isCallee(&U) ||
        Call->getFunctionType() != F.getFunctionType())
      continue;
    ++Sites;
    Hot = std::max(Hot, Costs.tripsPerEntry(*Call));
    for (unsigned Arg : Used.set_bits()) {
      const auto *Root =
          dyn_cast<AllocaInst>(getUnderlyingObject(Call->getArgOperand(Arg)));
      if (Root && Root->isStaticAlloca() &&
          Root->getFunction() == Call->getFunction()) {
        ++CallerStack;
        break;
      }
    }
  }
  if (Sites == 0)
    return;

  ORE.emit([&] {
    using NV = ore::NV;
    return OptimizationRemarkAnalysis(RemarkName, "IndirectScalarFields", &F)
           << NV("Function", F.getName()) << " indirect scalar fields "
           << NV("Fields", Fields.size()) << ", scalar loads "
           << NV("Loads", Reads) << ", direct sites " << NV("Sites", Sites)
           << ", caller stack environments " << NV("StackSites", CallerStack)
           << ", peak site frequency per caller entry "
           << NV("PeakFrequency", static_cast<float>(Hot));
  });
}

int forwardedShuffleCost(const CallInst &C) {
  const Function &Caller = *C.getFunction();
  if (!originalConvention(Caller.getCallingConv()))
    return SysVGP; // No assumption about the original physical argument map.

  // A value forwarded from the identical SysV prefix slot loses a possible
  // zero-copy edge. IR cannot know the placement of other values.
  DenseMap<const Argument *, unsigned> EntrySlots;
  unsigned Slot = 0;
  for (const Argument &Arg : Caller.args()) {
    const Bank B = bank(Arg.getType());
    if (B == Bank::Unknown ||
        !signatureAttrs(Caller.getAttributes().getParamAttrs(Arg.getArgNo()),
                        Arg.getArgNo(), Caller.getReturnType()))
      return SysVGP;
    if (B == Bank::GP && Slot < SysVGP)
      EntrySlots[&Arg] = Slot++;
  }

  // Encoded as architectural GP numbers, in argument order. Only matching
  // first-six positions can count as forwards under the original ABI.
  constexpr unsigned Old[SysVGP] = {7, 6, 2, 1, 8, 9};
  constexpr unsigned New[SysVGP] = {0, 1, 2, 7, 6, 8};
  Slot = 0;
  int Copies = 0;
  for (const Use &U : C.args()) {
    if (bank(U->getType()) != Bank::GP || Slot >= SysVGP)
      continue;
    if (const auto *Arg = dyn_cast<Argument>(U.get())) {
      auto Entry = EntrySlots.find(Arg);
      if (Entry != EntrySlots.end() && Old[Entry->second] == Old[Slot] &&
          Old[Slot] != New[Slot])
        ++Copies;
    }
    ++Slot;
  }
  return Copies;
}

unsigned requiredStackReads(Function &F, const DominatorTree &DT) {
  SmallVector<const BasicBlock *, 8> Returns;
  for (const BasicBlock &B : F)
    if (isa<ReturnInst>(B.getTerminator()) && DT.isReachableFromEntry(&B))
      Returns.push_back(&B);
  if (Returns.empty())
    return 0;

  unsigned GP = 0, Reads = 0;
  for (const Argument &A : F.args()) {
    if (bank(A.getType()) != Bank::GP || GP++ < SysVGP)
      continue;
    Reads += any_of(A.users(), [&](const User *U) {
      const auto *I = dyn_cast<Instruction>(U);
      return I && !isa<PHINode, IntrinsicInst>(I) &&
             all_of(Returns, [&](const BasicBlock *B) {
               return DT.dominates(I->getParent(), B);
             });
    });
  }
  return Reads;
}

struct Traffic {
  int ExtraCSR;
  int MaxLiveGP = 0;
  int WorstSite = std::numeric_limits<int>::max();
  double MinGain = 0;
};

std::optional<Traffic> estimate(const Layout &L, unsigned RequiredReads,
                                ArrayRef<CallInst *> Sites, Model &M) {
  const int Registers =
      std::max(0, static_cast<int>(L.GP) - static_cast<int>(RegCallVolatileGP));
  Traffic T{Registers};
  DenseMap<const Function *, double> Gains;
  for (CallInst *C : Sites) {
    auto Live = M.through(*C);
    if (!Live)
      return std::nullopt;
    T.MaxLiveGP = std::max(T.MaxLiveGP, Live->GP);
    const int Lost =
        std::clamp(Live->GP - (PreservedGPBudget - Registers), 0, Registers);
    const int PerCall = StackArgumentAccess *
                            (static_cast<int>(L.UsedStackGP + RequiredReads)) -
                        CrossCallSpillPair * Lost - forwardedShuffleCost(*C);
    const double Weight = M.tripsPerEntry(*C);
    T.WorstSite = std::min(T.WorstSite, PerCall);
    Gains[C->getFunction()] += Weight * PerCall;
  }
  T.MinGain = std::numeric_limits<double>::max();
  for (const auto &Entry : Gains)
    T.MinGain = std::min(T.MinGain, Entry.second - FrameSavePair * Registers);
  return T;
}

} // namespace

PreservedAnalyses CallABIPass::run(Module &M, ModuleAnalysisManager &AM) {
  const Triple &T = M.getTargetTriple();
  if (T.getArch() != Triple::x86_64 || (!T.isOSLinux() && !T.isOSWindows()) ||
      M.getDataLayout().getPointerSizeInBits(0) != 64 ||
      !M.getModuleInlineAsm().empty())
    return PreservedAnalyses::all();

  auto &FAM = AM.getResult<FunctionAnalysisManagerModuleProxy>(M).getManager();
  Model Costs(FAM);
  if (!Transform)
    for (Function &F : M)
      if (!F.isDeclaration() && !F.hasOptNone() && F.getName().contains(Filter))
        reportIndirectScalarFields(F, FAM, Costs);
  if (!T.isOSLinux())
    return PreservedAnalyses::all();

  SmallVector<StringRef, 16> Assembly;
  for (Function &F : M)
    for (BasicBlock &B : F)
      for (Instruction &I : B)
        if (const auto *C = dyn_cast<CallBase>(&I))
          if (const auto *Asm = dyn_cast<InlineAsm>(C->getCalledOperand()))
            if (!Asm->getAsmString().empty())
              Assembly.push_back(Asm->getAsmString());

  unsigned Applied = 0;
  for (Function &F : M) {
    if (!F.getName().contains(Filter) || F.isDeclaration() ||
        !F.hasLocalLinkage() || F.hasComdat() ||
        !originalConvention(F.getCallingConv()) || F.arg_size() <= SysVGP ||
        F.hasOptNone() || F.hasOptSize() || F.hasMinSize() || F.hasGC() ||
        F.hasFnAttribute(Attribute::Naked) ||
        F.hasFnAttribute(Attribute::ReturnsTwice))
      continue;
    const Layout L = layout(F);
    if (!L.Supported || L.GP <= SysVGP || L.GP > RegCallGP || L.SSE > 8 ||
        !L.UsedStackGP || !canUseSSE(F))
      continue;

    const char *Reason = "eligible";
    SmallVector<CallInst *, 4> Sites;
    std::optional<Traffic> Traffic;
    unsigned RequiredReads = 0;
    if (!collectSites(F, Sites, Assembly)) {
      Reason = "nonlocal, indirect, unwind, or tail boundary";
    } else if (!F.doesNotThrow() || !ordinaryLeaf(F)) {
      Reason = "nonleaf callee";
    } else if (!Costs.simdHeadroom(F)) {
      Reason = "SIMD save risk or unmodeled callee";
    } else {
      RequiredReads =
          requiredStackReads(F, FAM.getResult<DominatorTreeAnalysis>(F));
      Traffic = estimate(L, RequiredReads, Sites, Costs);
      if (!Traffic)
        Reason = "caller liveness unavailable";
      else if (Traffic->WorstSite <= 0 || Traffic->MinGain < 2.0)
        Reason = "call shuffles or save traffic exceed benefit";
      else if (Transform && Applied < 96) {
        F.setCallingConv(CallingConv::X86_RegCall);
        for (CallInst *C : Sites)
          C->setCallingConv(CallingConv::X86_RegCall);
        ++Applied;
        Reason = "applied x86_regcallcc";
      } else if (Transform) {
        Reason = "module rewrite budget";
      }
    }

    auto &ORE = FAM.getResult<OptimizationRemarkEmitterAnalysis>(F);
    ORE.emit([&] {
      using NV = ore::NV;
      return OptimizationRemarkAnalysis(RemarkName, "StackArguments", &F)
             << NV("Function", F.getName()) << " used stack GP "
             << NV("StackGP", L.UsedStackGP) << ", required reads "
             << NV("RequiredReads", RequiredReads) << ", direct sites "
             << NV("Sites", Sites.size()) << ", result " << NV("Result", Reason)
             << ", estimated extra CSR "
             << NV("ExtraCSR", Traffic ? Traffic->ExtraCSR : 0)
             << ", max live GP "
             << NV("LiveGP", Traffic ? Traffic->MaxLiveGP : 0)
             << ", worst direct site "
             << NV("WorstSite", Traffic ? Traffic->WorstSite : 0)
             << ", min caller score "
             << NV("Score",
                   Traffic ? static_cast<float>(Traffic->MinGain) : 0.0f);
    });
  }
  return Applied ? PreservedAnalyses::none() : PreservedAnalyses::all();
}

} // namespace rssn
