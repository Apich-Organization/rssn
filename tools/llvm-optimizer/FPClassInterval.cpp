#include "FPClassInterval.h"
#include "llvm/ADT/FloatingPointMode.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/Analysis/AssumptionCache.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/OptimizationRemarkEmitter.h"
#include "llvm/Analysis/TargetLibraryInfo.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/PatternMatch.h"
#include "llvm/Support/KnownFPClass.h"
#include "llvm/Transforms/Utils/Local.h"
#include <optional>

using namespace llvm;
using namespace llvm::PatternMatch;

namespace rssn {
namespace {
constexpr char PassName[] = "rssn-fp-class-interval";
constexpr uint64_t FiniteNonzeroEnd = UINT64_C(0xffe0000000000000) - 2;
constexpr FPClassTest FiniteNonzero = fcNormal | fcSubnormal;
constexpr FPClassTest NonfiniteOrZero = fcNan | fcInf | fcZero;
constexpr unsigned MaxSites = 16;
constexpr unsigned MaxClassNodes = 24;

bool excludesUnconditionalClasses(Value *Input, FPClassTest Excluded,
                                  const DataLayout &DL,
                                  const TargetLibraryInfo &TLI) {
  SmallVector<Value *, 8> Pending{Input};
  SmallPtrSet<Value *, 16> Seen;
  // A loop-carried phi can retain any preceding iteration's value. Follow
  // copies to their producers without using facts local to this iteration.
  while (!Pending.empty()) {
    Value *V = Pending.pop_back_val();
    if (!Seen.insert(V).second)
      continue;
    if (Seen.size() > MaxClassNodes)
      return false;
    if (auto *Phi = dyn_cast<PHINode>(V)) {
      for (Value *Incoming : Phi->incoming_values())
        Pending.push_back(Incoming);
    } else if (auto *Select = dyn_cast<SelectInst>(V)) {
      Pending.push_back(Select->getTrueValue());
      Pending.push_back(Select->getFalseValue());
    } else if (!computeKnownFPClass(V, DL, Excluded, &TLI)
                    .isKnownNever(Excluded)) {
      return false;
    }
  }
  return true;
}

struct ClassTest {
  IntrinsicInst *Intrinsic;
  bool Invalid;
};

std::optional<ClassTest> loopClassTest(CondBrInst &Branch, const Loop &L,
                                       const DominatorTree &DT,
                                       AssumptionCache &AC,
                                       const TargetLibraryInfo &TLI) {
  auto *Test = dyn_cast<IntrinsicInst>(Branch.getCondition());
  if (!Test || Test->getIntrinsicID() != Intrinsic::is_fpclass ||
      !Test->hasOneUse() || Test->hasOperandBundles() || Test->isConvergent() ||
      Test->hasFnAttr(Attribute::StrictFP) ||
      Test->paramHasAttr(0, Attribute::NoFPClass) ||
      Test->getMetadata(LLVMContext::MD_nofpclass))
    return std::nullopt;
  const auto Mask = static_cast<FPClassTest>(
      cast<ConstantInt>(Test->getArgOperand(1))->getZExtValue());
  Value *Input = Test->getArgOperand(0);
  if (!Input->getType()->isDoubleTy() || L.isLoopInvariant(Input))
    return std::nullopt;
  const FPClassTest Symmetric = unknown_sign(Mask);
  if (Symmetric != FiniteNonzero && Symmetric != NonfiniteOrZero)
    return std::nullopt;
  const FPClassTest Added = Symmetric & ~Mask;
  // InstCombine may have trimmed the original mask using a sign proof.
  // Recover the cheaper symmetric test only if the added classes are absent.
  const DataLayout &DL = Test->getDataLayout();
  if (Added &&
      !computeKnownFPClass(Input, DL, Added, &TLI, &AC, Test, &DT)
           .isKnownNever(Added) &&
      !excludesUnconditionalClasses(Input, Added, DL, TLI))
    return std::nullopt;
  return ClassTest{Test, Symmetric == NonfiniteOrZero};
}
} // namespace

Value *buildFiniteNonzeroTest(IRBuilder<> &B, Value *Weight, bool Invalid) {
  Value *Bits = B.CreateBitCast(Weight, B.getInt64Ty(), "weight.bits");
  // Doubling discards the sign. Subtracting two puts the smallest subnormal
  // at the start of the unsigned interval and moves both zeros past infinity.
  Value *Magnitude = B.CreateShl(Bits, 1);
  Value *Span = B.CreateSub(Magnitude, B.getInt64(2));
  return B.CreateICmp(Invalid ? ICmpInst::ICMP_UGE : ICmpInst::ICMP_ULT, Span,
                      B.getInt64(FiniteNonzeroEnd),
                      Invalid ? "weight.zero.or.nonfinite"
                              : "weight.nonzero.finite");
}

ICmpInst *findNonfiniteOrZeroTest(BitCastInst &Bits) {
  if (!Bits.getSrcTy()->isDoubleTy() || !Bits.getDestTy()->isIntegerTy(64) ||
      !Bits.hasOneUse())
    return nullptr;
  auto *Magnitude = dyn_cast<Instruction>(*Bits.user_begin());
  if (!Magnitude || !Magnitude->hasOneUse() ||
      !match(Magnitude, m_Shl(m_Specific(&Bits), m_SpecificInt(1))))
    return nullptr;
  auto *Span = dyn_cast<Instruction>(*Magnitude->user_begin());
  if (!Span || !Span->hasOneUse() ||
      !match(Span, m_Sub(m_Specific(Magnitude), m_SpecificInt(2))))
    return nullptr;
  auto *Cmp = dyn_cast<ICmpInst>(*Span->user_begin());
  if (!Cmp || Cmp->getPredicate() != ICmpInst::ICMP_UGE ||
      Cmp->getOperand(0) != Span ||
      !match(Cmp->getOperand(1), m_SpecificInt(FiniteNonzeroEnd)))
    return nullptr;
  return Cmp;
}

PreservedAnalyses FPClassIntervalPass::run(Function &F,
                                           FunctionAnalysisManager &AM) {
  const Triple &T = F.getParent()->getTargetTriple();
  if (T.getArch() != Triple::x86_64 ||
      (!T.isWindowsMSVCEnvironment() && (!T.isOSLinux() || T.isX32())) ||
      !F.getName().contains(Filter) || F.isDeclaration() || F.hasOptNone() ||
      F.hasOptSize() || F.hasMinSize() ||
      F.hasFnAttribute(Attribute::StrictFP) ||
      F.getDenormalMode(APFloat::IEEEdouble()) != DenormalMode::getIEEE() ||
      F.size() > 1024 || F.getInstructionCount() > 24000)
    return PreservedAnalyses::all();

  const auto &LI = AM.getResult<LoopAnalysis>(F);
  if (LI.empty())
    return PreservedAnalyses::all();
  const auto &DT = AM.getResult<DominatorTreeAnalysis>(F);
  auto &AC = AM.getResult<AssumptionAnalysis>(F);
  const auto &TLI = AM.getResult<TargetLibraryAnalysis>(F);
  auto &ORE = AM.getResult<OptimizationRemarkEmitterAnalysis>(F);
  unsigned Applied = 0;
  for (BasicBlock &B : F) {
    const auto *L = LI.getLoopFor(&B);
    auto *Branch = dyn_cast<CondBrInst>(B.getTerminator());
    if (!L || !Branch || B.isEHPad())
      continue;
    if (auto Class = loopClassTest(*Branch, *L, DT, AC, TLI)) {
      IntrinsicInst *Test = Class->Intrinsic;
      ORE.emit([&] {
        using NV = ore::NV;
        return OptimizationRemarkAnalysis(PassName, "FiniteNonzeroClass", Test)
               << NV("Function", F.getName()) << " "
               << NV("Action", Transform ? "compacted" : "can compact")
               << " a finite nonzero class check into one unsigned interval, "
                  "mask "
               << NV("Mask",
                     cast<ConstantInt>(Test->getArgOperand(1))->getZExtValue())
               << ", loop depth " << NV("LoopDepth", L->getLoopDepth());
      });
      if (Transform) {
        IRBuilder<> Builder(Test);
        Branch->setCondition(buildFiniteNonzeroTest(
            Builder, Test->getArgOperand(0), Class->Invalid));
        RecursivelyDeleteTriviallyDeadInstructions(Test);
        if (++Applied == MaxSites)
          break;
      }
    }
  }
  if (!Applied)
    return PreservedAnalyses::all();
  PreservedAnalyses PA;
  PA.preserveSet<CFGAnalyses>();
  return PA;
}
} // namespace rssn
