#include "FPEndpointGuards.h"
#include "FPClassInterval.h"
#include "llvm/ADT/FloatingPointMode.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/OptimizationRemarkEmitter.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/PatternMatch.h"
#include "llvm/Transforms/Utils/Local.h"
#include <optional>

using namespace llvm;
using namespace llvm::PatternMatch;

namespace rssn {
namespace {

constexpr char PassName[] = "rssn-fp-endpoint-guards";

struct Test {
  Value *Other;
  Value *Input;
};

bool exactFP(const Instruction &I) {
  return isa<FPMathOperator>(&I) && I.getFastMathFlags().none() &&
         !I.getMetadata(LLVMContext::MD_fpmath) &&
         !I.getMetadata(LLVMContext::MD_nofpclass);
}

Value *compareAgainst(Value *V, FCmpInst::Predicate Pred,
                      const APFloat &Bound) {
  const auto *Cmp = dyn_cast<FCmpInst>(V);
  if (!Cmp || !exactFP(*Cmp))
    return nullptr;
  if (Cmp->getPredicate() == Pred)
    if (const auto *C = dyn_cast<ConstantFP>(Cmp->getOperand(1));
        C && C->getValueAPF().bitwiseIsEqual(Bound))
      return Cmp->getOperand(0);
  if (Cmp->getSwappedPredicate() == Pred)
    if (const auto *C = dyn_cast<ConstantFP>(Cmp->getOperand(0));
        C && C->getValueAPF().bitwiseIsEqual(Bound))
      return Cmp->getOperand(1);
  return nullptr;
}

Value *finiteInput(Value *V) {
  if (auto *Class = dyn_cast<IntrinsicInst>(V);
      Class && Class->getIntrinsicID() == Intrinsic::is_fpclass)
    return cast<ConstantInt>(Class->getArgOperand(1))->getZExtValue() ==
                   fcFinite
               ? Class->getArgOperand(0)
               : nullptr;
  const APFloat Infinity = APFloat::getInf(APFloat::IEEEdouble());
  Value *Input = nullptr;
  auto *Abs = dyn_cast_or_null<IntrinsicInst>(
      compareAgainst(V, FCmpInst::FCMP_ONE, Infinity));
  if (Abs && exactFP(*Abs) &&
      match(Abs, m_Intrinsic<Intrinsic::fabs>(m_Value(Input))))
    return Input;
  return nullptr;
}

std::optional<Test> findTest(Value *V, unsigned Opcode,
                             function_ref<Value *(Value *)> Decode) {
  auto *Bin = dyn_cast<BinaryOperator>(V);
  if (!Bin || Bin->getOpcode() != Opcode)
    return std::nullopt;
  for (unsigned I : {0u, 1u})
    if (Value *Input = Decode(Bin->getOperand(I)))
      return Test{Bin->getOperand(1 - I), Input};
  return std::nullopt;
}

bool speculativeCorridor(const BasicBlock &B) {
  return !B.isEHPad() && all_of(B, [](const Instruction &I) {
    return I.isTerminator() || isa<DbgInfoIntrinsic>(I) ||
           (!I.mayReadOrWriteMemory() && isSafeToSpeculativelyExecute(&I));
  });
}

bool positiveFinite(const ConstantFP *C) {
  if (!C)
    return false;
  const APFloat &V = C->getValueAPF();
  return V.isFinite() && !V.isZero() && !V.isNegative();
}

bool sameFailure(const BasicBlock &Failure, const BasicBlock &A,
                 const BasicBlock &B) {
  return all_of(Failure.phis(), [&](const PHINode &Phi) {
    return Phi.getIncomingValueForBlock(&A) == Phi.getIncomingValueForBlock(&B);
  });
}

struct Plan {
  CondBrInst *First;
  CondBrInst *Last;
  SelectInst *Node;
  Value *Weight;
  Value *UpperBoundTest;
  Instruction *Positive;
  SmallVector<Instruction *, 3> Sink;
};

std::optional<Plan> planGuards(CondBrInst &First, const Loop &L) {
  const APFloat Zero = APFloat::getZero(APFloat::IEEEdouble());
  auto Finite = findTest(First.getCondition(), Instruction::And, finiteInput);
  if (!Finite || !Finite->Input->getType()->isDoubleTy())
    return std::nullopt;
  // The original test allows NaNs to reach the endpoint tests.
  Value *Distance = compareAgainst(Finite->Other, FCmpInst::FCMP_UGT, Zero);
  if (!Distance)
    return std::nullopt;
  auto *Positive = cast<FCmpInst>(Finite->Other);

  BasicBlock &Failure = *First.getSuccessor(1);
  BasicBlock &Lower = *First.getSuccessor(0);
  if (&Failure == &Lower || Lower.getSinglePredecessor() != First.getParent() ||
      !L.contains(&Lower) || !speculativeCorridor(Lower))
    return std::nullopt;
  auto *LowerGuard = dyn_cast<CondBrInst>(Lower.getTerminator());
  if (!LowerGuard || LowerGuard->getSuccessor(1) != &Failure)
    return std::nullopt;
  auto *Node = dyn_cast_or_null<SelectInst>(
      compareAgainst(LowerGuard->getCondition(), FCmpInst::FCMP_UGT, Zero));
  if (!Node || !exactFP(*Node) || Node->getParent() != &Lower)
    return std::nullopt;

  BasicBlock &Upper = *LowerGuard->getSuccessor(0);
  if (Upper.getSinglePredecessor() != &Lower || !L.contains(&Upper) ||
      !speculativeCorridor(Upper))
    return std::nullopt;
  auto *Last = dyn_cast<CondBrInst>(Upper.getTerminator());
  if (!Last || Last->getSuccessor(0) != &Failure ||
      Last->getSuccessor(1) == &Failure ||
      !sameFailure(Failure, *First.getParent(), Lower) ||
      !sameFailure(Failure, *First.getParent(), Upper))
    return std::nullopt;
  auto IsZero = findTest(Last->getCondition(), Instruction::Or, [&](Value *V) {
    return compareAgainst(V, FCmpInst::FCMP_OEQ, Zero);
  });
  if (!IsZero || IsZero->Input != Finite->Input)
    return std::nullopt;

  // Both arms use the same scaled distance. A nonpositive distance
  // reaches either endpoint (or extends past it), so the existing lower or
  // upper bound test already takes exactly the same failure edge.
  Value *Near = Node->getTrueValue(), *Far = Node->getFalseValue();
  auto *Subtract = dyn_cast<BinaryOperator>(Far);
  if (!Subtract || Subtract->getOpcode() != Instruction::FSub) {
    std::swap(Near, Far);
    Subtract = dyn_cast<BinaryOperator>(Far);
  }
  auto *Scale = dyn_cast<BinaryOperator>(Near);
  ConstantFP *Factor = nullptr;
  auto *Bound =
      Subtract ? dyn_cast<ConstantFP>(Subtract->getOperand(0)) : nullptr;
  if (!Subtract || Subtract->getOpcode() != Instruction::FSub ||
      Subtract->getOperand(1) != Near || !exactFP(*Subtract) || !Scale ||
      !exactFP(*Scale) ||
      !match(Scale, m_c_FMul(m_Specific(Distance), m_ConstantFP(Factor))) ||
      !positiveFinite(Factor) || !positiveFinite(Bound) ||
      compareAgainst(IsZero->Other, FCmpInst::FCMP_OGE, Bound->getValueAPF()) !=
          Node)
    return std::nullopt;

  return Plan{&First, Last, Node, Finite->Input, IsZero->Other, Positive, {}};
}

void planExpSink(Plan &P, const DominatorTree &DT) {
  if (!P.Positive->hasOneUse() || !P.First->getCondition()->hasOneUse())
    return;
  auto *Div = cast<FCmpInst>(P.Positive)->getOperand(0);
  if (isa<Constant>(Div))
    Div = cast<FCmpInst>(P.Positive)->getOperand(1);
  auto *Ratio = dyn_cast<BinaryOperator>(Div);
  auto *Add = Ratio ? dyn_cast<BinaryOperator>(Ratio->getOperand(1)) : nullptr;
  if (!Ratio || Ratio->getOpcode() != Instruction::FDiv ||
      !isa<ConstantFP>(Ratio->getOperand(0)) || !Add || !Add->hasOneUse() ||
      !match(Add,
             m_c_FAdd(m_Intrinsic<Intrinsic::exp>(m_Value()), m_ConstantFP())))
    return;
  auto *Exp = dyn_cast<IntrinsicInst>(Add->getOperand(0));
  if (!Exp)
    Exp = cast<IntrinsicInst>(Add->getOperand(1));
  if (!Exp->hasOneUse() || Exp->isConvergent() || Exp->hasOperandBundles() ||
      Exp->isMustTailCall() || Exp->hasFnAttr(Attribute::StrictFP) ||
      !Exp->doesNotAccessMemory() || !Exp->hasFnAttr(Attribute::Speculatable))
    return;
  for (Instruction *I : {cast<Instruction>(Exp), cast<Instruction>(Add),
                         cast<Instruction>(Ratio)})
    if (I->getParent() != P.First->getParent() || !exactFP(*I))
      return;

  Instruction &Insert = *P.First->getSuccessor(0)->getFirstInsertionPt();
  for (const Use &U : Ratio->uses())
    if (U.getUser() != P.Positive && U.getUser() != &Insert &&
        !DT.dominates(&Insert, U))
      return;
  P.Sink = {Exp, Add, Ratio};
}

void apply(const Plan &P) {
  IRBuilder<> B(P.First);
  Value *Valid = buildFiniteNonzeroTest(B, P.Weight, false);
  IRBuilder<> Direction(P.Node);
  // A previously skipped node can now reach this select. Both endpoint
  // choices fail for a nonpositive distance, even when the selector is poison.
  P.Node->setCondition(
      Direction.CreateFreeze(P.Node->getCondition(), "endpoint.direction"));
  auto *OldFirst = cast<Instruction>(P.First->getCondition());
  auto *OldLast = cast<Instruction>(P.Last->getCondition());
  P.First->setCondition(Valid);
  P.Last->setCondition(P.UpperBoundTest);
  // The conditions and the populations reaching each test have changed.
  P.First->setMetadata(LLVMContext::MD_prof, nullptr);
  P.First->getSuccessor(0)->getTerminator()->setMetadata(LLVMContext::MD_prof,
                                                         nullptr);
  P.Last->setMetadata(LLVMContext::MD_prof, nullptr);
  RecursivelyDeleteTriviallyDeadInstructions(OldFirst);
  RecursivelyDeleteTriviallyDeadInstructions(OldLast);
  auto Insert = P.First->getSuccessor(0)->getFirstInsertionPt();
  for (Instruction *I : P.Sink) {
    salvageDebugInfo(*I);
    I->moveBefore(Insert);
  }
}

} // namespace

PreservedAnalyses FPEndpointGuardsPass::run(Function &F,
                                            FunctionAnalysisManager &AM) {
  if (!F.getName().contains(Filter) || F.isDeclaration() || F.hasOptNone() ||
      F.hasOptSize() || F.hasMinSize() ||
      F.hasFnAttribute(Attribute::StrictFP) ||
      F.getDenormalMode(APFloat::IEEEdouble()) != DenormalMode::getIEEE() ||
      F.size() > 1024 || F.getInstructionCount() > 24000)
    return PreservedAnalyses::all();

  const auto &LI = AM.getResult<LoopAnalysis>(F);
  if (LI.empty())
    return PreservedAnalyses::all();
  const auto &DT = AM.getResult<DominatorTreeAnalysis>(F);
  auto &ORE = AM.getResult<OptimizationRemarkEmitterAnalysis>(F);
  unsigned Applied = 0;
  for (BasicBlock &B : F) {
    const auto *L = LI.getLoopFor(&B);
    auto *Branch = dyn_cast<CondBrInst>(B.getTerminator());
    if (!L || !Branch || B.isEHPad())
      continue;
    auto P = planGuards(*Branch, *L);
    if (!P)
      continue;
    planExpSink(*P, DT);
    ORE.emit([&] {
      using NV = ore::NV;
      return OptimizationRemarkAnalysis(PassName, "EndpointGuard", Branch)
             << NV("Function", F.getName()) << " "
             << NV("Action", Transform ? "merged" : "can merge")
             << " endpoint guards, " << NV("Sunk", P->Sink.size())
             << " exp chain operations guarded by finite nonzero weight";
    });
    if (Transform) {
      apply(*P);
      if (++Applied == 16)
        break;
    }
  }
  if (!Applied)
    return PreservedAnalyses::all();
  PreservedAnalyses PA;
  PA.preserveSet<CFGAnalyses>();
  return PA;
}
} // namespace rssn
