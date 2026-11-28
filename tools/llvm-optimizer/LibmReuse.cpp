#include "LibmReuse.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/OptimizationRemarkEmitter.h"
#include "llvm/Analysis/PostDominators.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/MDBuilder.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/PatternMatch.h"
#include "llvm/Transforms/Utils/BasicBlockUtils.h"
#include "llvm/Transforms/Utils/Local.h"
#include <optional>

using namespace llvm;
using namespace llvm::PatternMatch;

namespace rssn {
namespace {

constexpr char PassName[] = "rssn-libm-reuse";

struct ExpInput {
  Value *Seed;
  const ConstantFP *BeforeAbs;
  const ConstantFP *AfterAbs;
  bool Negated;
};

// Limit the profitability heuristic to sign-symmetric input expressions.
// Correctness relies on checking exp's actual operands at runtime.
std::optional<ExpInput> matchInput(Value *V) {
  ConstantFP *Outer = nullptr, *Inner = nullptr;
  Value *Operand;
  if (match(V, m_c_FMul(m_Value(Operand), m_ConstantFP(Outer))))
    V = Operand;
  if (!match(V, m_Intrinsic<Intrinsic::fabs>(m_Value(Operand))))
    return std::nullopt;
  if (match(Operand, m_c_FMul(m_Value(V), m_ConstantFP(Inner))))
    Operand = V;
  auto *Call = dyn_cast<CallInst>(Operand);
  if (!Call || Call->arg_size() != 1 || !Call->getCalledFunction() ||
      !Call->getArgOperand(0)->getType()->isDoubleTy())
    return std::nullopt;
  auto *II = dyn_cast<IntrinsicInst>(Call);
  if (II ? II->getIntrinsicID() != Intrinsic::sinh
         : (Call->isNoBuiltin() || !Call->getCalledFunction()->isDeclaration() ||
            Call->getCalledFunction()->getName() != "sinh"))
    return std::nullopt;

  Value *Seed = Call->getArgOperand(0);
  if (match(Seed, m_FNeg(m_Value(V))))
    return ExpInput{V, Inner, Outer, true};
  return ExpInput{Seed, Inner, Outer, false};
}

struct Candidate {
  IntrinsicInst *Exp;
  const Loop *LoopRegion;
  ExpInput Input;
};

bool exactFP(const Instruction &I) {
  return isa<FPMathOperator>(&I) && I.getFastMathFlags().none() &&
         !I.getMetadata(LLVMContext::MD_fpmath) &&
         !I.getMetadata(LLVMContext::MD_nofpclass);
}

bool isExp(const IntrinsicInst &I) {
  return I.getIntrinsicID() == Intrinsic::exp && I.getType()->isDoubleTy() &&
         !I.isMustTailCall() && !I.hasOperandBundles() && !I.isConvergent() &&
         !I.hasFnAttr(Attribute::StrictFP) && I.doesNotAccessMemory() &&
         !I.hasRetAttr(Attribute::NoFPClass) &&
         !I.paramHasAttr(0, Attribute::NoFPClass) && exactFP(I);
}

struct ReusePlan {
  IntrinsicInst *First;
  Instruction *Available;
  SmallVector<Instruction *, 4> Recompute;
};

BinaryOperator *constantConsumer(Instruction &I, const ConstantFP *&Other,
                                 unsigned &InputIndex) {
  if (!I.hasOneUse())
    return nullptr;
  auto *Op = dyn_cast<BinaryOperator>(*I.user_begin());
  if (!Op)
    return nullptr;
  switch (Op->getOpcode()) {
  case Instruction::FAdd:
  case Instruction::FSub:
  case Instruction::FMul:
  case Instruction::FDiv:
    break;
  default:
    return nullptr;
  }
  InputIndex = Op->getOperand(0) == &I ? 0 : 1;
  Other = dyn_cast<ConstantFP>(Op->getOperand(1 - InputIndex));
  return Other && exactFP(*Op) ? Op : nullptr;
}

ReusePlan planReuse(IntrinsicInst &First, IntrinsicInst &Second,
                    const DominatorTree &DT) {
  ReusePlan Plan{&First, &First, {&Second}};
  while (Plan.Recompute.size() < 4) {
    const ConstantFP *FirstConstant, *SecondConstant;
    unsigned FirstIndex, SecondIndex;
    auto *A = constantConsumer(*Plan.Available, FirstConstant, FirstIndex);
    auto *B = constantConsumer(*Plan.Recompute.back(), SecondConstant, SecondIndex);
    if (!A || !B || A->getOpcode() != B->getOpcode() ||
        FirstConstant != SecondConstant || FirstIndex != SecondIndex ||
        B->getParent() != Second.getParent() || !DT.dominates(A, &Second))
      break;
    Plan.Available = A;
    Plan.Recompute.push_back(B);
  }
  return Plan;
}

void reuseExp(const ReusePlan &Plan) {
  IntrinsicInst &First = *Plan.First;
  auto &Second = *cast<IntrinsicInst>(Plan.Recompute.front());
  IRBuilder<> Capture(&First);
  Value *Saved = Capture.CreateFreeze(First.getArgOperand(0), "exp.previous");
  First.setArgOperand(0, Saved);
  Value *SavedBits = Capture.CreateBitCast(Saved, Capture.getInt64Ty());

  IRBuilder<> Check(&Second);
  Value *Current = Check.CreateFreeze(Second.getArgOperand(0), "exp.current");
  Second.setArgOperand(0, Current);
  Value *CurrentBits = Check.CreateBitCast(Current, Check.getInt64Ty());
  Value *Different = Check.CreateICmpNE(SavedBits, CurrentBits, "exp.differs");

  BasicBlock *Head = Second.getParent();
  MDNode *Weights =
      MDBuilder(Second.getContext()).createBranchWeights(1, 15, /*IsExpected=*/true);
  Instruction *Jump = SplitBlockAndInsertIfThen(Different, &Second, false, Weights);
  BasicBlock *Evaluate = Jump->getParent();
  Evaluate->setName("exp.recompute");
  Instruction *Last = Plan.Recompute.back();
  auto *Result = PHINode::Create(Last->getType(), 2, "exp.reuse", Second.getIterator());
  Result->setDebugLoc(Last->getDebugLoc());
  Last->replaceAllUsesWith(Result);
  for (Instruction *I : Plan.Recompute) {
    salvageDebugInfo(*I);
    I->moveBefore(Jump->getIterator());
  }
  Result->addIncoming(Plan.Available, Head);
  Result->addIncoming(Last, Evaluate);
}

} // namespace

PreservedAnalyses LibmReusePass::run(Function &F, FunctionAnalysisManager &AM) {
  const Triple &T = F.getParent()->getTargetTriple();
  if (!F.getName().contains(Filter) || T.getArch() != Triple::x86_64 ||
      (!T.isWindowsMSVCEnvironment() && (!T.isOSLinux() || T.isX32())) ||
      F.isDeclaration() || F.hasOptNone() || F.hasOptSize() || F.hasMinSize() ||
      F.hasFnAttribute(Attribute::StrictFP) || F.hasGC() ||
      F.size() > 1024 || F.getInstructionCount() > 24000)
    return PreservedAnalyses::all();
  const auto &LI = AM.getResult<LoopAnalysis>(F);
  if (LI.empty())
    return PreservedAnalyses::all();

  DenseMap<Value *, SmallVector<Candidate, 2>> Candidates;
  SmallVector<Value *, 8> Seeds;
  for (BasicBlock &B : F) {
    const auto *LoopRegion = LI.getLoopFor(&B);
    if (!LoopRegion || B.isEHPad())
      continue;
    for (Instruction &I : B) {
      auto *II = dyn_cast<IntrinsicInst>(&I);
      if (!II || !isExp(*II))
        continue;
      if (auto Input = matchInput(II->getArgOperand(0))) {
        auto [It, New] = Candidates.try_emplace(Input->Seed);
        if (New)
          Seeds.push_back(Input->Seed);
        It->second.push_back({II, LoopRegion, *Input});
      }
    }
  }

  if (Seeds.empty())
    return PreservedAnalyses::all();
  const auto &DT = AM.getResult<DominatorTreeAnalysis>(F);
  const auto &PDT = AM.getResult<PostDominatorTreeAnalysis>(F);
  auto &ORE = AM.getResult<OptimizationRemarkEmitterAnalysis>(F);
  SmallVector<ReusePlan, 4> Plans;
  // In a single iteration require exactly two samples of the same seed.
  // Additional partners can make both the comparison and result long-lived.
  for (Value *Seed : Seeds) {
    const auto &Pair = Candidates.find(Seed)->second;
    if (Pair.size() != 2 || Pair[0].LoopRegion != Pair[1].LoopRegion ||
        Pair[0].Input.Negated == Pair[1].Input.Negated ||
        Pair[0].Input.BeforeAbs != Pair[1].Input.BeforeAbs ||
        Pair[0].Input.AfterAbs != Pair[1].Input.AfterAbs)
      continue;
    auto *First = Pair[0].Exp;
    auto *Second = Pair[1].Exp;
    if (DT.dominates(Second, First))
      std::swap(First, Second);
    if (!DT.dominates(First, Second) ||
        !PDT.dominates(Second->getParent(), First->getParent()))
      continue;

    auto Plan = planReuse(*First, *Second, DT);
    // The original argument and result can span several ABI clobbers.
    // Require an expensive consumer as well as the saved call.
    const bool Profitable = any_of(Plan.Recompute, [](const Instruction *I) {
      return I->getOpcode() == Instruction::FDiv;
    });
    const StringRef Action =
        !Profitable ? "kept" : (Transform ? "guarded" : "candidate");
    ORE.emit([&] {
      using NV = ore::NV;
      return OptimizationRemarkAnalysis(PassName, "MirroredExp", Second)
             << NV("Function", F.getName()) << " "
             << NV("Action", Action)
             << " reuse of scalar exp on the opposite-sign path with "
             << NV("ConsumerSteps", Plan.Recompute.size() - 1)
             << " matched consumers";
    });
    if (!Profitable)
      continue;
    Plans.push_back(std::move(Plan));
    if (Plans.size() == 4)
      break;
  }
  if (!Transform || Plans.empty())
    return PreservedAnalyses::all();
  for (const ReusePlan &Plan : Plans)
    reuseExp(Plan);
  return PreservedAnalyses::none();
}

} // namespace rssn
