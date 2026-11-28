#include "LibmPressure.h"
#include "Pressure.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/OptimizationRemarkEmitter.h"
#include "llvm/Analysis/TargetLibraryInfo.h"
#include "llvm/Analysis/TargetTransformInfo.h"
#include "llvm/Analysis/VectorUtils.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/ProfDataUtils.h"
#include "llvm/Transforms/Utils/BasicBlockUtils.h"
#include "llvm/Transforms/Utils/Local.h"

#include <optional>

using namespace llvm;

namespace rssn {

// Identify scalar double operations that can lower to a libm call on the
// supported targets. Constant-exponent combines and vector library selection
// have different lowering rules.
bool isScalarLibmCut(const Instruction &I) {
  const auto *Call = dyn_cast<CallInst>(&I);
  if (!Call || !I.getType()->isDoubleTy() || Call->isMustTailCall() ||
      Call->isConvergent() || Call->hasOperandBundles() ||
      Call->hasFnAttr(Attribute::StrictFP) ||
      Call->getCallingConv() != CallingConv::C)
    return false;
  if (const auto *II = dyn_cast<IntrinsicInst>(Call)) {
    switch (II->getIntrinsicID()) {
    case Intrinsic::exp:
    case Intrinsic::exp2:
    case Intrinsic::log:
      return !isa<Constant>(II->getArgOperand(0));
    case Intrinsic::pow:
      return !isa<Constant>(II->getArgOperand(0)) &&
             !isa<Constant>(II->getArgOperand(1));
    default:
      return false;
    }
  }
  const auto *Callee = Call->getCalledFunction();
  if (!Callee || !Callee->isDeclaration() || Callee->isVarArg())
    return false;
  const StringRef Name = Callee->getName();
  return ((Name == "pow" && Call->arg_size() == 2) ||
          ((Name == "exp" || Name == "exp2" || Name == "expm1" || Name == "log") &&
           Call->arg_size() == 1)) &&
         all_of(Call->args(), [](const Use &U) { return U->getType()->isDoubleTy(); });
}

namespace {

constexpr char PassName[] = "rssn-libm-pressure";
constexpr unsigned Window = 28;
constexpr unsigned MaxCone = 10;
constexpr unsigned MaxRuns = 16;
constexpr unsigned MaxTrials = 192;

bool replayable(const Instruction &I) {
  if (!I.getType()->isDoubleTy())
    return false;
  switch (I.getOpcode()) {
  case Instruction::FAdd:
  case Instruction::FSub:
  case Instruction::FMul:
  case Instruction::FDiv:
  case Instruction::FNeg:
    return isSafeToSpeculativelyExecute(&I);
  default:
    return false;
  }
}

using Nodes = SmallPtrSet<Instruction *, MaxCone>;

bool collect(Instruction &I, BasicBlock &B, const BlockPressure &P,
             unsigned Cut, Nodes &Cone) {
  if (I.getParent() != &B || P.Positions.lookup(&I) < Cut)
    return true; // Existing input that dominates the cut.
  if (Cone.contains(&I))
    return true;
  if (!replayable(I) || Cone.size() == MaxCone)
    return false; // Includes the result of the libcall itself.
  Cone.insert(&I);
  for (const Use &U : I.operands())
    if (auto *Def = dyn_cast<Instruction>(U.get()))
      if (!collect(*Def, B, P, Cut, Cone))
        return false;
  return true;
}

int arithmeticLatency(const Instruction &I, const TargetTransformInfo &TTI) {
  const auto Cost = TTI.getInstructionCost(&I, TargetTransformInfo::TCK_Latency);
  if (I.getOpcode() == Instruction::FDiv)
    return Cost.isValid() ? std::clamp<int64_t>(Cost.getValue(), 4, 24) : 12;
  return Cost.isValid() ? std::clamp<int64_t>(Cost.getValue(), 1, 12) : 4;
}

struct Choice {
  Instruction *Cut = nullptr;
  SmallVector<Instruction *, MaxCone> Order;
  int Gain = 0;
  int RemovedSIMD = 0;
};

void evaluate(BasicBlock &B, const Liveness &Live, const BlockPressure &Old,
              unsigned Cut, const Nodes &Cone, const TargetTransformInfo &TTI,
              Choice &Best) {
  SmallVector<Instruction *, MaxCone> Order(Cone.begin(), Cone.end());
  sort(Order, [&](const Instruction *A, const Instruction *B) {
    return Old.Positions.lookup(A) < Old.Positions.lookup(B);
  });
  SmallVector<Instruction *, 96> Code;
  Code.reserve(Old.Code.size());
  for (auto [Index, I] : enumerate(Old.Code)) {
    if (Index == Cut)
      append_range(Code, Order);
    if (!Cone.contains(I))
      Code.push_back(I);
  }
  BlockPressure Trial(B, Live, /*IncludeScalarLibm=*/true, Code);
  const Pressure Before = Old.Through[Cut];
  const Pressure After = Trial.Through[Cut + Order.size()];
  const int Saved = Before.SIMD - After.SIMD;
  if (Saved < 2 || After.GP > Before.GP)
    return;

  for (auto [Point, I] : enumerate(Trial.Code)) {
    const auto Existing = Old.Positions.lookup(I);
    if (Trial.Before[Point].worsens(Old.Before[Existing], Live.Limit))
      return;
  }

  // Existing values that have crossed an earlier call may already reside
  // in a backing spill slot. Reward the cut shrink, but price the critical
  // chain and every newly introduced live result separately.
  int PathLatency = 0;
  int Outputs = 0;
  DenseMap<const Value *, int> Ready;
  for (Instruction *I : Order) {
    int Latest = 0;
    for (const Use &U : I->operands())
      Latest = std::max(Latest, Ready.lookup(U.get()));
    Ready[I] = Latest + arithmeticLatency(*I, TTI);
    PathLatency = std::max(PathLatency, Ready[I]);
    Outputs += any_of(I->users(), [&](const User *U) {
      return !Cone.contains(dyn_cast<Instruction>(U));
    });
  }
  const int Gain = Before.callCost(Liveness::CallLimit) -
                   After.callCost(Liveness::CallLimit) -
                   Outputs - (PathLatency + 3) / 4;
  if (Gain > Best.Gain)
    Best = {Old.Code[Cut], std::move(Order), Gain, Saved};
}

Choice choose(BasicBlock &B, const Liveness &Live,
              const TargetTransformInfo &TTI, unsigned &Trials) {
  Choice Best;
  BlockPressure P(B, Live, /*IncludeScalarLibm=*/true);
  for (unsigned Cut : P.Calls.set_bits()) {
    if (!isScalarLibmCut(*P.Code[Cut]) || P.Through[Cut].SIMD < 3)
      continue;
    const unsigned End = std::min<unsigned>(Cut + 1 + Window, P.Code.size());
    for (unsigned J = Cut + 1; J < End; ++J) {
      Instruction *I = P.Code[J];
      if (isa<CallBase>(I))
        break;
      if (!replayable(*I) || I->use_empty())
        continue;
      Nodes Cone;
      if (!collect(*I, B, P, Cut, Cone) || Cone.size() < 2)
        continue;
      if (++Trials > MaxTrials)
        return Best;
      evaluate(B, Live, P, Cut, Cone, TTI, Best);
    }
  }
  return Best;
}

struct Guard {
  IntrinsicInst *Call;
  SelectInst *Select;
  Value *Condition;
  bool OnTrue;
  unsigned Lanes;
  MDNode *Weights;
};

std::optional<Guard> findGuard(SelectInst &Select, const TargetLibraryInfo &TLI) {
  auto *VT = dyn_cast<FixedVectorType>(Select.getType());
  if (!VT || !VT->getElementType()->isDoubleTy())
    return std::nullopt;
  const unsigned Lanes = VT->getNumElements();
  if (Lanes != 2 && Lanes != 4 && Lanes != 8)
    return std::nullopt;
  // The original vector condition must select the same edge in all defined
  // lanes. Any poison lanes are refined by the frozen scalar branch below.
  Value *Cond = Select.getCondition();
  if (!Cond->getType()->isIntegerTy(1))
    Cond = getSplatValue(Cond);
  if (!Cond || !Cond->getType()->isIntegerTy(1) || isa<Constant>(Cond))
    return std::nullopt;

  for (unsigned Edge = 1; Edge <= 2; ++Edge) {
    auto *Call = dyn_cast<IntrinsicInst>(Select.getOperand(Edge));
    if (!Call || !Call->hasOneUse() || Call->getParent() != Select.getParent() ||
        !Call->comesBefore(&Select) || Call->isConvergent() ||
        Call->hasOperandBundles() || !Call->doesNotAccessMemory() ||
        !Call->hasFnAttr(Attribute::Speculatable))
      continue;
    StringRef Scalar;
    switch (Call->getIntrinsicID()) {
    case Intrinsic::exp:
      Scalar = "exp";
      break;
    case Intrinsic::log:
      Scalar = "log";
      break;
    default:
      continue;
    }
    if (isa<Constant>(Call->getArgOperand(0)) || Call->hasApproxFunc())
      continue;
    const std::string IntrinsicName = (Twine("llvm.") + Scalar + ".f64").str();
    bool HasVectorLibrary = false;
    for (unsigned Width = 2; Width <= Lanes; Width *= 2) {
      const auto VF = ElementCount::getFixed(Width);
      HasVectorLibrary |= TLI.isFunctionVectorizable(Scalar, VF) ||
                          TLI.isFunctionVectorizable(IntrinsicName, VF);
    }
    if (HasVectorLibrary)
      continue;
    uint64_t TrueWeight, FalseWeight;
    MDNode *Weights = nullptr;
    if (extractBranchWeights(Select, TrueWeight, FalseWeight)) {
      const uint64_t SkipWeight = Edge == 1 ? FalseWeight : TrueWeight;
      // On explicitly profiled edges, require at least one avoided scalar
      // libcall per two visits in expectation to pay for the new branch.
      if (2 * Lanes * SkipWeight < TrueWeight + FalseWeight || !SkipWeight)
        continue;
      Weights = Select.getMetadata(LLVMContext::MD_prof);
    }
    return Guard{Call, &Select, Cond, Edge == 1, Lanes, Weights};
  }
  return std::nullopt;
}

bool canShareGuard(ArrayRef<Guard> Group, const Guard &Next) {
  const Guard &Previous = Group.back();
  if (Previous.Condition != Next.Condition || Previous.OnTrue != Next.OnTrue ||
      Previous.Weights != Next.Weights ||
      Previous.Select->getParent() != Next.Select->getParent())
    return false;

  // Keep the packing local. Between the selected results only allow mask and
  // result assembly, along with the new guarded call.
  for (const Instruction *I = Previous.Select->getNextNode(); I != Next.Select;
       I = I->getNextNode())
    if (I != Next.Call && !isa<InsertElementInst, ShuffleVectorInst>(I))
      return false;

  // A result used before the new join cannot be moved there. PHI users in
  // other blocks are valid because the split updates the predecessor edges.
  const BasicBlock *Block = Next.Select->getParent();
  for (const Guard &G : Group)
    for (const User *U : G.Select->users())
      if (const auto *I = dyn_cast<Instruction>(U);
          I && I->getParent() == Block && !Next.Select->comesBefore(I))
        return false;
  return true;
}

bool guardVectorLibm(Function &F, FunctionAnalysisManager &AM, bool Transform) {
  auto &TLI = AM.getResult<TargetLibraryAnalysis>(F);
  SmallVector<SelectInst *, 8> Candidates;
  for (BasicBlock &B : F)
    if (!B.isEHPad())
      for (Instruction &I : B)
        if (auto *S = dyn_cast<SelectInst>(&I))
          if (S->getType()->isVectorTy())
            Candidates.push_back(S);

  SmallVector<Guard, 8> Plans;
  auto &ORE = AM.getResult<OptimizationRemarkEmitterAnalysis>(F);
  for (SelectInst *S : Candidates) {
    if (Plans.size() == 12)
      break;
    auto G = findGuard(*S, TLI);
    if (!G)
      continue;
    ORE.emit([&] {
      using NV = ore::NV;
      return OptimizationRemarkAnalysis(PassName, "UniformVectorGuard", S)
             << NV("Function", F.getName()) << " "
             << NV("Action", Transform ? "guarded" : "candidate") << " "
             << NV("Libm", G->Call->getCalledFunction()->getName())
             << " with " << NV("ScalarLanes", G->Lanes) << " scalarized lanes";
    });
    Plans.push_back(*G);
  }
  if (!Transform)
    return false;

  ArrayRef<Guard> Remaining = Plans;
  while (!Remaining.empty()) {
    unsigned Count = 1;
    while (Count < Remaining.size() &&
           canShareGuard(Remaining.take_front(Count), Remaining[Count]))
      ++Count;
    ArrayRef<Guard> Group = Remaining.take_front(Count);
    Remaining = Remaining.drop_front(Count);

    const Guard &Last = Group.back();
    SelectInst *S = Last.Select;
    BasicBlock *Head = S->getParent();
    IRBuilder<> Builder(S);
    Value *Cond = Builder.CreateFreeze(Last.Condition, "libm.condition");
    if (!Last.OnTrue)
      Cond = Builder.CreateNot(Cond, "libm.active");
    Instruction *Jump = SplitBlockAndInsertIfThen(Cond, S, false, Last.Weights);
    if (Last.Weights && !Last.OnTrue)
      Head->getTerminator()->swapProfMetadata();
    BasicBlock *Evaluate = Jump->getParent();
    Evaluate->setName("libm.evaluate");
    BasicBlock *Join = S->getParent();
    for (const Guard &G : Group) {
      Value *Other = G.Select->getOperand(G.OnTrue ? 2 : 1);
      salvageDebugInfo(*G.Call);
      G.Call->moveBefore(Jump->getIterator());
      auto *Result = PHINode::Create(G.Select->getType(), 2, "", Join->begin());
      Result->takeName(G.Select);
      Result->setDebugLoc(G.Select->getDebugLoc());
      Result->addIncoming(Other, Head);
      Result->addIncoming(G.Call, Evaluate);
      G.Select->replaceAllUsesWith(Result);
      G.Select->eraseFromParent();
    }
  }
  return !Plans.empty();
}

} // namespace

PreservedAnalyses LibmPressurePass::run(Function &F, FunctionAnalysisManager &AM) {
  const Triple &T = F.getParent()->getTargetTriple();
  const bool SysV = T.isOSLinux() && T.getEnvironment() != Triple::GNUX32;
  const bool Win64 = T.isWindowsMSVCEnvironment();
  if (!F.getName().contains(Filter) || (!SysV && !Win64) ||
      T.getArch() != Triple::x86_64 ||
      F.isDeclaration() || F.hasOptNone() || F.hasOptSize() || F.hasMinSize() ||
      F.hasFnAttribute(Attribute::StrictFP) || F.hasGC())
    return PreservedAnalyses::all();

  if (!TrialCones && F.size() <= 1024 && F.getInstructionCount() <= 24000 &&
      guardVectorLibm(F, AM, Transform))
    return PreservedAnalyses::none();

  if (!SysV || (Transform && !TrialCones) || !Liveness::supports(F))
    return PreservedAnalyses::all();
  auto &LI = AM.getResult<LoopAnalysis>(F);
  if (LI.empty() || !any_of(F, [&](const BasicBlock &B) {
        return LI.getLoopFor(&B) && any_of(B, isScalarLibmCut);
      }))
    return PreservedAnalyses::all();
  auto &TTI = AM.getResult<TargetIRAnalysis>(F);
  Liveness Live(TTI, F);
  if (!Live.compute(F))
    return PreservedAnalyses::all();

  auto &ORE = AM.getResult<OptimizationRemarkEmitterAnalysis>(F);
  unsigned Changes = 0, Trials = 0;
  for (BasicBlock &B : F) {
    if (!LI.getLoopFor(&B) || B.isEHPad() || !any_of(B, isScalarLibmCut))
      continue;
    for (unsigned J = 0; J < 4 && Changes < MaxRuns && Trials < MaxTrials; ++J) {
      Choice Plan = choose(B, Live, TTI, Trials);
      if (!Plan.Cut)
        break;
      ORE.emit([&] {
        using NV = ore::NV;
        const auto *Callee = cast<CallInst>(Plan.Cut)->getCalledFunction();
        return OptimizationRemarkAnalysis(PassName, "ConsumerCone", Plan.Cut)
               << NV("Function", F.getName()) << " "
               << NV("Action", Transform ? "scheduled" : "candidate") << " "
               << NV("Operations", Plan.Order.size()) << " operations before "
               << NV("Libm", Callee->getName()) << ", SIMD cut reduced by "
               << NV("RemovedSIMD", Plan.RemovedSIMD) << ", score "
               << NV("Score", Plan.Gain);
      });
      if (!Transform)
        break;
      for (Instruction *I : Plan.Order) {
        salvageDebugInfo(*I);
        I->moveBefore(Plan.Cut->getIterator());
      }
      ++Changes;
    }
  }
  if (!Changes)
    return PreservedAnalyses::all();
  PreservedAnalyses PA;
  PA.preserveSet<CFGAnalyses>();
  return PA;
}

} // namespace rssn
