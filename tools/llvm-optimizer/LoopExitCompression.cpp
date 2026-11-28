#include "LoopExitCompression.h"
#include "Pressure.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/Analysis/BlockFrequencyInfo.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/OptimizationRemarkEmitter.h"
#include "llvm/Analysis/TargetTransformInfo.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/Module.h"
#include "llvm/Transforms/Utils/Local.h"
#include <optional>

using namespace llvm;

namespace rssn {
namespace {

constexpr char PassName[] = "rssn-loop-exit-compression";
constexpr unsigned MaxOps = 6;
constexpr unsigned MaxInputs = 7;
constexpr unsigned MaxTrials = 12;
constexpr unsigned MaxExitScan = 32;

struct Plan {
  Instruction *Insert = nullptr;
  SmallVector<Instruction *, MaxOps> Ops;
  unsigned ImprovedCalls = 0;
  int MaxReduction = 0;
  double Score = 0;
};

bool foldable(const Instruction &I) {
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

using Inputs = SmallPtrSet<Value *, MaxInputs>;

bool collect(Instruction &I, BasicBlock &Exit,
             SmallVectorImpl<Instruction *> &Ops, Inputs &Leaves) {
  if (is_contained(Ops, &I))
    return true;
  for (Value *V : I.operand_values()) {
    if (isa<Constant>(V))
      continue;
    auto *Inner = dyn_cast<Instruction>(V);
    if (Inner && Inner->getParent() == &Exit && foldable(*Inner)) {
      if (!collect(*Inner, Exit, Ops, Leaves))
        return false;
    } else {
      Leaves.insert(V);
      if (Leaves.size() > MaxInputs)
        return false;
    }
  }
  if (Ops.size() == MaxOps)
    return false;
  Ops.push_back(&I);
  return true;
}

std::optional<Plan> buildPlan(Instruction &Root, BasicBlock &Preheader,
                              const DominatorTree &DT) {
  if (!foldable(Root) || Root.use_empty() || Root.getParent() == &Preheader ||
      !DT.dominates(&Preheader, Root.getParent()))
    return std::nullopt;
  Plan P;
  Inputs Leaves;
  if (!collect(Root, *Root.getParent(), P.Ops, Leaves) || P.Ops.size() < 2 ||
      Leaves.size() < 2)
    return std::nullopt;

  for (Instruction *I : ArrayRef(P.Ops).drop_back())
    if (!all_of(I->users(), [&](const User *U) {
          return is_contained(P.Ops, dyn_cast<Instruction>(U));
        }))
      return std::nullopt;

  P.Insert = &*Preheader.getFirstInsertionPt();
  for (Value *V : Leaves) {
    if (!DT.dominates(V, Preheader.getTerminator()))
      return std::nullopt;
    auto *Def = dyn_cast<Instruction>(V);
    if (Def && Def->getParent() == &Preheader && !Def->comesBefore(P.Insert))
      P.Insert = Def->getNextNode();
  }
  return P;
}

int spillCost(const Instruction &I, Pressure P, bool Win64) {
  const auto *Call = cast<CallInst>(&I);
  const CallingConv::ID CC = Call->getCallingConv();
  if (CC == CallingConv::C || CC == CallingConv::Fast)
    return Win64 ? 5 * std::max(0, P.SIMD - 10) : P.SIMDSpill;
  if (CC == CallingConv::X86_64_SysV)
    return P.SIMDSpill;
  return 0;
}

struct Baseline {
  DenseMap<const Instruction *, Pressure> AtCall;
  DenseMap<const BasicBlock *, int> PeakSIMD;
};

Baseline baseline(Function &F, const Liveness &Live) {
  Baseline Result;
  for (BasicBlock &B : F) {
    BlockPressure P(B, Live, /*IncludeScalarLibm=*/true);
    for (unsigned Site : P.Calls.set_bits())
      Result.AtCall[P.Code[Site]] = P.Through[Site];
    int Peak = 0;
    for (Pressure Before : P.Before)
      Peak = std::max(Peak, Before.SIMD);
    Result.PeakSIMD[&B] = Peak;
  }
  return Result;
}

bool evaluate(Plan &P, Function &F, const Baseline &Before,
              const TargetTransformInfo &TTI, const BlockFrequencyInfo &BFI,
              bool Win64) {
  BasicBlock &Dest = *P.Insert->getParent();
  BasicBlock &Exit = *P.Ops.back()->getParent();
  SmallVector<Instruction *, 32> Early, Late;
  for (Instruction &I : Dest) {
    if (&I == P.Insert)
      append_range(Early, P.Ops);
    Early.push_back(&I);
  }
  for (Instruction &I : Exit)
    if (!is_contained(P.Ops, &I))
      Late.push_back(&I);
  const Liveness::Schedule Trial[] = {{&Dest, Early}, {&Exit, Late}};
  Liveness Live(TTI, F);
  if (!Live.compute(F, Trial))
    return false;

  for (BasicBlock &B : F) {
    ArrayRef<Instruction *> Code;
    if (&B == &Dest)
      Code = Early;
    else if (&B == &Exit)
      Code = Late;
    BlockPressure Pk(B, Live, /*IncludeScalarLibm=*/true, Code);
    if (&B == &Dest && any_of(Pk.Before, [&](Pressure Point) {
          return Point.SIMD >
                 std::max(Before.PeakSIMD.lookup(&B), Live.Limit.SIMD);
        }))
      return false;
    const double Frequency = BFI.getBlockFreq(&B).getFrequency();
    for (unsigned Site : Pk.Calls.set_bits()) {
      const Instruction &Call = *Pk.Code[Site];
      const Pressure Old = Before.AtCall.lookup(&Call);
      const Pressure New = Pk.Through[Site];
      const int Saved =
          spillCost(Call, Old, Win64) - spillCost(Call, New, Win64);
      if (Saved < 0)
        return false;
      if (Saved > 0) {
        ++P.ImprovedCalls;
        P.MaxReduction = std::max(P.MaxReduction, Old.SIMD - New.SIMD);
        P.Score += Saved * Frequency;
      }
    }
  }

  const double EntryFreq = BFI.getBlockFreq(&Dest).getFrequency();
  const double ExitFreq = BFI.getBlockFreq(&Exit).getFrequency();
  if (EntryFreq > ExitFreq) {
    InstructionCost Cost = 0;
    for (Instruction *I : P.Ops)
      Cost +=
          TTI.getInstructionCost(I, TargetTransformInfo::TCK_RecipThroughput);
    if (!Cost.isValid())
      return false;
    P.Score -= (EntryFreq - ExitFreq) * std::max<int64_t>(1, Cost.getValue());
  }
  return P.ImprovedCalls >= 2 && P.MaxReduction > 0 && P.Score > 0;
}
} // namespace

PreservedAnalyses LoopExitCompressionPass::run(Function &F,
                                               FunctionAnalysisManager &AM) {
  const Triple &T = F.getParent()->getTargetTriple();
  const bool Win64 = T.isWindowsMSVCEnvironment();
  if (T.getArch() != Triple::x86_64 ||
      (!Win64 && (!T.isOSLinux() || T.isX32())) ||
      !F.getName().contains(Filter) || F.isDeclaration() || F.hasOptNone() ||
      F.hasOptSize() || F.hasMinSize() ||
      F.hasFnAttribute(Attribute::StrictFP) ||
      F.hasFnAttribute(Attribute::NoImplicitFloat) || !Liveness::supports(F))
    return PreservedAnalyses::all();

  const auto &LI = AM.getResult<LoopAnalysis>(F);
  if (LI.empty())
    return PreservedAnalyses::all();
  const auto &DT = AM.getResult<DominatorTreeAnalysis>(F);
  SmallVector<Plan, MaxTrials> Plans;
  for (const Loop *L : LI) {
    BasicBlock *Preheader = L->getLoopPreheader();
    if (!Preheader || Preheader->isEHPad())
      continue;
    SmallVector<BasicBlock *, 4> Exits;
    L->getUniqueExitBlocks(Exits);
    for (BasicBlock *Exit : Exits) {
      if (Exit->isEHPad() || LI.getLoopFor(Exit))
        continue;
      unsigned Scanned = 0;
      for (Instruction &I : make_range(Exit->getFirstNonPHIIt(), Exit->end())) {
        if (auto P = buildPlan(I, *Preheader, DT))
          Plans.push_back(std::move(*P));
        if (Plans.size() == MaxTrials || ++Scanned == MaxExitScan)
          break;
      }
      if (Plans.size() == MaxTrials)
        break;
    }
    if (Plans.size() == MaxTrials)
      break;
  }
  if (Plans.empty())
    return PreservedAnalyses::all();

  const auto &TTI = AM.getResult<TargetIRAnalysis>(F);
  if (!TTI.isTypeLegal(Type::getDoubleTy(F.getContext())))
    return PreservedAnalyses::all();
  const auto &BFI = AM.getResult<BlockFrequencyAnalysis>(F);
  Liveness Live(TTI, F);
  if (!Live.compute(F))
    return PreservedAnalyses::all();
  const Baseline Before = baseline(F, Live);
  Plan *Best = nullptr;
  for (Plan &P : Plans)
    if (evaluate(P, F, Before, TTI, BFI, Win64) &&
        (!Best || P.Score > Best->Score))
      Best = &P;
  if (!Best)
    return PreservedAnalyses::all();
  AM.getResult<OptimizationRemarkEmitterAnalysis>(F).emit([&] {
    using NV = ore::NV;
    return OptimizationRemarkAnalysis(PassName, "PackedLiveness",
                                      Best->Ops.back())
           << NV("Function", F.getName()) << " "
           << NV("Action", Transform ? "moved" : "can move") << " "
           << NV("Ops", Best->Ops.size())
           << " exit operations before the loop, "
           << NV("Calls", Best->ImprovedCalls) << " improved call sites, up to "
           << NV("MaxReduction", Best->MaxReduction)
           << " fewer live SIMD values";
  });
  if (!Transform)
    return PreservedAnalyses::all();
  for (Instruction *I : Best->Ops) {
    salvageDebugInfo(*I);
    I->moveBefore(Best->Insert->getIterator());
  }
  PreservedAnalyses PA;
  PA.preserveSet<CFGAnalyses>();
  return PA;
}
} // namespace rssn
