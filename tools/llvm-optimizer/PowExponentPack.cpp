#include "PowExponentPack.h"
#include "LibmPressure.h"
#include "Pressure.h"
#include "llvm/ADT/MapVector.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/ADT/bit.h"
#include "llvm/Analysis/BlockFrequencyInfo.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/OptimizationRemarkEmitter.h"
#include "llvm/Analysis/TargetTransformInfo.h"
#include "llvm/IR/CFG.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/InlineAsm.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/Module.h"
#include <optional>

using namespace llvm;

namespace rssn {
namespace {

constexpr char PassName[] = "rssn-pow-exponent-pack";
constexpr unsigned MaxValues = 12;
constexpr unsigned MaxBetween = 12;

struct Window {
  IntrinsicInst *First;
  IntrinsicInst *Second;
  BasicBlock *Preheader;
  unsigned SharedCuts = 0;
  double Frequency = 0;
};

using Exponents = std::pair<Value *, Value *>;
using Windows = SmallVector<Window, 4>;

struct Plan {
  Exponents Values;
  Instruction *Insert;
  MutableArrayRef<Window> Uses;
  SmallVector<const Instruction *, 16> Relief;
  double SetupFrequency;
  double BypassedLowUses = 0;
  double Score = 0;
  unsigned Mask = 0;
};

using CallBudget = DenseMap<const Instruction *, int>;

struct PackingCosts {
  double Setup;
  double HighLane;
  double Reload;
};

std::optional<PackingCosts> packingCosts(const TargetTransformInfo &TTI, LLVMContext &C) {
  auto *Double = Type::getDoubleTy(C);
  auto *VT = FixedVectorType::get(Double, 2);
  const auto Kind = TargetTransformInfo::TCK_RecipThroughput;
  const auto Insert = TTI.getVectorInstrCost(Instruction::InsertElement, VT, Kind, 0) +
                      TTI.getVectorInstrCost(Instruction::InsertElement, VT, Kind, 1);
  const auto High = TTI.getShuffleCost(TargetTransformInfo::SK_PermuteSingleSrc,
                                      VT, VT, {1, 1}, Kind);
  const auto Read = TTI.getMemoryOpCost(Instruction::Load, Double, Align(8), 0, Kind);
  if (!Insert.isValid() || !High.isValid() || !Read.isValid())
    return std::nullopt;
  // An estimated dependent scalar reload has both issue and latency costs.
  return PackingCosts{static_cast<double>(std::max<int64_t>(2, Insert.getValue())),
                      static_cast<double>(std::max<int64_t>(1, High.getValue())),
                      2.0 * std::max<int64_t>(1, Read.getValue())};
}

bool scalarPow(const Instruction &I) {
  const auto *II = dyn_cast<IntrinsicInst>(&I);
  return II && II->getIntrinsicID() == Intrinsic::pow &&
         isScalarLibmCut(I) && !II->hasApproxFunc();
}

bool endsWindow(const Instruction &I) {
  return isa<CallBase>(I) || I.mayReadOrWriteMemory() || I.mayHaveSideEffects();
}

IntrinsicInst *firstPow(BasicBlock &B, const Loop &L, unsigned Distance) {
  for (Instruction &I : B) {
    if (isa<DbgInfoIntrinsic>(I))
      continue;
    if (scalarPow(I) && L.isLoopInvariant(cast<CallInst>(I).getArgOperand(1)))
      return cast<IntrinsicInst>(&I);
    if (I.isTerminator() || endsWindow(I) || ++Distance > MaxBetween)
      break;
  }
  return nullptr;
}

MapVector<Exponents, Windows> collectWindows(Function &F, const LoopInfo &LI) {
  MapVector<Exponents, Windows> Groups;
  for (BasicBlock &B : F) {
    const Loop *L = LI.getLoopFor(&B);
    BasicBlock *Preheader = L ? L->getLoopPreheader() : nullptr;
    if (!Preheader || B.isEHPad())
      continue;
    const auto Add = [&](IntrinsicInst *First, IntrinsicInst *Second) {
      Value *A = First->getArgOperand(1);
      Value *B = Second->getArgOperand(1);
      if (A != B)
        Groups[{A, B}].push_back({First, Second, Preheader});
    };
    IntrinsicInst *Previous = nullptr;
    unsigned Distance = 0;
    for (Instruction &I : B) {
      if (isa<DbgInfoIntrinsic>(I))
        continue;
      if (I.isTerminator())
        break;
      if (scalarPow(I) && L->isLoopInvariant(cast<CallInst>(I).getArgOperand(1))) {
        auto *Next = cast<IntrinsicInst>(&I);
        if (Previous)
          Add(Previous, Next);
        Previous = Next;
        Distance = 0;
      } else if (endsWindow(I) || ++Distance > MaxBetween) {
        Previous = nullptr;
      }
    }
    if (!Previous || !isa<CondBrInst, UncondBrInst>(B.getTerminator()))
      continue;
    for (BasicBlock *Next : successors(&B)) {
      if (Next == &B || Next == L->getHeader() || Next->isEHPad() ||
          Next->getSinglePredecessor() != &B || LI.getLoopFor(Next) != L)
        continue;
      if (auto *First = firstPow(*Next, *L, Distance))
        Add(Previous, First);
    }
  }
  return Groups;
}

std::optional<Plan> planPair(Exponents Values, MutableArrayRef<Window> Sites,
                             const DominatorTree &DT, const BlockFrequencyInfo &BFI) {
  if (Sites.size() < 2 || Sites.size() > 32)
    return std::nullopt;
  BasicBlock *Materialize = Sites.front().Preheader;
  SmallPtrSet<const BasicBlock *, 8> Blocks;
  SmallPtrSet<const Use *, 32> Replaced;
  for (const Window &W : Sites) {
    Materialize = DT.findNearestCommonDominator(Materialize, W.Preheader);
    if (!Materialize)
      return std::nullopt;
    Blocks.insert(W.First->getParent());
    Replaced.insert(&W.First->getArgOperandUse(1));
    Replaced.insert(&W.Second->getArgOperandUse(1));
  }
  if (Materialize->isEHPad() || Blocks.size() < 2)
    return std::nullopt;
  Instruction *Insert = Materialize->getTerminator();
  for (Value *V : {Values.first, Values.second}) {
    if (!DT.dominates(V, Insert))
      return std::nullopt;
    // An exponent with a use after materialization would keep its original
    // scalar interval alongside the packed interval. Count every SSA use.
    for (const Use &U : V->uses()) {
      if (Replaced.contains(&U))
        continue;
      const auto *I = dyn_cast<Instruction>(U.getUser());
      if (!I || I == Insert || !DT.dominates(I, Insert))
        return std::nullopt;
    }
  }
  return Plan{Values, Insert, Sites, {},
              static_cast<double>(BFI.getBlockFreq(Materialize).getFrequency())};
}

bool liveThrough(Exponents E, const BlockPressure &P, const Liveness &Live,
                 const BasicBlock &B, unsigned Point) {
  const BitVector &Out = Live.liveOut(B);
  const auto Crosses = [&](const Value *V) {
    const unsigned ID = Live.id(V);
    // Exponents dominate this call. One can be live through the entire
    // block without occurring in the block's local use table.
    return P.isLiveBefore(V, Point + 1) || (ID < Out.size() && Out.test(ID));
  };
  return Crosses(E.first) && Crosses(E.second);
}

CallBudget measureCuts(MutableArrayRef<Plan> Plans, const Liveness &Live,
                       const BlockFrequencyInfo &BFI, const DominatorTree &DT) {
  struct WindowCut {
    Window *Owner;
    IntrinsicInst *Call;
    unsigned Bit;
  };
  DenseMap<BasicBlock *, SmallVector<WindowCut, 8>> Blocks;
  for (Plan &Candidate : Plans) {
    DenseMap<IntrinsicInst *, double> Bypasses;
    for (Window &W : Candidate.Uses) {
      auto *First = W.First->getParent();
      auto *Second = W.Second->getParent();
      const double FirstFrequency = BFI.getBlockFreq(First).getFrequency();
      const double SecondFrequency = BFI.getBlockFreq(Second).getFrequency();
      W.Frequency = std::min(FirstFrequency, SecondFrequency);
      if (First != Second) {
        auto It = Bypasses.try_emplace(W.First, FirstFrequency).first;
        It->second -= W.Frequency;
      }
      Blocks[First].push_back({&W, W.First, 1U});
      Blocks[Second].push_back({&W, W.Second, 2U});
    }
    // Charge a potential copy into the first call's argument register when
    // it executes without a high-lane consumer on the selected edges.
    for (const auto &[Call, Frequency] : Bypasses)
      Candidate.BypassedLowUses += std::max(0.0, Frequency);
  }

  CallBudget Budget;
  for (const auto &[B, Sites] : Blocks) {
    BlockPressure P(*B, Live, /*IncludeScalarLibm=*/true, {}, CallPressureABI::Win64);
    for (const auto &[W, Call, Bit] : Sites) {
      const Exponents E{W->First->getArgOperand(1), W->Second->getArgOperand(1)};
      const unsigned Point = P.Positions.lookup(Call);
      if (!P.Calls.test(Point))
        continue;
      if (liveThrough(E, P, Live, *B, Point))
        W->SharedCuts |= Bit;
      const int Excess = P.Through[Point].PreservableSIMD - Pressure::Win64PreservedXMM;
      if (Excess <= 0 || !Budget.try_emplace(Call, Excess).second)
        continue;
      for (Plan &Candidate : Plans)
        if (DT.dominates(Candidate.Insert, Call) &&
            liveThrough(Candidate.Values, P, Live, *B, Point))
          Candidate.Relief.push_back(Call);
    }
  }
  return Budget;
}

double score(const Plan &P, const CallBudget &Budget, PackingCosts Costs) {
  double Gain = -Costs.Setup * P.SetupFrequency - P.BypassedLowUses;
  for (const Window &W : P.Uses) {
    const bool SavesReload = ((W.SharedCuts & 1U) && Budget.lookup(W.First) > 0) ||
                             ((W.SharedCuts & 2U) && Budget.lookup(W.Second) > 0);
    // Both calls can share the same carrier reload. Credit at most one scalar
    // reload per window, and charge the high-lane shuffle at every use.
    Gain += W.Frequency * ((SavesReload ? Costs.Reload : 0) - Costs.HighLane);
  }
  return Gain;
}

unsigned assignMasks(MutableArrayRef<Plan> Plans) {
  DenseMap<Value *, unsigned> Values;
  for (Plan &P : Plans) {
    const unsigned A = Values.try_emplace(P.Values.first, Values.size()).first->second;
    const unsigned B = Values.try_emplace(P.Values.second, Values.size()).first->second;
    if (Values.size() > MaxValues)
      return 0;
    P.Mask = (1U << A) | (1U << B);
  }
  return 1U << Values.size();
}

SmallVector<unsigned, 6> bestMatching(ArrayRef<Plan> Plans, unsigned States) {
  SmallVector<double, 64> Scores(States);
  SmallVector<int, 64> Choices(States, -1);
  for (unsigned Mask = 1; Mask < States; ++Mask) {
    const unsigned Bit = Mask & -Mask;
    Scores[Mask] = Scores[Mask ^ Bit];
    for (auto [I, P] : enumerate(Plans)) {
      if ((P.Mask & Mask) != P.Mask || !(P.Mask & Bit))
        continue;
      const double Total = P.Score + Scores[Mask ^ P.Mask];
      if (Total > Scores[Mask]) {
        Scores[Mask] = Total;
        Choices[Mask] = I;
      }
    }
  }
  SmallVector<unsigned, 6> Selected;
  for (unsigned Mask = States - 1; Mask;) {
    if (Choices[Mask] < 0) {
      Mask ^= Mask & -Mask;
      continue;
    }
    const unsigned Index = Choices[Mask];
    Selected.push_back(Index);
    Mask ^= Plans[Index].Mask;
  }
  return Selected;
}

SmallVector<unsigned, 6> choosePairs(MutableArrayRef<Plan> Plans, CallBudget Budget,
                                     PackingCosts Costs) {
  SmallVector<unsigned, 6> Selected;
  const unsigned States = assignMasks(Plans);
  if (States <= 1)
    return Selected;

  unsigned Used = 0;
  for (;;) {
    for (Plan &P : Plans)
      P.Score = (P.Mask & Used) ? 0 : score(P, Budget, Costs);
    const auto Matching = bestMatching(Plans, States);
    if (Matching.empty())
      break;
    const unsigned Index = *std::max_element(Matching.begin(), Matching.end(),
        [&](unsigned A, unsigned B) { return Plans[A].Score < Plans[B].Score; });
    Selected.push_back(Index);
    const Plan &P = Plans[Index];
    Used |= P.Mask;
    // A second disjoint pair may relieve the very same oversubscribed cut.
    // Rescore before committing it against the residual preserved XMM budget.
    for (const Instruction *Call : P.Relief) {
      int &Remaining = Budget.find(Call)->second;
      if (Remaining > 0)
        --Remaining;
    }
  }
  return Selected;
}

void pack(const Plan &P) {
  IRBuilder<> Entry(P.Insert);
  auto *VT = FixedVectorType::get(Entry.getDoubleTy(), 2);
  Value *Carrier = Entry.CreateInsertElement(PoisonValue::get(VT), P.Values.first,
                                             uint64_t(0));
  Carrier = Entry.CreateInsertElement(Carrier, P.Values.second, uint64_t(1));
  Carrier = Entry.CreateFreeze(Carrier, "pow.exponents");
  auto *Signature = FunctionType::get(VT, {VT}, false);
  // Keep the register-class constraint at each consuming call. SelectionDAG
  // can otherwise split an insert/extract-only value back into scalar ranges.
  InlineAsm *Low = InlineAsm::get(Signature, "", "=x,0", false);
  InlineAsm *High = InlineAsm::get(Signature, "vunpckhpd $1, $1, $0", "=x,x", false);
  SmallPtrSet<IntrinsicInst *, 16> Rewritten;
  for (const Window &W : P.Uses) {
    for (auto [Call, Assembly] :
         {std::pair(W.First, Low), std::pair(W.Second, High)}) {
      if (!Rewritten.insert(Call).second)
        continue;
      IRBuilder<> AtCall(Call);
      auto *Register = AtCall.CreateCall(Assembly, {Carrier}, "pow.exponent.reg");
      Register->setDoesNotAccessMemory();
      Register->setDoesNotThrow();
      Register->addFnAttr(Attribute::WillReturn);
      Call->setArgOperand(1, AtCall.CreateExtractElement(Register, uint64_t(0)));
    }
  }
}

} // namespace

PreservedAnalyses PowExponentPackPass::run(Function &F, FunctionAnalysisManager &AM) {
  const Triple &T = F.getParent()->getTargetTriple();
  if (T.getArch() != Triple::x86_64 || !T.isWindowsMSVCEnvironment() ||
      !F.getName().contains(Filter) || F.isDeclaration() || F.hasOptNone() ||
      F.hasOptSize() || F.hasMinSize() || F.hasFnAttribute(Attribute::StrictFP) ||
      F.hasFnAttribute(Attribute::NoImplicitFloat) || !Liveness::supports(F))
    return PreservedAnalyses::all();
  const auto &TTI = AM.getResult<TargetIRAnalysis>(F);
  if (TTI.getRegisterBitWidth(TargetTransformInfo::RGK_FixedWidthVector)
          .getKnownMinValue() < 256)
    return PreservedAnalyses::all();

  auto Groups = collectWindows(F, AM.getResult<LoopAnalysis>(F));
  if (Groups.empty() || Groups.size() > 32)
    return PreservedAnalyses::all();
  const auto Costs = packingCosts(TTI, F.getContext());
  if (!Costs)
    return PreservedAnalyses::all();
  Liveness Live(TTI, F);
  if (!Live.compute(F))
    return PreservedAnalyses::all();

  const auto &DT = AM.getResult<DominatorTreeAnalysis>(F);
  const auto &BFI = AM.getResult<BlockFrequencyAnalysis>(F);
  SmallVector<Plan, 8> Plans;
  for (auto &[Key, Sites] : Groups)
    if (auto P = planPair(Key, Sites, DT, BFI))
      Plans.push_back(*P);
  auto Budget = measureCuts(Plans, Live, BFI, DT);
  const auto Selected = choosePairs(Plans, std::move(Budget), *Costs);
  auto &ORE = AM.getResult<OptimizationRemarkEmitterAnalysis>(F);
  for (unsigned I : Selected) {
    const Plan &P = Plans[I];
    ORE.emit([&] {
      using NV = ore::NV;
      unsigned Shared = 0;
      unsigned CrossEdge = 0;
      for (const Window &W : P.Uses) {
        Shared += llvm::popcount(W.SharedCuts);
        CrossEdge += W.First->getParent() != W.Second->getParent();
      }
      return OptimizationRemarkAnalysis(PassName, "PackedExponents", P.Insert)
             << NV("Function", F.getName()) << " "
             << NV("Action", Transform ? "packed" : "can pack")
             << " two pow exponents at " << NV("Windows", P.Uses.size()) << " windows"
             << ", cross-edge windows " << NV("CrossEdge", CrossEdge)
             << ", overlapping calls " << NV("SharedCuts", Shared);
    });
    if (Transform)
      pack(P);
  }
  if (!Transform || Selected.empty())
    return PreservedAnalyses::all();
  PreservedAnalyses PA;
  PA.preserveSet<CFGAnalyses>();
  return PA;
}

} // namespace rssn
