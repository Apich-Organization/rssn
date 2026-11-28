#include "Planning.h"
#include "llvm/ADT/MapVector.h"
#include "llvm/ADT/STLExtras.h"
#include "llvm/Analysis/TargetTransformInfo.h"
#include "llvm/Transforms/Utils/Local.h"

using namespace llvm;

namespace rssn {
namespace {

constexpr unsigned MaxMembers = MaxPlanWidth;
constexpr unsigned MaxPool = 12;
constexpr unsigned BeamWidth = 8;
constexpr unsigned MaxRecipes = 96;
constexpr unsigned MaxWindow = 2 * MaxDistance + 16;

struct Motion {
  Instruction *Def;
  GetElementPtrInst *Address;
  unsigned OldPoint;
  unsigned NewPoint;
  unsigned Start;
  int Bucket;
  int LatencyCost;

  const Instruction &source() const { return Address ? *Address : *Def; }
};

using Plan = SmallVector<const Motion *, MaxMembers>;

struct InputTail {
  Pressure Weight;
  unsigned OldEnd;
  unsigned NewEnd;
  unsigned LastMember;
};

struct Evaluation {
  int Gain = 0;
  int Overflow = 0;

  bool accepted() const { return Overflow == 0 && Gain >= 8; }
  int rank() const { return Gain - 12 * Overflow; }
};

int excessIncrease(Pressure New, Pressure Old, Pressure Limit) {
  return std::max(0, New.GP - std::max(Old.GP, Limit.GP)) +
         std::max(0, New.SIMD - std::max(Old.SIMD, Limit.SIMD));
}

void addInterval(MutableArrayRef<Pressure> Events, unsigned Base,
                 unsigned Begin, unsigned End, Pressure Weight, bool Add) {
  Begin = std::max(Base, Begin);
  if (Begin > End)
    return;
  if (Add) {
    Events[Begin - Base] += Weight;
    Events[End - Base + 1] -= Weight;
  } else {
    Events[Begin - Base] -= Weight;
    Events[End - Base + 1] += Weight;
  }
}

bool independent(const Motion &A, const Motion &B, const BlockPressure &P) {
  if (A.Def == B.Def || A.Def == B.Address || B.Def == A.Address ||
      (A.Address && A.Address == B.Address))
    return false;
  auto touches = [&](const Motion &X, const Motion &Y) {
    auto fromOther = [&](const Value *V) {
      return V == Y.Def || (Y.Address && V == Y.Address);
    };
    return fromOther(P.Code[X.NewPoint]) ||
           any_of(X.source().operands(),
                  [&](const Use &U) { return fromOther(U.get()); });
  };
  return !touches(A, B) && !touches(B, A);
}

Evaluation evaluate(const Plan &Members, const BlockPressure &P,
                    const Liveness &Live) {
  const unsigned Begin =
      (*min_element(Members,
                    [](auto *A, auto *B) { return A->Start < B->Start; }))
          ->Start +
      1;
  const unsigned End = (*max_element(Members, [](auto *A, auto *B) {
                         return A->NewPoint < B->NewPoint;
                       }))->NewPoint;
  if (End - Begin >= MaxWindow)
    return {-1, 1};

  SmallVector<Pressure, MaxWindow + 2> Events(End - Begin + 2);
  SmallVector<Pressure, MaxWindow + 2> CallEvents(End - Begin + 2);
  SmallDenseMap<const Value *, InputTail, 8> Inputs;
  Evaluation E;
  for (auto [Index, M] : enumerate(Members)) {
    const Pressure Result = Live.weight(M->Def);
    addInterval(Events, Begin, M->OldPoint + 1, M->NewPoint, Result, false);
    addInterval(CallEvents, Begin, M->OldPoint + 1, M->NewPoint - 1, Result,
                false);
    if (M->Address) {
      const Pressure Addr = Live.weight(M->Address);
      addInterval(Events, Begin, M->Start + 1, M->OldPoint, Addr, false);
      addInterval(CallEvents, Begin, M->Start + 1, M->OldPoint, Addr, false);
    }
    E.Gain -= M->LatencyCost;
    for (const Use &U : M->source().operands()) {
      const Value *V = U.get();
      auto Old = P.LastUse.find(V);
      if (Old == P.LastUse.end() || Old->second >= M->NewPoint)
        continue;
      auto [It, New] = Inputs.try_emplace(
          V, InputTail{Live.weight(V), Old->second, M->NewPoint,
                       static_cast<unsigned>(Index)});
      if (!New && M->NewPoint >= It->second.NewEnd) {
        It->second.NewEnd = M->NewPoint;
        It->second.LastMember = Index;
      }
    }
  }
  for (const auto &Entry : Inputs) {
    const InputTail &T = Entry.second;
    addInterval(Events, Begin, T.OldEnd + 1, T.NewEnd, T.Weight, true);
    addInterval(CallEvents, Begin, T.OldEnd, T.NewEnd - 1, T.Weight, true);
  }

  Pressure Delta, CallDelta;
  for (unsigned K = Begin; K <= End; ++K) {
    Delta += Events[K - Begin];
    CallDelta += CallEvents[K - Begin];
    const Pressure Old = P.Before[K];
    Pressure New = Old;
    New += Delta;
    E.Overflow += excessIncrease(New, Old, Live.Limit);
    E.Gain += Old.congestion(Live.Limit) - New.congestion(Live.Limit);

    if (P.Calls.test(K)) {
      const Pressure CallOld = P.Through[K];
      Pressure CallNew = CallOld;
      CallNew += CallDelta;
      E.Overflow += excessIncrease(CallNew, CallOld, Liveness::CallLimit);
      E.Gain += CallWeight * (CallOld.callCost(Liveness::CallLimit) -
                              CallNew.callCost(Liveness::CallLimit));
    }

    // Replay the actual insertion order at a shared gap. Multiple definitions
    // can overlap the same extended operand or a transient GEP result.
    for (auto [Index, M] : enumerate(Members)) {
      if (M->NewPoint != K)
        continue;
      for (const auto &Entry : Inputs) {
        const InputTail &Tail = Entry.second;
        if (Tail.NewEnd == K && Tail.LastMember == Index)
          New -= Tail.Weight;
      }
      if (M->Address) {
        Pressure BeforeLoad = New;
        BeforeLoad += Live.weight(M->Address);
        E.Overflow += excessIncrease(BeforeLoad, Old, Live.Limit);
      }
      New += Live.weight(M->Def);
      E.Overflow += excessIncrease(New, Old, Live.Limit);
    }
  }
  return E;
}

// TTI costs are target units rather than calibrated processor cycles. Count
// independent IR operations in a bounded gap as scheduling opportunities.
int exposedLatency(const BlockPressure &P, unsigned At, unsigned FirstUse,
                   const Instruction &I, int Latency) {
  unsigned Slots = 0;
  for (unsigned K = std::max(At, FirstUse > 24 ? FirstUse - 24 : 0u);
       K < FirstUse; ++K) {
    if (P.Calls.test(K)) {
      Slots = 0;
      continue;
    }
    const Instruction *Other = P.Code[K];
    if (!isa<GetElementPtrInst, CastInst>(Other) &&
        !is_contained(Other->operand_values(), &I))
      ++Slots;
  }
  return std::max(0, Latency - static_cast<int>(Slots / 2));
}

int targetLatency(const Instruction &I, const TargetTransformInfo &TTI) {
  const auto Latency =
      TTI.getInstructionCost(&I, TargetTransformInfo::TCK_Latency);
  const auto Rate =
      TTI.getInstructionCost(&I, TargetTransformInfo::TCK_RecipThroughput);
  const int Default = isa<LoadInst>(I) ? 4 : 1;
  if (!Latency.isValid() || !Rate.isValid())
    return Default;
  return std::clamp<int64_t>(std::max(Latency.getValue(), Rate.getValue()),
                             Default, 14);
}

void addRecipes(BlockPressure &P, const Liveness &Live,
                const TargetTransformInfo &TTI, MemoryModel &Memory,
                const DenseSet<const Instruction *> &Moved,
                SmallVectorImpl<Motion> &Recipes) {
  for (auto [Index, I] : enumerate(P.Code)) {
    const unsigned Def = Index;
    auto First = P.FirstUse.find(I);
    if (Moved.contains(I) || First == P.FirstUse.end() ||
        First->second <= Def + 1)
      continue;
    const unsigned Use = First->second;
    const unsigned Bound = std::min(Use, Def + MaxDistance);
    const int LastCall = P.Calls.find_prev(Bound);
    const bool CrossesCall = LastCall > static_cast<int>(Def);
    if (Bound - Def < 8 && !CrossesCall)
      continue;
    auto *Load = dyn_cast<LoadInst>(I);
    if (Load ? !simpleLoad(*Load, Live) : !cheapPure(*I, Live, TTI))
      continue;

    auto *Addr =
        Load ? dyn_cast<GetElementPtrInst>(Load->getPointerOperand()) : nullptr;
    if (Addr &&
        (!Addr->hasOneUse() || Addr->getParent() != I->getParent() ||
         Def - P.Positions.lookup(Addr) > 8 || !cheapPure(*Addr, Live, TTI)))
      Addr = nullptr;

    const int Latency = targetLatency(*I, TTI);
    const unsigned Earliest = CrossesCall ? LastCall + 1 : Def + 2;
    const unsigned Padding = std::min<unsigned>(2 * Latency, Bound - Earliest);
    const unsigned Preferred = Bound - Padding;
    // Keep a late alternative in regions where scheduling slack and register
    // congestion favor different insertion positions.
    SmallVector<unsigned, 2> Points{Preferred};
    if (!Load && Preferred != Bound)
      Points.push_back(Bound);
    for (unsigned Insert : Points) {
      if (Load && !Memory.unchanged(*Load, P.Code, Def + 1, Insert))
        continue;
      const int Extra = exposedLatency(P, Insert, Use, *I, Latency) -
                        exposedLatency(P, Def + 1, Use, *I, Latency);
      const int Bucket =
          CrossesCall ? LastCall : -2 - static_cast<int>(Def / 64);
      if (Addr) {
        Recipes.push_back({I, Addr, Def, Insert, P.Positions.lookup(Addr),
                           Bucket, 4 * std::max(0, Extra)});
      }
      Recipes.push_back(
          {I, nullptr, Def, Insert, Def, Bucket, 4 * std::max(0, Extra)});
      if (Recipes.size() >= MaxRecipes)
        return;
    }
  }
}

struct PoolEntry {
  const Motion *Recipe;
  Evaluation Score;
};

struct BeamEntry {
  Plan Members;
  Evaluation Score;
  unsigned Next;
};

void consider(Plan &Best, Evaluation &BestScore, const Plan &Candidate,
              Evaluation Score) {
  if (Score.accepted() && Score.Gain > BestScore.Gain) {
    Best = Candidate;
    BestScore = Score;
  }
}

void extendPool(ArrayRef<PoolEntry> Pool, unsigned Limit,
                const BlockPressure &P, const Liveness &Live, Plan &Best,
                Evaluation &BestScore) {
  SmallVector<BeamEntry, MaxPool> Beam;
  for (auto [I, Entry] : enumerate(Pool))
    Beam.push_back(
        {Plan{Entry.Recipe}, Entry.Score, static_cast<unsigned>(I + 1)});
  for (unsigned Depth = 2; Depth <= Limit && !Beam.empty(); ++Depth) {
    SmallVector<BeamEntry, 32> Next;
    for (const auto &Partial : Beam) {
      for (unsigned I = Partial.Next; I < Pool.size(); ++I) {
        const Motion *M = Pool[I].Recipe;
        if (any_of(Partial.Members, [&](const Motion *Old) {
              return !independent(*M, *Old, P);
            }))
          continue;
        Plan Candidate = Partial.Members;
        Candidate.push_back(M);
        const Evaluation Score = evaluate(Candidate, P, Live);
        consider(Best, BestScore, Candidate, Score);
        if (Depth < Limit)
          Next.push_back({std::move(Candidate), Score, I + 1});
      }
    }
    stable_sort(Next, [](const BeamEntry &A, const BeamEntry &B) {
      return A.Score.rank() > B.Score.rank();
    });
    if (Next.size() > BeamWidth)
      Next.resize(BeamWidth);
    Beam.assign(Next.begin(), Next.end());
  }
}

} // namespace

bool trySink(BlockPressure &P, const Liveness &Live,
             const TargetTransformInfo &TTI, MemoryModel &Memory,
             DenseSet<const Instruction *> &Moved, Changes &Result,
             unsigned RewriteBudget, unsigned PlanWidth) {
  if (!P.Pressured)
    return false;
  SmallVector<Motion, 24> Recipes;
  addRecipes(P, Live, TTI, Memory, Moved, Recipes);
  Plan Best;
  Evaluation BestScore;
  // Keep pool traversal deterministic independently of pointer values.
  MapVector<int, SmallVector<PoolEntry, MaxPool>> Pools;
  for (const Motion &M : Recipes) {
    const Plan Single{&M};
    const Evaluation Score = evaluate(Single, P, Live);
    consider(Best, BestScore, Single, Score);
    auto &Pool = Pools[M.Bucket];
    auto SameDef = find_if(Pool, [&](const PoolEntry &E) {
      return E.Recipe->Def == M.Def && E.Recipe->Address == M.Address;
    });
    if (SameDef != Pool.end()) {
      if (Score.rank() > SameDef->Score.rank())
        *SameDef = {&M, Score};
    } else if (Pool.size() < MaxPool) {
      Pool.push_back({&M, Score});
    }
  }
  const unsigned Limit = std::min(PlanWidth, RewriteBudget);
  if (Limit > 1)
    for (const auto &Pool : Pools)
      extendPool(Pool.second, Limit, P, Live, Best, BestScore);

  for (const Motion *M : Best) {
    const auto Anchor = P.Code[M->NewPoint]->getIterator();
    if (M->Address) {
      salvageDebugInfo(*M->Address);
      M->Address->moveBefore(Anchor);
      Moved.insert(M->Address);
    }
    salvageDebugInfo(*M->Def);
    M->Def->moveBefore(Anchor);
    Moved.insert(M->Def);
    if (M->Address || isa<LoadInst>(M->Def))
      ++Result.LoadsSunk;
    else
      ++Result.Sunk;
  }
  if (Best.size() > 1) {
    ++Result.MotionGroups;
    Result.GroupMembers += Best.size();
  }
  return !Best.empty();
}

} // namespace rssn
