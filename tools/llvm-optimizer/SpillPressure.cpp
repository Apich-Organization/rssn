#include "CallABI.h"
#include "LibmPressure.h"
#include "LibmReuse.h"
#include "LibmTailCutoff.h"
#include "LibmUnitLog.h"
#include "FPClassInterval.h"
#include "FPEndpointGuards.h"
#include "FPRegisterBank.h"
#include "LoopExitCompression.h"
#include "LoopPreheaderPressure.h"
#include "PowExponentPack.h"
#include "PowUnitBase.h"
#include "Planning.h"
#include "ScalarLoadPack.h"
#include "SwitchLoadReuse.h"
#include "llvm/ADT/DenseSet.h"
#include "llvm/Analysis/AliasAnalysis.h"
#include "llvm/Analysis/BlockFrequencyInfo.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/OptimizationRemarkEmitter.h"
#include "llvm/Analysis/TargetTransformInfo.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/Passes/PassBuilder.h"
#include "llvm/Plugins/PassPlugin.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Transforms/Utils/Local.h"

#include <cstdlib>

using namespace llvm;

namespace rssn {
namespace {

constexpr char PassName[] = "rssn-spill-pressure";
constexpr unsigned MaxLocalChanges = 64;
constexpr unsigned MaxBlockChanges = 16;
constexpr unsigned MaxCrossChanges = 3;
constexpr unsigned MaxCrossCandidates = 64;

enum class Mode { Off, Report, Transform };

} // namespace

bool modeled(Pressure P) { return P.GP || P.SIMD; }

bool cheapPure(const Instruction &I, const Liveness &Live,
               const TargetTransformInfo &TTI) {
  switch (I.getOpcode()) {
  case Instruction::Add:
  case Instruction::Sub:
  case Instruction::Mul:
  case Instruction::FAdd:
  case Instruction::FSub:
  case Instruction::FMul:
  case Instruction::And:
  case Instruction::Or:
  case Instruction::Xor:
  case Instruction::Shl:
  case Instruction::LShr:
  case Instruction::AShr:
  case Instruction::FNeg:
  case Instruction::GetElementPtr:
  case Instruction::Trunc:
  case Instruction::ZExt:
  case Instruction::SExt:
  case Instruction::BitCast:
  case Instruction::ExtractElement:
  case Instruction::InsertElement:
  case Instruction::ShuffleVector:
    break;
  default:
    return false;
  }
  if (!modeled(Live.weight(&I)) ||
      any_of(I.operands(), [&](const Use &U) {
        return !isa<Constant>(U) && !modeled(Live.weight(U.get()));
      }))
    return false;
  auto Cost = TTI.getInstructionCost(&I, TargetTransformInfo::TCK_RecipThroughput);
  return Cost.isValid() && Cost <= 1;
}

bool simpleLoad(const LoadInst &L, const Liveness &Live) {
  return L.isSimple() && L.getPointerAddressSpace() == 0 &&
         !L.hasMetadata(LLVMContext::MD_invariant_group) &&
         modeled(Live.weight(&L));
}

std::optional<AddressKey> accessAddress(const LoadInst &Load) {
  const DataLayout &DL = Load.getDataLayout();
  APInt Offset(DL.getIndexSizeInBits(Load.getPointerAddressSpace()), 0);
  const Value *Base = Load.getPointerOperand()
                          ->stripAndAccumulateInBoundsConstantOffsets(DL, Offset);
  if (!Offset.isSignedIntN(64))
    return std::nullopt;
  return AddressKey{Base, Offset.getSExtValue()};
}

namespace {

bool hasMatchingMetadata(const LoadInst &A, const LoadInst &B) {
  SmallVector<std::pair<unsigned, MDNode *>, 4> AM, BM;
  A.getAllMetadataOtherThanDebugLoc(AM);
  B.getAllMetadataOtherThanDebugLoc(BM);
  return AM == BM;
}

bool mergeFits(const LoadInst &First, const LoadInst &Second,
               const BlockPressure &P, const Liveness &Live,
               const TargetTransformInfo &TTI, unsigned To) {
  const unsigned Last = P.LastUse.lookup(&First);
  const Pressure Added = Live.weight(&Second);
  int Area = 0;
  for (unsigned K = Last + 1; K <= To; ++K) {
    Pressure New = P.Before[K];
    New += Added;
    if (New.worsens(P.Before[K], Live.Limit))
      return false;
    Area += New.congestion(Live.Limit) - P.Before[K].congestion(Live.Limit);
  }
  // A use at a call is live before it but not live through it.
  for (int K = P.Calls.find_next(Last ? Last - 1 : 0);
       K >= 0 && static_cast<unsigned>(K) < To; K = P.Calls.find_next(K)) {
    Pressure New = P.Through[K];
    New += Added;
    if (New.worsens(P.Through[K], Liveness::CallLimit))
      return false;
    Area += CallWeight * (New.callCost(Liveness::CallLimit) -
                         P.Through[K].callCost(Liveness::CallLimit));
  }
  // After the second definition, overlapping copies collapse to one value.
  const unsigned OverlapEnd =
      std::min({Last, P.LastUse.lookup(&Second), To + MaxDistance,
                static_cast<unsigned>(P.Code.size() - 1)});
  for (unsigned K = To + 1; K <= OverlapEnd; ++K) {
    Pressure New = P.Before[K];
    New -= Added;
    Area += New.congestion(Live.Limit) - P.Before[K].congestion(Live.Limit);
  }
  const auto Saved = TTI.getInstructionCost(&Second,
                                           TargetTransformInfo::TCK_RecipThroughput);
  return Saved.isValid() &&
         Area <= 8 * std::clamp<int64_t>(Saved.getValue(), 1, 16);
}

bool tryMergeLoads(BlockPressure &P, const Liveness &Live,
                   const TargetTransformInfo &TTI, MemoryModel &Memory, Changes &Result) {
  DenseMap<std::pair<AddressKey, Type *>, LoadInst *> Previous;
  for (auto [K, I] : enumerate(P.Code)) {
    auto *Later = dyn_cast<LoadInst>(I);
    if (!Later || !simpleLoad(*Later, Live))
      continue;
    const auto Address = accessAddress(*Later);
    if (!Address)
      continue;
    const auto Key = std::make_pair(*Address, Later->getType());
    LoadInst *Earlier = Previous.lookup(Key);
    Previous[Key] = Later;
    if (!Earlier || Later->use_empty() || !P.FirstUse.count(Earlier) ||
        P.FirstUse.lookup(Earlier) >= K ||
        K - P.Positions.lookup(Earlier) > MaxDistance)
      continue;
    // Keep the dataflow at the block boundary independent of deletions.
    if (any_of(Later->users(), [&](const User *U) {
          const auto *Use = dyn_cast<Instruction>(U);
          return !Use || Use->getParent() != Later->getParent() ||
                 isa<PHINode>(Use);
        }))
      continue;
    if (!hasMatchingMetadata(*Earlier, *Later) ||
        !mergeFits(*Earlier, *Later, P, Live, TTI, K) ||
        !Memory.unchanged(*Earlier, P.Code, P.Positions.lookup(Earlier) + 1, K))
      continue;
    Later->replaceAllUsesWith(Earlier);
    Later->eraseFromParent();
    ++Result.LoadsMerged;
    return true;
  }
  return false;
}

struct Destination {
  BasicBlock *Block;
  SmallVector<Use *, 4> Uses;
};

bool collectDestinations(Instruction &I, const LoopInfo &LI,
                         SmallVectorImpl<Destination> &Dests) {
  if (I.use_empty() || I.hasNUsesOrMore(13))
    return false;
  for (Use &U : I.uses()) {
    auto *User = dyn_cast<Instruction>(U.getUser());
    if (!User || isa<PHINode>(User))
      return false;
    BasicBlock *B = User->getParent();
    if (B == I.getParent() || B->isEHPad() ||
        LI.getLoopFor(I.getParent()) != LI.getLoopFor(B))
      return false;
    auto It = find_if(Dests, [B](const Destination &D) { return D.Block == B; });
    if (It == Dests.end()) {
      if (Dests.size() == 3)
        return false;
      Dests.push_back({B, {}});
      It = std::prev(Dests.end());
    }
    It->Uses.push_back(&U);
  }
  return true;
}

bool findInsertions(Instruction &I, ArrayRef<Destination> Dests,
                    const Liveness &Live, const BlockFrequencyInfo &BFI,
                    SmallVectorImpl<Instruction *> &Insertions, int &Gain) {
  const uint64_t SourceFreq = BFI.getBlockFreq(I.getParent()).getFrequency();
  uint64_t Remaining = SourceFreq;
  for (const Destination &D : Dests) {
    const uint64_t Frequency = BFI.getBlockFreq(D.Block).getFrequency();
    if (!Frequency || Frequency > Remaining)
      return false;
    Remaining -= Frequency;
    BlockPressure P(*D.Block, Live);
    const unsigned Point = P.FirstUse.lookup(&I);
    if (any_of(I.operands(), [&](const Use &U) {
          return !isa<Constant>(U) && !P.isLiveBefore(U.get(), Point);
        }))
      return false;
    const Pressure Result = Live.weight(&I);
    for (unsigned K = 0; K < Point; ++K) {
      Pressure After = P.Before[K];
      After -= Result;
      Gain += P.Before[K].congestion(Live.Limit) - After.congestion(Live.Limit);
      if (P.Calls.test(K)) {
        After = P.Through[K];
        After -= Result;
        Gain += CallWeight * (P.Through[K].callCost(Liveness::CallLimit) -
                              After.callCost(Liveness::CallLimit));
      }
    }
    Insertions.push_back(P.Code[Point]);
  }
  return true;
}

bool tryCrossBlock(Function &F, const Liveness &Live,
                   const TargetTransformInfo &TTI, const DominatorTree &DT,
                   const LoopInfo &LI, const BlockFrequencyInfo &BFI,
                   Changes &Result) {
  unsigned Candidates = 0;
  for (BasicBlock &B : F) {
    if (B.isEHPad() || !DT.isReachableFromEntry(&B))
      continue;
    BlockPressure P(B, Live);
    if (!P.Pressured)
      continue;
    for (unsigned Def = 0; Def < P.Code.size(); ++Def) {
      Instruction *I = P.Code[Def];
      if (P.FirstUse.count(I) || I->use_empty() || !cheapPure(*I, Live, TTI))
        continue;
      SmallVector<Destination, 3> Dests;
      if (!collectDestinations(*I, LI, Dests))
        continue;
      if (++Candidates > MaxCrossCandidates)
        return false;
      if (Dests.size() > 1 && any_of(I->operands(), [&](const Use &U) {
            return !isGuaranteedNotToBeUndefOrPoison(U.get(), nullptr, I, &DT);
          }))
        continue;

      const Pressure W = Live.weight(I);
      int Gain = 0;
      for (unsigned K = Def + 1; K < P.Code.size(); ++K) {
        Pressure After = P.Before[K];
        After -= W;
        Gain += P.Before[K].congestion(Live.Limit) - After.congestion(Live.Limit);
        if (P.Calls.test(K)) {
          After = P.Through[K];
          After -= W;
          Gain += CallWeight * (P.Through[K].callCost(Liveness::CallLimit) -
                                After.callCost(Liveness::CallLimit));
        }
      }
      SmallVector<Instruction *, 3> Insertions;
      const auto Size = TTI.getInstructionCost(I, TargetTransformInfo::TCK_CodeSize);
      const auto Rate = TTI.getInstructionCost(I, TargetTransformInfo::TCK_RecipThroughput);
      const auto ExtraCost = Size + Rate;
      const unsigned Copies = Dests.size() - 1;
      if ((Copies && (!ExtraCost.isValid() || ExtraCost > 4)) ||
          !findInsertions(*I, Dests, Live, BFI, Insertions, Gain) ||
          Gain < 6 + static_cast<int>(Copies) *
                         (ExtraCost.isValid() ? std::max<int64_t>(2, ExtraCost.getValue()) : 0))
        continue;

      salvageDebugInfo(*I);
      if (Dests.size() == 1) {
        I->moveBefore(Insertions.front()->getIterator());
      } else {
        for (auto [D, Insert] : zip_equal(Dests, Insertions)) {
          Instruction *Copy = I->clone();
          Copy->insertBefore(Insert->getIterator());
          Copy->setName(I->getName() + ".rssn");
          for (Use *U : D.Uses)
            U->set(Copy);
        }
        Result.Copies += Dests.size() - 1;
        I->eraseFromParent();
      }
      ++Result.CrossBlock;
      return true;
    }
  }
  return false;
}

bool selected(Function &F, StringRef Filter, bool ReportWin64) {
  if (!F.getName().contains(Filter) || F.isDeclaration() || F.hasOptNone() ||
      F.hasOptSize() || F.hasMinSize() || F.hasGC() ||
      F.hasFnAttribute(Attribute::StrictFP) ||
      F.hasFnAttribute(Attribute::Naked))
    return false;
  const Triple &T = F.getParent()->getTargetTriple();
  return T.getArch() == Triple::x86_64 &&
         (ReportWin64 ||
          (T.isOSLinux() && T.getEnvironment() != Triple::GNUX32));
}

void reportWin64Pressure(Function &F, FunctionAnalysisManager &AM,
                         const Liveness &Live, OptimizationRemarkEmitter &ORE) {
  struct CallSummary {
    const Instruction *Hot = nullptr;
    Pressure AtHot;
    double Score = 0;
    unsigned Sites = 0;
    unsigned OverBudget = 0;
    int Unknown = 0;
  };
  CallSummary Explicit, PossibleLibm;
  const auto &BFI = AM.getResult<BlockFrequencyAnalysis>(F);
  const auto &DT = AM.getResult<DominatorTreeAnalysis>(F);
  for (BasicBlock &B : F) {
    if (!DT.isReachableFromEntry(&B))
      continue;
    const BlockPressure P(B, Live, true, {}, CallPressureABI::Win64);
    const double Frequency = BFI.getBlockFreq(&B).getFrequency();
    for (unsigned Point : P.Calls.set_bits()) {
      const Instruction *I = P.Code[Point];
      auto &S = isa<IntrinsicInst>(I) ? PossibleLibm : Explicit;
      const Pressure Through = P.Through[Point];
      ++S.Sites;
      S.Unknown = std::max(S.Unknown, Through.SIMD - Through.PreservableSIMD -
                                          Through.VolatileSIMD);
      const int Excess = Through.win64CallExcess();
      S.OverBudget += Excess > 0;
      const double Score = Frequency * Excess;
      if (!S.Hot || Score > S.Score) {
        S.Hot = I;
        S.AtHot = Through;
        S.Score = Score;
      }
    }
  }
  const auto Emit = [&](const CallSummary &S, StringRef Kind) {
    if (!S.Hot)
      return;
    const auto *Callee = cast<CallBase>(S.Hot)->getCalledFunction();
    ORE.emit([&] {
      using NV = ore::NV;
      std::string ResultName;
      raw_string_ostream ResultOS(ResultName);
      if (S.Hot->getType()->isVoidTy())
        ResultOS << "<void>";
      else
        S.Hot->printAsOperand(ResultOS, false, F.getParent());
      return OptimizationRemarkAnalysis(PassName, "Win64CallPressure", S.Hot)
             << NV("Function", F.getName()) << " Win64 " << NV("Kind", Kind)
             << " sites " << NV("Sites", S.Sites) << ", over preservation budget "
             << NV("OverBudget", S.OverBudget) << ", at selected site GP "
             << NV("GP", S.AtHot.GP) << " / " << NV("PreservedGP", Pressure::Win64PreservedGP)
             << ", XMM-sized " << NV("XMM", S.AtHot.PreservableSIMD) << " / "
             << NV("PreservedXMM", Pressure::Win64PreservedXMM) << ", native wide SIMD "
             << NV("NativeWide", S.AtHot.VolatileSIMD) << ", unclassified SIMD peak "
             << NV("Unclassified", S.Unknown) << ", callee "
             << NV("Callee", Callee ? Callee->getName() : StringRef("<indirect>"))
             << ", IR result " << NV("IRResult", ResultName);
    });
  };
  Emit(Explicit, "ABI call");
  Emit(PossibleLibm, "possible libm");
}

void reportWin64ArgumentSlots(Function &F, FunctionAnalysisManager &AM,
                              OptimizationRemarkEmitter &ORE) {
  if (F.getCallingConv() == CallingConv::C ||
      F.getCallingConv() == CallingConv::Fast ||
      F.getCallingConv() == CallingConv::Win64) {
    const auto Slots = win64ScalarStackSlots(*F.getFunctionType(), F.getAttributes());
    if (Slots && *Slots)
      ORE.emit([&] {
        using NV = ore::NV;
        const unsigned Used = count_if(drop_begin(F.args(), 4), [](const Argument &A) {
          return !A.use_empty();
        });
        return OptimizationRemarkAnalysis(PassName, "Win64IncomingStackSlots", &F)
               << NV("Function", F.getName()) << " Win64 incoming scalar stack slots "
               << NV("Slots", *Slots) << ", used parameters " << NV("Used", Used);
      });
  }
  const auto &BFI = AM.getResult<BlockFrequencyAnalysis>(F);
  const auto &DT = AM.getResult<DominatorTreeAnalysis>(F);
  const CallInst *HotStackCall = nullptr;
  unsigned HotStackSlots = 0;
  unsigned StackCallSites = 0;
  double StackScore = 0;
  for (const BasicBlock &B : F) {
    if (!DT.isReachableFromEntry(&B))
      continue;
    const double Frequency = BFI.getBlockFreq(&B).getFrequency();
    for (const Instruction &I : B) {
      if (!isWin64CallPressurePoint(I))
        continue;
      const auto &Call = cast<CallInst>(I);
      if (Call.arg_size() <= 4)
        continue;
      const auto Slots = win64ScalarStackSlots(*Call.getFunctionType(),
                                               Call.getAttributes());
      if (!Slots || *Slots == 0)
        continue;
      ++StackCallSites;
      const double Score = Frequency * *Slots;
      if (!HotStackCall || Score > StackScore) {
        HotStackCall = &Call;
        HotStackSlots = *Slots;
        StackScore = Score;
      }
    }
  }
  if (HotStackCall)
    ORE.emit([&] {
      using NV = ore::NV;
      const Function *Callee = HotStackCall->getCalledFunction();
      return OptimizationRemarkAnalysis(PassName, "Win64OutgoingStackSlots", HotStackCall)
             << NV("Function", F.getName()) << " Win64 calls with scalar stack slots "
             << NV("Sites", StackCallSites) << ", at selected site slots "
             << NV("Slots", HotStackSlots) << ", callee "
             << NV("Callee", Callee ? Callee->getName() : StringRef("<indirect>"));
    });
}

class SpillPressurePass : public PassInfoMixin<SpillPressurePass> {
  Mode RunMode;
  std::string Filter;
  unsigned PlanWidth;

public:
  SpillPressurePass(Mode RunMode, StringRef Filter, unsigned PlanWidth)
      : RunMode(RunMode), Filter(Filter.str()), PlanWidth(PlanWidth) {}

  PreservedAnalyses run(Function &F, FunctionAnalysisManager &AM) {
    const bool ReportWin64 = RunMode == Mode::Report &&
                             F.getParent()->getTargetTriple().isOSWindows();
    if (!selected(F, Filter, ReportWin64))
      return PreservedAnalyses::all();

    auto &ORE = AM.getResult<OptimizationRemarkEmitterAnalysis>(F);
    if (ReportWin64) {
      if (!ORE.allowExtraAnalysis(PassName))
        return PreservedAnalyses::all();
      reportWin64ArgumentSlots(F, AM, ORE);
    }
    if (!Liveness::supports(F))
      return PreservedAnalyses::all();
    auto &TTI = AM.getResult<TargetIRAnalysis>(F);
    auto Live = std::make_unique<Liveness>(TTI, F);
    if (!Live->compute(F))
      return PreservedAnalyses::all();
    if (ReportWin64) {
      reportWin64Pressure(F, AM, *Live, ORE);
      return PreservedAnalyses::all();
    }
    const bool WantRemarks = ORE.allowExtraAnalysis(PassName) ||
                             RunMode == Mode::Report;
    Summary Before;
    if (WantRemarks)
      Before = summarize(F, *Live);

    Changes Result;
    unsigned PrivateBoundaries = 0;
    if (RunMode == Mode::Transform) {
      MemoryModel Memory(AM.getResult<AAManager>(F));
      // Fusion introduces a new vector SSA definition. Refresh the table
      // once after the bounded prepass, before any interval-based planning.
      for (BasicBlock &B : F) {
        if (B.isEHPad())
          continue;
        if (Result.PacksFused + Result.LaneGroupsFused >= 24)
          break;
        if (!tryFusePackedLoads(B, *Live, TTI, Memory, Result))
          tryFuseLaneLoads(B, *Live, TTI, Memory, Result);
      }
      if (Result.PacksFused || Result.LaneGroupsFused) {
        Live = std::make_unique<Liveness>(TTI, F);
        if (!Live->compute(F))
          return PreservedAnalyses::none();
      }

      unsigned Total = 0;
      for (BasicBlock &B : F) {
        if (B.isEHPad() || B.size() < 4)
          continue;
        DenseSet<const Instruction *> Moved;
        unsigned BlockChanges = 0;
        while (BlockChanges < MaxBlockChanges && Total < MaxLocalChanges) {
          const unsigned Previous = Result.Sunk + Result.LoadsSunk +
                                    Result.LoadsMerged;
          const unsigned Budget = std::min(MaxBlockChanges - BlockChanges,
                                           MaxLocalChanges - Total);
          BlockPressure P(B, *Live);
          if (!tryMergeLoads(P, *Live, TTI, Memory, Result) &&
              !trySink(P, *Live, TTI, Memory, Moved, Result, Budget, PlanWidth))
            break;
          const unsigned Applied = Result.Sunk + Result.LoadsSunk +
                                   Result.LoadsMerged - Previous;
          BlockChanges += Applied;
          Total += Applied;
        }
        if (Total == MaxLocalChanges)
          break;
      }
      PrivateBoundaries = Memory.privateBoundaries();
      // Local changes preserve all live-out sets. Rebuild the value table
      // before global motion so erased values no longer have IDs.
      if (Result.LoadsMerged) {
        Live = std::make_unique<Liveness>(TTI, F);
        if (!Live->compute(F))
          return PreservedAnalyses::none();
      }
      if (F.size() > 1) {
        auto &DT = AM.getResult<DominatorTreeAnalysis>(F);
        auto &LI = AM.getResult<LoopAnalysis>(F);
        if (any_of(F, [&](BasicBlock &B) { return !Live->liveOut(B).none(); })) {
          auto &BFI = AM.getResult<BlockFrequencyAnalysis>(F);
          for (unsigned K = 0; K < MaxCrossChanges; ++K) {
            if (!tryCrossBlock(F, *Live, TTI, DT, LI, BFI, Result))
              break;
            Live = std::make_unique<Liveness>(TTI, F);
            if (!Live->compute(F))
              return PreservedAnalyses::none();
          }
        }
      }
    }

    if (WantRemarks) {
      const Summary After = Result.any() ? summarize(F, *Live) : Before;
      ORE.emit([&] {
        using NV = ore::NV;
        return OptimizationRemarkAnalysis(PassName, "Pressure", &F)
               << NV("Function", F.getName()) << " GP peak "
               << NV("GPPeakBefore", Before.Peak.GP) << " -> "
               << NV("GPPeakAfter", After.Peak.GP) << ", SIMD peak "
               << NV("SIMDPeakBefore", Before.Peak.SIMD) << " -> "
               << NV("SIMDPeakAfter", After.Peak.SIMD) << ", over-limit points "
               << NV("ExcessBefore", Before.ExcessPoints) << " -> "
               << NV("ExcessAfter", After.ExcessPoints) << ", SIMD live through "
               << NV("CallSites", Before.CallSites) << " call sites "
               << NV("SIMDCallsBefore", Before.SIMDThroughCalls) << " -> "
               << NV("SIMDCallsAfter", After.SIMDThroughCalls) << ", sunk "
               << NV("Sunk", Result.Sunk) << ", delayed loads "
               << NV("LoadsSunk", Result.LoadsSunk) << ", merged loads "
               << NV("LoadsMerged", Result.LoadsMerged) << ", cross-block "
               << NV("CrossBlock", Result.CrossBlock) << ", extra copies "
               << NV("Copies", Result.Copies) << ", congestion "
               << NV("CongestionBefore", Before.Congestion) << " -> "
               << NV("CongestionAfter", After.Congestion) << ", SIMD call cost "
               << NV("SIMDCallCostBefore", Before.SIMDSpillThroughCalls) << " -> "
               << NV("SIMDCallCostAfter", After.SIMDSpillThroughCalls)
               << ", joint plans " << NV("MotionGroups", Result.MotionGroups)
               << ", joint members " << NV("GroupMembers", Result.GroupMembers)
               << ", fused packs " << NV("PacksFused", Result.PacksFused)
               << ", lane groups " << NV("LaneGroupsFused", Result.LaneGroupsFused)
               << ", source loads fused " << NV("ReadsFused", Result.ReadsFused)
               << ", private call boundaries "
               << NV("PrivateBoundaries", PrivateBoundaries)
               << ", GP live through calls "
               << NV("GPThroughBefore", Before.GPThroughCalls) << " -> "
               << NV("GPThroughAfter", After.GPThroughCalls)
               << ", GP call excess "
               << NV("GPCallExcessBefore", Before.GPCallExcess) << " -> "
               << NV("GPCallExcessAfter", After.GPCallExcess);
      });
    }
    if (!Result.any())
      return PreservedAnalyses::all();
    PreservedAnalyses PA;
    PA.preserveSet<CFGAnalyses>();
    return PA;
  }
};

StringRef env(const char *Name) {
  const char *Value = std::getenv(Name);
  return Value ? StringRef(Value) : StringRef();
}

Mode automaticMode() {
  const StringRef Setting = env("RSSN_SPILL_MODE");
  if (Setting.empty() || Setting == "transform")
    return Mode::Transform;
  if (Setting == "report")
    return Mode::Report;
  if (Setting == "off")
    return Mode::Off;
  report_fatal_error("RSSN_SPILL_MODE must be transform, report, or off");
}

Mode abiMode() {
  const StringRef Setting = env("RSSN_SPILL_ABI");
  if (Setting.empty() || Setting == "off")
    return Mode::Off;
  if (Setting == "report")
    return Mode::Report;
  if (Setting == "transform")
    return Mode::Transform;
  report_fatal_error("RSSN_SPILL_ABI must be off, report, or transform");
}

Mode libmMode(const char *Name) {
  const StringRef Setting = env(Name);
  if (Setting.empty() || Setting == "off")
    return Mode::Off;
  if (Setting == "report")
    return Mode::Report;
  if (Setting == "transform")
    return Mode::Transform;
  report_fatal_error(Twine(Name) + " must be off, report, or transform");
}

unsigned motionPlanWidth() {
  const StringRef Setting = env("RSSN_SPILL_PLAN_WIDTH");
  unsigned Width = MaxPlanWidth;
  if (!Setting.empty() &&
      (Setting.getAsInteger(10, Width) || !Width || Width > MaxPlanWidth))
    report_fatal_error("RSSN_SPILL_PLAN_WIDTH must be between 1 and 4");
  return Width;
}

void registerPasses(PassBuilder &PB) {
  const std::string Filter = env("RSSN_SPILL_FILTER").str();
  const unsigned Width = motionPlanWidth();
  PB.registerPipelineParsingCallback(
      [Filter](StringRef Name, ModulePassManager &MPM,
               ArrayRef<PassBuilder::PipelineElement> Inner) {
        if (!Inner.empty() || (Name != "rssn-call-abi" &&
                               Name != "rssn-call-abi-report"))
          return false;
        MPM.addPass(CallABIPass(Name == "rssn-call-abi", Filter));
        return true;
      });
  PB.registerPipelineParsingCallback(
      [Filter, Width](StringRef Name, FunctionPassManager &FPM,
               ArrayRef<PassBuilder::PipelineElement> Inner) {
        if (!Inner.empty())
          return false;
        if (Name == PassName)
          FPM.addPass(SpillPressurePass(Mode::Transform, Filter, Width));
        else if (Name == "rssn-spill-pressure-report")
          FPM.addPass(SpillPressurePass(Mode::Report, Filter, Width));
        else if (Name == "rssn-libm-pressure" || Name == "rssn-libm-pressure-report" ||
                 Name == "rssn-libm-cone-trial")
          FPM.addPass(LibmPressurePass(Name != "rssn-libm-pressure-report", Filter,
                                       Name == "rssn-libm-cone-trial"));
        else if (Name == "rssn-libm-reuse" || Name == "rssn-libm-reuse-report")
          FPM.addPass(LibmReusePass(Name == "rssn-libm-reuse", Filter));
        else if (Name == "rssn-libm-tail-cutoff" || Name == "rssn-libm-tail-cutoff-report")
          FPM.addPass(LibmTailCutoffPass(Name == "rssn-libm-tail-cutoff", Filter));
        else if (Name == "rssn-libm-unit-log" || Name == "rssn-libm-unit-log-report")
          FPM.addPass(LibmUnitLogPass(Name == "rssn-libm-unit-log", Filter));
        else if (Name == "rssn-pow-exponent-pack" || Name == "rssn-pow-exponent-pack-report")
          FPM.addPass(PowExponentPackPass(Name == "rssn-pow-exponent-pack", Filter));
        else if (Name == "rssn-loop-exit-compression" ||
                 Name == "rssn-loop-exit-compression-report")
          FPM.addPass(LoopExitCompressionPass(Name == "rssn-loop-exit-compression", Filter));
        else if (Name == "rssn-loop-preheader-pressure" ||
                 Name == "rssn-loop-preheader-pressure-report")
          FPM.addPass(LoopPreheaderPressurePass(
              Name == "rssn-loop-preheader-pressure", Filter));
        else if (Name == "rssn-fp-endpoint-guards" ||
                 Name == "rssn-fp-endpoint-guards-report")
          FPM.addPass(FPEndpointGuardsPass(Name == "rssn-fp-endpoint-guards", Filter));
        else if (Name == "rssn-fp-class-interval" || Name == "rssn-fp-class-interval-report")
          FPM.addPass(FPClassIntervalPass(Name == "rssn-fp-class-interval", Filter));
        else if (Name == "rssn-scalar-load-pack" || Name == "rssn-scalar-load-pack-report")
          FPM.addPass(ScalarLoadPackPass(Name == "rssn-scalar-load-pack", Filter));
        else if (Name == "rssn-switch-load-reuse" ||
                 Name == "rssn-switch-load-reuse-report")
          FPM.addPass(SwitchLoadReusePass(Name == "rssn-switch-load-reuse", Filter));
        else if (Name == "rssn-fp-register-bank" ||
                 Name == "rssn-fp-register-bank-report")
          FPM.addPass(FPRegisterBankPass(Name == "rssn-fp-register-bank", Filter));
        else if (Name == "rssn-pow-unit-base" || Name == "rssn-pow-unit-base-report")
          FPM.addPass(PowUnitBasePass(Name == "rssn-pow-unit-base", Filter));
        else
          return false;
        return true;
      });

  const Mode Auto = automaticMode();
  const Mode ABI = abiMode();
  const Mode Libm = libmMode("RSSN_SPILL_LIBM");
  const Mode LibmReuse = libmMode("RSSN_SPILL_LIBM_REUSE");
  const Mode LibmTailCutoff = libmMode("RSSN_LIBM_TAIL_CUTOFF");
  const Mode LibmUnitLog = libmMode("RSSN_LIBM_UNIT_LOG");
  const Mode ExponentPack = libmMode("RSSN_POW_EXPONENT_PACK");
  const Mode ExitCompression = libmMode("RSSN_LOOP_EXIT_COMPRESSION");
  const Mode PreheaderPressure = libmMode("RSSN_LOOP_PREHEADER_PRESSURE");
  const Mode EndpointGuards = libmMode("RSSN_FP_ENDPOINT_GUARDS");
  const Mode ClassInterval = libmMode("RSSN_FP_CLASS_INTERVAL");
  const Mode LoadPack = libmMode("RSSN_SCALAR_LOAD_PACK");
  const Mode SwitchReuse = libmMode("RSSN_SWITCH_LOAD_REUSE");
  const Mode FPBank = libmMode("RSSN_FP_REGISTER_BANK");
  const Mode PowUnit = libmMode("RSSN_POW_UNIT_BASE");
  auto Add = [Auto, ABI, Libm, LibmReuse, LibmTailCutoff, LibmUnitLog,
              ExponentPack, ExitCompression, PreheaderPressure, EndpointGuards,
              ClassInterval, LoadPack, SwitchReuse, FPBank, PowUnit, Filter, Width](
                 ModulePassManager &MPM, OptimizationLevel Level) {
    if (Level != OptimizationLevel::O3)
      return;
    if (ABI != Mode::Off)
      MPM.addPass(CallABIPass(ABI == Mode::Transform, Filter));
    if (Libm != Mode::Off)
      MPM.addPass(createModuleToFunctionPassAdaptor(
          LibmPressurePass(Libm == Mode::Transform, Filter)));
    if (LibmReuse != Mode::Off)
      MPM.addPass(createModuleToFunctionPassAdaptor(
          LibmReusePass(LibmReuse == Mode::Transform, Filter)));
    if (LibmTailCutoff != Mode::Off)
      MPM.addPass(createModuleToFunctionPassAdaptor(
          LibmTailCutoffPass(LibmTailCutoff == Mode::Transform, Filter)));
    if (LibmUnitLog != Mode::Off)
      MPM.addPass(createModuleToFunctionPassAdaptor(
          LibmUnitLogPass(LibmUnitLog == Mode::Transform, Filter)));
    if (EndpointGuards != Mode::Off)
      MPM.addPass(createModuleToFunctionPassAdaptor(
          FPEndpointGuardsPass(EndpointGuards == Mode::Transform, Filter)));
    if (ClassInterval != Mode::Off)
      MPM.addPass(createModuleToFunctionPassAdaptor(
          FPClassIntervalPass(ClassInterval == Mode::Transform, Filter)));
    if (LoadPack != Mode::Off)
      MPM.addPass(createModuleToFunctionPassAdaptor(
          ScalarLoadPackPass(LoadPack == Mode::Transform, Filter)));
    if (SwitchReuse != Mode::Off)
      MPM.addPass(createModuleToFunctionPassAdaptor(
          SwitchLoadReusePass(SwitchReuse == Mode::Transform, Filter)));
    if (PowUnit != Mode::Off)
      MPM.addPass(createModuleToFunctionPassAdaptor(
          PowUnitBasePass(PowUnit == Mode::Transform, Filter)));
    if (ExitCompression != Mode::Off)
      MPM.addPass(createModuleToFunctionPassAdaptor(
          LoopExitCompressionPass(ExitCompression == Mode::Transform, Filter)));
    if (PreheaderPressure != Mode::Off)
      MPM.addPass(createModuleToFunctionPassAdaptor(
          LoopPreheaderPressurePass(PreheaderPressure == Mode::Transform,
                                    Filter)));
    if (Auto != Mode::Off)
      MPM.addPass(
          createModuleToFunctionPassAdaptor(SpillPressurePass(Auto, Filter, Width)));
    if (FPBank != Mode::Off)
      MPM.addPass(createModuleToFunctionPassAdaptor(
          FPRegisterBankPass(FPBank == Mode::Transform, Filter)));
    if (ExponentPack != Mode::Off)
      MPM.addPass(createModuleToFunctionPassAdaptor(
          PowExponentPackPass(ExponentPack == Mode::Transform, Filter)));
  };
  PB.registerPipelineParsingCallback(
      [Add](StringRef Name, ModulePassManager &MPM,
            ArrayRef<PassBuilder::PipelineElement> Inner) {
        if (Name != "rssn-release" || !Inner.empty())
          return false;
        Add(MPM, OptimizationLevel::O3);
        return true;
      });
  if (Auto == Mode::Off && ABI == Mode::Off && Libm == Mode::Off &&
      LibmReuse == Mode::Off && LibmTailCutoff == Mode::Off &&
      LibmUnitLog == Mode::Off && ExponentPack == Mode::Off &&
      ExitCompression == Mode::Off && PreheaderPressure == Mode::Off &&
      EndpointGuards == Mode::Off && ClassInterval == Mode::Off &&
      LoadPack == Mode::Off && SwitchReuse == Mode::Off && FPBank == Mode::Off &&
      PowUnit == Mode::Off)
    return;
  PB.registerOptimizerLastEPCallback(
      [Add](ModulePassManager &MPM, OptimizationLevel Level,
            ThinOrFullLTOPhase Phase) {
        if (Phase == ThinOrFullLTOPhase::None ||
            Phase == ThinOrFullLTOPhase::ThinLTOPostLink)
          Add(MPM, Level);
      });
  PB.registerFullLinkTimeOptimizationLastEPCallback(Add);
}

} // namespace
} // namespace rssn

#if defined(_WIN32)
extern "C" __declspec(dllexport) PassPluginLibraryInfo
#else
extern "C" LLVM_ATTRIBUTE_WEAK LLVM_ATTRIBUTE_VISIBILITY_DEFAULT PassPluginLibraryInfo
#endif
llvmGetPassPluginInfo() {
  return {LLVM_PLUGIN_API_VERSION, "RssnLLVMOptimizer", "0.1",
          rssn::registerPasses, nullptr};
}
