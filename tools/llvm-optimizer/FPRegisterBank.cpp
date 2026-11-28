#include "FPRegisterBank.h"
#include "LibmPressure.h"
#include "Pressure.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/Analysis/BlockFrequencyInfo.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/OptimizationRemarkEmitter.h"
#include "llvm/Analysis/TargetTransformInfo.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/InlineAsm.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/Module.h"
#include <limits>
#include <optional>

using namespace llvm;

namespace rssn {
namespace {

constexpr char PassName[] = "rssn-fp-register-bank";
constexpr int PreservedXMM = 10;
// Win64 has eight nonvolatile GPRs other than RSP. Leave two for the frame
// pointer, address materialization, and integer values introduced in lowering.
constexpr int GPRBudget = 6;
constexpr unsigned MaxCandidates = 24;
constexpr unsigned MaxMoves = 4;

struct Block {
  BlockPressure P;
  unsigned Offset;
  double Frequency;
};

struct Site {
  Instruction *Consumer;
  unsigned Cut;
  double Frequency;
};

struct Plan {
  Value *Scalar;
  Instruction *Insert;
  double SetupCost;
  SmallVector<Site, 6> Sites;
  BitVector AsInteger;
  BitVector AcrossCall;
};

bool supportedCut(const Instruction &I) {
  const auto *C = dyn_cast<CallBase>(&I);
  if (!C || isa<IntrinsicInst>(C))
    return true;
  return isa<CallInst>(C) && !C->isInlineAsm() && !C->isMustTailCall() &&
         !C->hasFnAttr(Attribute::ReturnsTwice) &&
         (C->getCallingConv() == CallingConv::C ||
          C->getCallingConv() == CallingConv::Fast);
}

bool reachesLibm(const Instruction &Start) {
  const Instruction *I = &Start;
  for (unsigned Step = 0; Step != 3; ++Step) {
    if (isScalarLibmCut(*I))
      return true;
    if (!I->hasOneUse() || (I->getOpcode() != Instruction::FMul &&
                            I->getOpcode() != Instruction::FSub &&
                            I->getOpcode() != Instruction::FAdd &&
                            I->getOpcode() != Instruction::FNeg))
      return false;
    const auto *Next = dyn_cast<Instruction>(*I->user_begin());
    if (!Next || Next->getParent() != I->getParent())
      return false;
    I = Next;
  }
  return false;
}

bool stackArgument(const Argument &A) {
  const Function &F = *A.getParent();
  if (A.getArgNo() < 4 || F.getCallingConv() != CallingConv::C ||
      (!F.getReturnType()->isVoidTy() && !F.getReturnType()->isDoubleTy()))
    return false;
  return win64ScalarStackSlots(*F.getFunctionType(), F.getAttributes())
      .has_value();
}

Instruction *placement(Value &V, ArrayRef<Instruction *> Users,
                       const DominatorTree &DT, const LoopInfo &LI,
                       const DenseMap<const BasicBlock *, Block> &Blocks) {
  BasicBlock *Common = Users.front()->getParent();
  for (Instruction *U : Users) {
    Common = DT.findNearestCommonDominator(Common, U->getParent());
    if (!Common)
      return nullptr;
  }
  while (const Loop *L = LI.getLoopFor(Common)) {
    Common = L->getLoopPreheader();
    if (!Common)
      return nullptr;
  }
  // Delay initialization until the first common path to the consumers.
  // Keep uses before the first call in the defining block in the FP bank.
  Instruction *Insert = Common->getTerminator();
  const auto *Def = dyn_cast<Instruction>(&V);
  if (Def && Def->getParent() == Common) {
    const BlockPressure &P = Blocks.find(Common)->second.P;
    for (unsigned Cut : P.Calls.set_bits())
      if (Cut > P.Positions.lookup(Def) && P.isLiveBefore(&V, Cut + 1))
        return P.Code[Cut];
  } else {
    for (Instruction *U : Users)
      if (U->getParent() == Common && U->comesBefore(Insert))
        Insert = U;
  }
  return Insert;
}

std::optional<Plan> plan(Value &V, const DominatorTree &DT, const LoopInfo &LI,
                         const DenseMap<const BasicBlock *, Block> &Blocks,
                         unsigned NumPoints) {
  if (!V.getType()->isDoubleTy() || V.getNumUses() < 2 || V.getNumUses() > 16)
    return std::nullopt;
  const auto *A = dyn_cast<Argument>(&V);
  const auto *Def = dyn_cast<Instruction>(&V);
  constexpr unsigned ComputedOps[] = {Instruction::FAdd, Instruction::FSub,
                                      Instruction::FMul, Instruction::FDiv,
                                      Instruction::FNeg};
  if (A ? !stackArgument(*A)
        : !Def || !is_contained(ComputedOps, Def->getOpcode()))
    return std::nullopt;

  SmallVector<Instruction *, 8> Users;
  SmallPtrSet<Instruction *, 8> Seen;
  for (User *U : V.users()) {
    auto *I = dyn_cast<Instruction>(U);
    if (!I || I->isTerminator() || isa<PHINode>(I) ||
        I->getParent()->isEHPad() || !DT.isReachableFromEntry(I->getParent()))
      return std::nullopt;
    if (Seen.insert(I).second)
      Users.push_back(I);
  }
  Instruction *Insert = placement(V, Users, DT, LI, Blocks);
  if (!Insert || Insert->getParent()->isEHPad() || !DT.dominates(&V, Insert))
    return std::nullopt;
  const auto &Home = Blocks.find(Insert->getParent())->second;
  // One setup transfer (or an integer load of an incoming stack argument)
  // and one save/restore at entry and exit if a new GPR is used.
  Plan P{&V,
         Insert,
         Home.Frequency + 2.0,
         {},
         BitVector(NumPoints),
         BitVector(NumPoints)};
  SmallPtrSet<const BasicBlock *, 4> Hot;
  constexpr unsigned NoCut = std::numeric_limits<unsigned>::max();
  for (Instruction *I : Users) {
    if (I == Insert ? Def != nullptr : !DT.dominates(Insert, I))
      continue;
    BasicBlock *B = I->getParent();
    const auto &Info = Blocks.find(B)->second;
    const auto &BP = Info.P;
    const unsigned Position = BP.Positions.lookup(I);
    unsigned Anchor = NoCut;
    const Loop *L = LI.getLoopFor(B);
    if (L && !L->contains(Insert->getParent()) && reachesLibm(*I))
      for (unsigned Cut : BP.Calls.set_bits()) {
        if (Cut > Position)
          break;
        if (BP.isLiveBefore(&V, Cut + 1)) {
          Anchor = Info.Offset + Cut;
          Hot.insert(B);
        }
      }
    P.Sites.push_back({I, Anchor, Info.Frequency});
  }
  if (Hot.size() < 2)
    return std::nullopt;

  for (const auto &[B, Info] : Blocks) {
    if (!DT.dominates(Insert->getParent(), B))
      continue;
    const auto &BP = Info.P;
    const unsigned Start =
        B == Insert->getParent() ? BP.Positions.lookup(Insert) : 0;
    for (unsigned I = Start, End = BP.Code.size(); I != End; ++I) {
      if (BP.isLiveBefore(&V, I))
        P.AsInteger.set(Info.Offset + I);
      if (BP.Calls.test(I) && BP.isLiveBefore(&V, I + 1))
        P.AcrossCall.set(Info.Offset + I);
    }
  }
  return P;
}

double evaluate(const Plan &P, ArrayRef<Pressure> Before,
                ArrayRef<Pressure> Across, ArrayRef<unsigned> GP,
                ArrayRef<unsigned> Moved, int LocalBudget) {
  for (unsigned I : P.AcrossCall.set_bits())
    if (Across[I].GP + static_cast<int>(Moved[I]) >= GPRBudget)
      return 0;
  for (unsigned I : P.AsInteger.set_bits())
    if (Before[I].GP + static_cast<int>(GP[I]) >= LocalBudget)
      return 0;

  double Score = -P.SetupCost;
  for (const Site &S : P.Sites) {
    Score -= S.Frequency; // VMOVQ at each rewritten consumer.
    if (S.Cut < Across.size() && P.AcrossCall.test(S.Cut) &&
        Across[S.Cut].PreservableSIMD - static_cast<int>(Moved[S.Cut]) >
            PreservedXMM)
      // The estimate distinguishes a reload into a libm argument register
      // from a source that may fold into an ordinary FP operation.
      Score += S.Frequency * (isScalarLibmCut(*S.Consumer) ? 3 : 2);
  }
  return Score;
}

void moveToIntegerBank(const Plan &P) {
  IRBuilder<> Entry(P.Insert);
  Value *Scalar = P.Scalar;
  if (!isGuaranteedNotToBeUndefOrPoison(Scalar))
    Scalar = Entry.CreateFreeze(Scalar, "fp.bank.defined");
  Value *Bits = Entry.CreateBitCast(Scalar, Entry.getInt64Ty(), "fp.bank.bits");
  auto *FT =
      FunctionType::get(Entry.getDoubleTy(), {Entry.getInt64Ty()}, false);
  // A bare round-trip bitcast is eliminated before register allocation.
  // This input constraint preserves the integer interval up to its consumer.
  InlineAsm *Restore = InlineAsm::get(FT, "vmovq $1, $0", "=x,r", false);
  for (const Site &S : P.Sites) {
    IRBuilder<> At(S.Consumer);
    auto *FP = At.CreateCall(Restore, {Bits}, "fp.bank.restore");
    FP->setDoesNotAccessMemory();
    FP->setDoesNotThrow();
    FP->addFnAttr(Attribute::WillReturn);
    S.Consumer->replaceUsesOfWith(P.Scalar, FP);
  }
}

} // namespace

PreservedAnalyses FPRegisterBankPass::run(Function &F,
                                          FunctionAnalysisManager &AM) {
  const Triple &T = F.getParent()->getTargetTriple();
  if (T.getArch() != Triple::x86_64 || !T.isWindowsMSVCEnvironment() ||
      !F.getName().contains(Filter) || F.isDeclaration() || F.hasOptNone() ||
      F.hasOptSize() || F.hasMinSize() ||
      F.hasFnAttribute(Attribute::StrictFP) ||
      F.hasFnAttribute(Attribute::NoImplicitFloat) || !Liveness::supports(F))
    return PreservedAnalyses::all();
  const auto &TTI = AM.getResult<TargetIRAnalysis>(F);
  if (TTI.getRegisterBitWidth(TargetTransformInfo::RGK_FixedWidthVector)
          .getKnownMinValue() < 256)
    return PreservedAnalyses::all();
  const auto &LI = AM.getResult<LoopAnalysis>(F);
  if (LI.empty() ||
      llvm::count_if(F,
                     [&](const BasicBlock &B) {
                       return LI.getLoopFor(&B) && any_of(B, isScalarLibmCut);
                     }) < 2 ||
      any_of(F, [](const BasicBlock &B) {
        return B.isEHPad() || !all_of(B, supportedCut);
      }))
    return PreservedAnalyses::all();
  Liveness Live(TTI, F);
  if (!Live.compute(F))
    return PreservedAnalyses::all();

  DenseMap<const BasicBlock *, Block> Blocks;
  SmallVector<Pressure, 128> Before, Across;
  const auto &BFI = AM.getResult<BlockFrequencyAnalysis>(F);
  const double Entry = std::max<uint64_t>(1, BFI.getEntryFreq().getFrequency());
  for (BasicBlock &B : F) {
    auto &Info =
        Blocks
            .try_emplace(
                &B, Block{BlockPressure(B, Live, /*IncludeScalarLibm=*/true),
                          static_cast<unsigned>(Before.size()),
                          BFI.getBlockFreq(&B).getFrequency() / Entry})
            .first->second;
    append_range(Before, Info.P.Before);
    append_range(Across, Info.P.Through);
  }

  const auto &DT = AM.getResult<DominatorTreeAnalysis>(F);
  SmallVector<Plan, 8> Plans;
  auto Consider = [&](Value &V) {
    if (Plans.size() <= MaxCandidates)
      if (auto P = plan(V, DT, LI, Blocks, Before.size()))
        Plans.push_back(std::move(*P));
  };
  for (Argument &A : F.args())
    Consider(A);
  for (BasicBlock &B : F)
    if (!LI.getLoopFor(&B))
      for (Instruction &I : B)
        Consider(I);
  if (Plans.empty() || Plans.size() > MaxCandidates)
    return PreservedAnalyses::all();

  SmallVector<unsigned, 128> GP(Before.size()), Moved(Before.size());
  SmallVector<unsigned, MaxMoves> Chosen;
  auto &ORE = AM.getResult<OptimizationRemarkEmitterAnalysis>(F);
  while (Chosen.size() != MaxMoves) {
    double BestScore = 0;
    unsigned Best = Plans.size();
    for (auto [Index, P] : enumerate(Plans)) {
      if (is_contained(Chosen, Index))
        continue;
      const double Score =
          evaluate(P, Before, Across, GP, Moved, Live.Limit.GP);
      if (Score > BestScore) {
        BestScore = Score;
        Best = Index;
      }
    }
    if (Best == Plans.size())
      break;
    const Plan &P = Plans[Best];
    Chosen.push_back(Best);
    for (unsigned I : P.AsInteger.set_bits())
      ++GP[I];
    for (unsigned I : P.AcrossCall.set_bits())
      ++Moved[I];
    ORE.emit([&] {
      using NV = ore::NV;
      return OptimizationRemarkAnalysis(PassName, "FPRegisterBank", P.Insert)
             << NV("Function", F.getName()) << " "
             << NV("Action", Transform ? "carried" : "can carry") << " one "
             << NV("Source",
                   isa<Argument>(P.Scalar) ? "stack argument" : "FP temporary")
             << " in the GPR bank across " << NV("Calls", P.AcrossCall.count())
             << " calls, restoring it at " << NV("Uses", P.Sites.size())
             << " uses";
    });
  }
  if (!Transform || Chosen.empty())
    return PreservedAnalyses::all();
  for (unsigned Index : Chosen)
    moveToIntegerBank(Plans[Index]);
  PreservedAnalyses PA;
  PA.preserveSet<CFGAnalyses>();
  return PA;
}
} // namespace rssn
