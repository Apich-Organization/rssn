#include "LoopPreheaderPressure.h"
#include "Pressure.h"
#include "llvm/ADT/MapVector.h"
#include "llvm/ADT/STLExtras.h"
#include "llvm/Analysis/AssumptionCache.h"
#include "llvm/Analysis/BlockFrequencyInfo.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/OptimizationRemarkEmitter.h"
#include "llvm/Analysis/TargetTransformInfo.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/Instructions.h"
#include "llvm/IR/Module.h"
#include "llvm/Transforms/Utils/Local.h"

using namespace llvm;

namespace rssn {
namespace {

constexpr char PassName[] = "rssn-loop-preheader-pressure";
constexpr unsigned MaxMoves = 512;
constexpr unsigned MaxUses = 32;
using Placements = MapVector<Instruction *, BasicBlock *>;
using BlockSchedules = DenseMap<BasicBlock *, SmallVector<Instruction *, 16>>;

struct PlacementPlan {
  Placements Moves;
  Placements BackendCompares;
};

std::optional<CallPressureABI> preheaderABI(const Module &M) {
  const Triple &Target = M.getTargetTriple();
  if (Target.getArch() == Triple::x86_64) {
    if (Target.isOSWindows())
      return CallPressureABI::Win64;
    if (Target.isOSLinux() && Target.getEnvironment() != Triple::GNUX32)
      return CallPressureABI::SysV;
  }
  if (Target.getArch() == Triple::aarch64 && !Target.isWindowsArm64EC() &&
      Target.getEnvironment() != Triple::GNUILP32 &&
      M.getDataLayout().getPointerSizeInBits() == 64 &&
      (Target.isOSLinux() || Target.isOSDarwin() || Target.isOSWindows() ||
       Target.isOSFreeBSD()))
    return CallPressureABI::AAPCS64;
  return std::nullopt;
}

int callCost(const Pressure &Need, CallPressureABI ABI) {
  switch (ABI) {
  case CallPressureABI::SysV:
    return Need.callCost(Liveness::CallLimit);
  case CallPressureABI::Win64:
    return Need.win64CallCost();
  case CallPressureABI::AAPCS64:
    return Need.aapcs64CallCost();
  }
  llvm_unreachable("unknown call pressure ABI");
}

bool isPlacementCandidate(Instruction &I, const DominatorTree &DT,
                           AssumptionCache &AC) {
  if (I.use_empty() || I.hasNUsesOrMore(MaxUses))
    return false;
  constexpr unsigned IntegerOps[] = {
      Instruction::Add, Instruction::Sub, Instruction::Mul, Instruction::And,
      Instruction::Or,  Instruction::Xor, Instruction::Shl, Instruction::LShr,
      Instruction::AShr};
  if (!isa<ICmpInst, GetElementPtrInst, BitCastInst, ZExtInst, SExtInst,
           TruncInst, ShuffleVectorInst, InsertElementInst>(I) &&
      !(I.getType()->isIntOrIntVectorTy() &&
        is_contained(IntegerOps, I.getOpcode())))
    return false;
  // Repeating the operation can propagate the same poison. Its operands
  // must not supply a new undef choice after crossing the loop boundary.
  return all_of(I.operand_values(), [&](const Value *V) {
    return isGuaranteedNotToBeUndef(V, &AC, &I, &DT);
  });
}

BasicBlock *destination(Instruction &I, const Loop &L, const LoopInfo &LI,
                        const DominatorTree &DT, const Placements &Planned) {
  BasicBlock *Target = nullptr;
  for (Use &U : I.uses()) {
    auto *User = dyn_cast<Instruction>(U.getUser());
    if (!User)
      return nullptr;
    BasicBlock *B = User->getParent();
    if (auto *Phi = dyn_cast<PHINode>(User))
      B = Phi->getIncomingBlock(U.getOperandNo());
    else if (auto It = Planned.find(User); It != Planned.end())
      B = It->second;
    if (!L.contains(B) || !DT.isReachableFromEntry(B))
      return nullptr;
    Target = Target ? DT.findNearestCommonDominator(Target, B) : B;
  }
  // Keep the recipe out of nested cycles and off the unconditional header.
  while (Target && L.contains(Target) && LI.getLoopFor(Target) != &L)
    Target = DT.getNode(Target)->getIDom()->getBlock();
  if (!Target || !L.contains(Target) || Target == L.getHeader() ||
      Target->isEHPad() || Target->getFirstInsertionPt() == Target->end())
    return nullptr;
  return Target;
}

bool isBackendLocalCompare(const Instruction &I, const BasicBlock &UseBlock) {
  const auto *Cmp = dyn_cast<ICmpInst>(&I);
  if (!Cmp || !Cmp->getType()->isIntegerTy(1) || !Cmp->hasOneUse())
    return false;
  const auto *OperandTy = dyn_cast<IntegerType>(Cmp->getOperand(0)->getType());
  const auto *Branch = dyn_cast<CondBrInst>(*Cmp->user_begin());
  return OperandTy && OperandTy->getBitWidth() <= 64 && Branch &&
         Branch->getParent() == &UseBlock;
}

PlacementPlan planPlacements(const LoopInfo &LI, const DominatorTree &DT,
                             const BlockFrequencyInfo &BFI,
                             AssumptionCache &AC) {
  PlacementPlan Plan;
  for (Loop *L : LI.getLoopsInPreorder()) {
    BasicBlock *Preheader = L->getLoopPreheader();
    if (!Preheader || Preheader->isEHPad() ||
        !DT.isReachableFromEntry(Preheader))
      continue;
    const auto SourceFrequency = BFI.getBlockFreq(Preheader);
    // Consider use-side instructions first so shared address and broadcast
    // chains can follow their consumers into the same arm.
    for (Instruction &I : reverse(*Preheader)) {
      if (!isPlacementCandidate(I, DT, AC))
        continue;
      BasicBlock *Target = destination(I, *L, LI, DT, Plan.Moves);
      if (!Target)
        continue;
      const auto TargetFrequency = BFI.getBlockFreq(Target);
      if (!TargetFrequency.getFrequency() || TargetFrequency > SourceFrequency)
        continue;
      Plan.Moves.insert({&I, Target});
      // CodeGenPrepare already places these compares in their branch's
      // block. Count that placement in the baseline as well as the trial.
      if (isBackendLocalCompare(I, *Target))
        Plan.BackendCompares.insert({&I, Target});
      if (Plan.Moves.size() == MaxMoves)
        return Plan;
    }
  }
  return Plan;
}

Placements requiredMoves(const PlacementPlan &Plan, const DominatorTree &DT) {
  Placements Required;
  for (const auto &[I, Dest] : Plan.Moves) {
    // A modeled comparison can stay in place if the actual move of an
    // operand will still dominate it. Keep it in a joint move otherwise.
    if (Plan.BackendCompares.contains(I) && none_of(I->operand_values(), [&](Value *V) {
          const auto *NewDef = Plan.Moves.lookup(dyn_cast<Instruction>(V));
          return NewDef && !DT.dominates(NewDef, I->getParent());
        }))
      continue;
    Required.insert({I, Dest});
  }
  return Required;
}

BlockSchedules buildSchedules(const Placements &Planned) {
  BlockSchedules Blocks, Insertions;
  for (auto [I, Target] : reverse(Planned)) {
    Insertions[Target].push_back(I);
    Blocks.try_emplace(Target);
    Blocks.try_emplace(I->getParent());
  }
  for (auto &[B, Code] : Blocks) {
    const auto At = B->getFirstInsertionPt();
    const auto Added = Insertions.find(B);
    for (auto &I : *B) {
      if (Added != Insertions.end() && I.getIterator() == At)
        append_range(Code, Added->second);
      if (Planned.find(&I) == Planned.end())
        Code.push_back(&I);
    }
  }
  return Blocks;
}

double preservationPressure(Function &F, const Liveness &Live,
                             const BlockFrequencyInfo &BFI,
                             const DominatorTree &DT,
                             const BlockSchedules &Blocks,
                             CallPressureABI ABI) {
  const double Entry = std::max<uint64_t>(1, BFI.getEntryFreq().getFrequency());
  double Cost = 0;
  for (BasicBlock &B : F) {
    if (!DT.isReachableFromEntry(&B))
      continue;
    const auto It = Blocks.find(&B);
    const ArrayRef<Instruction *> Code =
        It == Blocks.end() ? ArrayRef<Instruction *>()
                           : ArrayRef<Instruction *>(It->second);
    const BlockPressure P(B, Live, true, Code, ABI);
    const double Frequency = BFI.getBlockFreq(&B).getFrequency() / Entry;
    for (unsigned Point : P.Calls.set_bits())
      Cost += Frequency * callCost(P.Through[Point], ABI);
  }
  return Cost;
}

} // namespace

PreservedAnalyses LoopPreheaderPressurePass::run(Function &F,
                                                 FunctionAnalysisManager &AM) {
  const auto ABI = preheaderABI(*F.getParent());
  if (!ABI || !F.getName().contains(Filter) || F.isDeclaration() || F.hasOptNone() ||
      F.hasOptSize() || F.hasMinSize() || F.hasGC() || !Liveness::supports(F))
    return PreservedAnalyses::all();
  const auto &LI = AM.getResult<LoopAnalysis>(F);
  if (LI.empty())
    return PreservedAnalyses::all();
  const auto &DT = AM.getResult<DominatorTreeAnalysis>(F);
  const auto &BFI = AM.getResult<BlockFrequencyAnalysis>(F);
  auto &AC = AM.getResult<AssumptionAnalysis>(F);
  const auto Plan = planPlacements(LI, DT, BFI, AC);
  const auto Moves = requiredMoves(Plan, DT);
  if (Moves.empty())
    return PreservedAnalyses::all();

  const auto &TTI = AM.getResult<TargetIRAnalysis>(F);
  auto score = [&](const Placements &Placements) -> std::optional<double> {
    const auto Blocks = buildSchedules(Placements);
    SmallVector<Liveness::Schedule, 16> Schedule;
    for (const auto &[B, Instructions] : Blocks) {
      if (Instructions.size() > Liveness::MaxBlockSize)
        return std::nullopt;
      Schedule.push_back({B, Instructions});
    }
    Liveness Live(TTI, F);
    if (!Live.compute(F, Schedule))
      return std::nullopt;
    return preservationPressure(F, Live, BFI, DT, Blocks, *ABI);
  };

  const auto Original = score(Plan.BackendCompares);
  if (!Original || *Original == 0)
    return PreservedAnalyses::all();
  const auto New = score(Plan.Moves);
  if (!New)
    return PreservedAnalyses::all();
  const bool Eligible = *New < *Original;
  const bool Apply = Transform && Eligible;
  auto &ORE = AM.getResult<OptimizationRemarkEmitterAnalysis>(F);
  ORE.emit([&] {
    using NV = ore::NV;
    return OptimizationRemarkAnalysis(PassName, "Placement", &F)
           << NV("Function", F.getName()) << ", preheader placements "
           << NV("Placements", Moves.size()) << ", estimated call pressure "
           << NV("Before", static_cast<float>(*Original)) << " -> "
           << NV("After", static_cast<float>(*New)) << ", "
           << (Apply ? "applied" : Eligible ? "eligible" : "rejected")
           << ", backend-local branch comparisons "
           << NV("BackendCompares", Plan.BackendCompares.size());
  });
  if (!Apply)
    return PreservedAnalyses::all();
  for (auto [I, B] : Moves) {
    salvageDebugInfo(*I);
    I->moveBefore(*B, B->getFirstInsertionPt());
  }
  PreservedAnalyses PA;
  PA.preserveSet<CFGAnalyses>();
  return PA;
}

} // namespace rssn
