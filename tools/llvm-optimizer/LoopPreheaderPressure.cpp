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

Placements planPlacements(const LoopInfo &LI, const DominatorTree &DT,
                           const BlockFrequencyInfo &BFI,
                           AssumptionCache &AC) {
  Placements Planned;
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
      BasicBlock *Target = destination(I, *L, LI, DT, Planned);
      if (!Target)
        continue;
      const auto TargetFrequency = BFI.getBlockFreq(Target);
      if (!TargetFrequency.getFrequency() || TargetFrequency > SourceFrequency)
        continue;
      Planned.insert({&I, Target});
      if (Planned.size() == MaxMoves)
        return Planned;
    }
  }
  return Planned;
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
                             const BlockSchedules &Blocks) {
  const double Entry = std::max<uint64_t>(1, BFI.getEntryFreq().getFrequency());
  double Cost = 0;
  for (BasicBlock &B : F) {
    if (!DT.isReachableFromEntry(&B))
      continue;
    const auto It = Blocks.find(&B);
    const ArrayRef<Instruction *> Code =
        It == Blocks.end() ? ArrayRef<Instruction *>()
                           : ArrayRef<Instruction *>(It->second);
    const BlockPressure P(B, Live, true, Code, CallPressureABI::Win64);
    const double Frequency = BFI.getBlockFreq(&B).getFrequency() / Entry;
    for (unsigned Point : P.Calls.set_bits()) {
      const auto &Need = P.Through[Point];
      Cost += Frequency * Need.win64CallCost();
    }
  }
  return Cost;
}

} // namespace

PreservedAnalyses LoopPreheaderPressurePass::run(Function &F,
                                                 FunctionAnalysisManager &AM) {
  const Triple &Target = F.getParent()->getTargetTriple();
  if (Target.getArch() != Triple::x86_64 || !Target.isOSWindows() ||
      !F.getName().contains(Filter) || F.isDeclaration() || F.hasOptNone() ||
      F.hasOptSize() || F.hasMinSize() || F.hasGC() || !Liveness::supports(F))
    return PreservedAnalyses::all();
  const auto &LI = AM.getResult<LoopAnalysis>(F);
  if (LI.empty())
    return PreservedAnalyses::all();
  const auto &DT = AM.getResult<DominatorTreeAnalysis>(F);
  const auto &BFI = AM.getResult<BlockFrequencyAnalysis>(F);
  auto &AC = AM.getResult<AssumptionAnalysis>(F);
  const auto Planned = planPlacements(LI, DT, BFI, AC);
  if (Planned.empty())
    return PreservedAnalyses::all();

  const auto &TTI = AM.getResult<TargetIRAnalysis>(F);
  Liveness Before(TTI, F);
  if (!Before.compute(F))
    return PreservedAnalyses::all();
  const double Original = preservationPressure(F, Before, BFI, DT, BlockSchedules());
  if (Original == 0)
    return PreservedAnalyses::all();

  const auto Blocks = buildSchedules(Planned);
  SmallVector<Liveness::Schedule, 16> Trial;
  for (const auto &[B, Instructions] : Blocks) {
    if (Instructions.size() > Liveness::MaxBlockSize)
      return PreservedAnalyses::all();
    Trial.push_back({B, Instructions});
  }
  Liveness After(TTI, F);
  if (!After.compute(F, Trial))
    return PreservedAnalyses::all();
  const double New = preservationPressure(F, After, BFI, DT, Blocks);
  const bool Eligible = New < Original;
  const bool Apply = Transform && Eligible;
  auto &ORE = AM.getResult<OptimizationRemarkEmitterAnalysis>(F);
  ORE.emit([&] {
    using NV = ore::NV;
    return OptimizationRemarkAnalysis(PassName, "Placement", &F)
           << NV("Function", F.getName()) << ", preheader placements "
           << NV("Placements", Planned.size()) << ", estimated call pressure "
           << NV("Before", static_cast<float>(Original)) << " -> "
           << NV("After", static_cast<float>(New)) << ", "
           << (Apply ? "applied" : Eligible ? "eligible" : "rejected");
  });
  if (!Apply)
    return PreservedAnalyses::all();
  for (auto [I, B] : Planned) {
    salvageDebugInfo(*I);
    I->moveBefore(*B, B->getFirstInsertionPt());
  }
  PreservedAnalyses PA;
  PA.preserveSet<CFGAnalyses>();
  return PA;
}

} // namespace rssn
