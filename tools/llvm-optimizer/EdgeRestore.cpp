#include "llvm/ADT/DenseMap.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/CodeGen/LivePhysRegs.h"
#include "llvm/CodeGen/LiveRegUnits.h"
#include "llvm/CodeGen/MachineBlockFrequencyInfo.h"
#include "llvm/CodeGen/MachineBranchProbabilityInfo.h"
#include "llvm/CodeGen/MachineFrameInfo.h"
#include "llvm/CodeGen/MachineFunctionPass.h"
#include "llvm/CodeGen/MachineInstr.h"
#include "llvm/CodeGen/MachineMemOperand.h"
#include "llvm/CodeGen/MachineOptimizationRemarkEmitter.h"
#include "llvm/CodeGen/MachineRegisterInfo.h"
#include "llvm/CodeGen/TargetInstrInfo.h"
#include "llvm/CodeGen/TargetSubtargetInfo.h"
#include "llvm/IR/Function.h"
#include "llvm/Pass.h"
#include "llvm/Target/TargetMachine.h"
#include <cstdlib>
#include <optional>

using namespace llvm;

namespace {

constexpr char PassName[] = "rssn-edge-restore";
constexpr unsigned MaxRecipes = 16;
constexpr unsigned MaxSuffix = 64;

enum class Source { ConstantPool, IncomingArgument, SpillSlot, Zero };

struct Restore {
  MachineInstr *Instruction;
  Register Reg;
  Source Kind;
  unsigned Cost;
};

struct Edge {
  MachineBasicBlock *From;
  MachineBasicBlock *To;
  SmallVector<Restore, 8> Recipes;
};

std::optional<Restore> getRestore(MachineInstr &MI, const MachineFrameInfo &MFI,
                                  const TargetInstrInfo &TII,
                                  const TargetRegisterInfo &TRI,
                                  const MachineRegisterInfo &MRI) {
  if (MI.isBundled() || MI.getFlags() || MI.peekDebugInstrNum() ||
      MI.getNumExplicitDefs() != 1 || MI.mayStore() || MI.isTerminator() ||
      MI.hasUnmodeledSideEffects())
    return std::nullopt;
  const MachineOperand &Def = MI.getOperand(0);
  if (!Def.isReg() || !Def.getReg().isPhysical() || Def.getSubReg() ||
      Def.isDead() || Def.isInternalRead() || MRI.isReserved(Def.getReg()))
    return std::nullopt;
  const Register Reg = Def.getReg();
  if (any_of(MI.uses(), [&](const MachineOperand &Op) {
        return Op.isRegMask() ||
               (Op.isReg() && Op.getReg() &&
                (!Op.isUse() || !Op.getReg().isPhysical() ||
                 StringRef(TRI.getName(Op.getReg())) != "RIP"));
      }))
    return std::nullopt;

  if (!MI.mayLoad()) {
    // X86's independent vector-zero pseudos have no other operands or flags.
    const StringRef Opcode = TII.getName(MI.getOpcode());
    if (MI.getNumOperands() == 1 &&
        (Opcode == "AVX512_512_SET0" || Opcode == "AVX_SET0" ||
         Opcode == "V_SET0"))
      return Restore{&MI, Reg, Source::Zero, 1};
    return std::nullopt;
  }
  if (!MI.hasOneMemOperand() || MI.hasOrderedMemoryRef())
    return std::nullopt;
  const auto &Memory = **MI.memoperands_begin();
  const auto *PSV = Memory.getPseudoValue();
  if (!Memory.isLoad() || !PSV || !Memory.getSize().hasValue() ||
      Memory.getSize().isScalable())
    return std::nullopt;
  const StringRef Opcode = TII.getName(MI.getOpcode());
  if (PSV->isConstantPool() &&
      (Opcode.starts_with("VBROADCAST") || Opcode.starts_with("VPBROADCAST") ||
       TII.isTriviallyReMaterializable(MI)) &&
      any_of(MI.uses(), [](const MachineOperand &Op) { return Op.isCPI(); }))
    return Restore{&MI, Reg, Source::ConstantPool, 2};

  int FI = 0;
  const auto *Stack = dyn_cast<FixedStackPseudoSourceValue>(PSV);
  if (!Stack || TII.isLoadFromStackSlot(MI, FI) != Reg ||
      FI != Stack->getFrameIndex() || Memory.getOffset() != 0 ||
      MFI.isDeadObjectIndex(FI))
    return std::nullopt;
  const uint64_t Bytes = Memory.getSize().getValue().getFixedValue();
  if (Bytes > static_cast<uint64_t>(MFI.getObjectSize(FI)))
    return std::nullopt;
  if (MFI.isFixedObjectIndex(FI) && MFI.isImmutableObjectIndex(FI) &&
      !MFI.isAliasedObjectIndex(FI))
    return Restore{&MI, Reg, Source::IncomingArgument, 2};
  if (!MFI.isSpillSlotObjectIndex(FI) ||
      MFI.isStatepointSpillSlotObjectIndex(FI))
    return std::nullopt;
  // A folded wide transfer adds predictable traffic in addition to the load.
  return Restore{&MI, Reg, Source::SpillSlot, Bytes > 16 ? 3u : 2u};
}

bool crossesStore(const Restore &R, ArrayRef<const MachineInstr *> Stores,
                  const MachineFrameInfo &MFI) {
  if (R.Kind == Source::Zero || R.Kind == Source::ConstantPool)
    return false;
  const auto *Load = cast<FixedStackPseudoSourceValue>(
      (*R.Instruction->memoperands_begin())->getPseudoValue());
  for (const MachineInstr *Store : Stores) {
    if (none_of(Store->memoperands(),
                [](const MachineMemOperand *Mem) { return Mem->isStore(); }))
      return true;
    for (const MachineMemOperand *Mem : Store->memoperands()) {
      if (!Mem->isStore())
        continue;
      if (const auto *Other = dyn_cast_or_null<FixedStackPseudoSourceValue>(
              Mem->getPseudoValue())) {
        const int FI = Other->getFrameIndex();
        // Distinct local frame objects occupy disjoint storage here.
        if (FI == Load->getFrameIndex() ||
            (MFI.isFixedObjectIndex(FI) &&
             MFI.isFixedObjectIndex(Load->getFrameIndex())))
          return true;
      } else if (!Mem->getValue() ||
                 (R.Kind == Source::SpillSlot && Load->mayAlias(&MFI))) {
        return true;
      }
    }
  }
  return false;
}

SmallVector<Restore, 8> collectRecipes(MachineBasicBlock &From,
                                       MachineBasicBlock &To,
                                       MachineBasicBlock *Other) {
  auto &MF = *From.getParent();
  const auto &TII = *MF.getSubtarget().getInstrInfo();
  const auto &TRI = *MF.getSubtarget().getRegisterInfo();
  const auto &MFI = MF.getFrameInfo();
  LiveRegUnits Suffix(TRI), Exit(TRI);
  Exit.addLiveIns(To);
  if (Other)
    Suffix.addLiveIns(*Other);
  SmallVector<const MachineInstr *, 8> Stores;
  SmallVector<Restore, 8> Recipes;
  unsigned Visited = 0;
  for (MachineInstr &MI : reverse(From.instrs())) {
    if (++Visited > MaxSuffix || MI.isBundled() || MI.isInlineAsm() ||
        MI.isLabel() || MI.isCFIInstruction())
      break;
    if (MI.isCall())
      return Recipes;
    if (MI.hasUnmodeledSideEffects() &&
        MI.getOpcode() != TII.getCallFrameDestroyOpcode() &&
        MI.getOpcode() != TII.getCallFrameSetupOpcode())
      break;
    if (auto R = getRestore(MI, MFI, TII, TRI, MF.getRegInfo());
        R && Suffix.available(R->Reg) && !Exit.available(R->Reg) &&
        !crossesStore(*R, Stores, MFI) && Recipes.size() < MaxRecipes)
      Recipes.push_back(*R);
    if (MI.mayStore())
      Stores.push_back(&MI);
    Suffix.accumulate(MI);
    if (MI.isDebugInstr())
      for (const MachineOperand &MO : MI.operands())
        if (MO.isReg() && MO.getReg().isPhysical())
          Suffix.addReg(MO.getReg());
  }
  return {};
}

struct Group {
  MachineBasicBlock *Destination;
  SmallVector<Edge, 8> Inputs;
};

enum class Placement { Destination, Fallthrough, Detached };

bool sameRecipes(ArrayRef<Restore> A, ArrayRef<Restore> B) {
  return A.size() == B.size() &&
         llvm::equal(A, B, [](const Restore &L, const Restore &R) {
           return L.Kind == R.Kind && L.Cost == R.Cost &&
                  L.Instruction->isIdenticalTo(*R.Instruction);
         });
}

SmallVector<Group, 8> findGroups(MachineFunction &MF) {
  const auto &TII = *MF.getSubtarget().getInstrInfo();
  SmallVector<Group, 8> Groups;
  DenseMap<MachineBasicBlock *, SmallVector<unsigned, 2>> ByDestination;
  for (MachineBasicBlock &B : MF) {
    if (B.succ_empty() || B.succ_size() > 2 || B.isEHPad() ||
        B.isEHFuncletEntry() || B.hasAddressTaken())
      continue;
    MachineBasicBlock *Taken = nullptr, *Else = nullptr;
    SmallVector<MachineOperand, 4> Cond;
    if (TII.analyzeBranch(B, Taken, Else, Cond) || Taken == Else)
      continue;
    MachineBasicBlock *Fallthrough = B.getLogicalFallThrough();
    for (MachineBasicBlock *To : B.successors()) {
      // An implicit edge can be used without any branch rewrites if its
      // destination is private to this predecessor.
      const bool Private = To->pred_size() == 1;
      if (To == &B || (!Private && To == Fallthrough) ||
          (Private && B.succ_size() != 2) || To->isEHPad() ||
          To->isEHFuncletEntry() || To->hasAddressTaken() ||
          (To != Taken && To != Else && To != Fallthrough))
        continue;
      MachineBasicBlock *Other = nullptr;
      for (MachineBasicBlock *Succ : B.successors())
        if (Succ != To)
          Other = Succ;
      auto Recipes = collectRecipes(B, *To, Other);
      if (Recipes.empty())
        continue;
      llvm::sort(Recipes, [](const Restore &L, const Restore &R) {
        return L.Reg.id() < R.Reg.id();
      });
      auto &Indices = ByDestination[To];
      auto It = find_if(Indices, [&](unsigned I) {
        return sameRecipes(Groups[I].Inputs.front().Recipes, Recipes);
      });
      if (It == Indices.end()) {
        Indices.push_back(Groups.size());
        Groups.push_back({To, {}});
        It = std::prev(Indices.end());
      }
      Groups[*It].Inputs.push_back({&B, To, std::move(Recipes)});
    }
  }
  return Groups;
}

const Function *knownCallee(const MachineInstr &MI, const TargetInstrInfo &TII,
                            const TargetRegisterInfo &TRI) {
  if (MI.isBundled())
    return nullptr;
  const StringRef Opcode = TII.getName(MI.getOpcode());
  const MachineOperand *Callee = nullptr;
  if (Opcode == "CALL64pcrel32") {
    Callee = &MI.getOperand(0);
  } else if (Opcode == "CALL64m" && MI.getNumExplicitOperands() == 5) {
    const auto &Base = MI.getOperand(0);
    const auto &Scale = MI.getOperand(1);
    const auto &Index = MI.getOperand(2);
    const auto &Displacement = MI.getOperand(3);
    const auto &Segment = MI.getOperand(4);
    if (!Base.isReg() || !Base.getReg().isPhysical() ||
        StringRef(TRI.getName(Base.getReg())) != "RIP" || !Scale.isImm() ||
        Scale.getImm() != 1 || !Index.isReg() || Index.getReg() ||
        !Segment.isReg() || Segment.getReg() || !MI.hasOneMemOperand())
      return nullptr;
    const auto *Memory = (*MI.memoperands_begin())->getPseudoValue();
    if (!Memory || !Memory->isGOT() ||
        !any_of(TII.getSerializableDirectMachineOperandTargetFlags(),
                [&](const auto &Flag) {
                  return Flag.first == Displacement.getTargetFlags() &&
                         StringRef(Flag.second) == "x86-gotpcrel";
                }))
      return nullptr;
    Callee = &Displacement;
  }
  return Callee && Callee->isGlobal() && Callee->getOffset() == 0
             ? dyn_cast<Function>(Callee->getGlobal())
             : nullptr;
}

bool canPlaceBefore(MachineBasicBlock &Destination,
                    const TargetInstrInfo &TII) {
  const auto &MF = *Destination.getParent();
  if (&Destination == &MF.front() || Destination.isBeginSection() ||
      MF.getTarget().getBBSectionsType() != BasicBlockSection::None)
    return false;
  auto &Previous = *std::prev(Destination.getIterator());
  if (Previous.isEndSection() || !Previous.sameSection(&Destination))
    return false;
  const auto Last = Previous.getLastNonDebugInstr(/*SkipPseudoOp=*/false);
  if (Last != Previous.end() && Last->isBarrier() && !Last->isInlineAsm() &&
      !TII.isPredicated(*Last) && !Previous.getLogicalFallThrough())
    return true;
  if (!Previous.succ_empty())
    return false;
  // The backend can leave a direct noreturn call (and its frame-destroy
  // pseudo) with no barrier at the end of a cold block.
  for (const MachineInstr &MI : reverse(Previous.instrs())) {
    if (MI.isDebugInstr() || MI.getOpcode() == TII.getCallFrameDestroyOpcode())
      continue;
    const auto *F = knownCallee(MI, TII, *MF.getSubtarget().getRegisterInfo());
    return F && F->doesNotReturn() && F->doesNotThrow();
  }
  return false;
}

Placement choosePlacement(const Group &G, const TargetInstrInfo &TII) {
  if (&G.Destination->getParent()->front() != G.Destination &&
      G.Inputs.size() == G.Destination->pred_size() &&
      none_of(*G.Destination, [](const MachineInstr &MI) {
        return MI.isCFIInstruction() || MI.isEHLabel() ||
               MI.getFlag(MachineInstr::FrameSetup) ||
               MI.getFlag(MachineInstr::FrameDestroy);
      }))
    return Placement::Destination;
  return canPlaceBefore(*G.Destination, TII) ? Placement::Fallthrough
                                             : Placement::Detached;
}

double selectEdges(Group &G, Placement Where,
                   const MachineBlockFrequencyInfo &BFI,
                   const MachineBranchProbabilityInfo &BPI) {
  unsigned RecipeCost = 0;
  for (const Restore &R : G.Inputs.front().Recipes)
    RecipeCost += R.Cost;
  double Score = 0;
  const unsigned JumpCost = Where == Placement::Detached ? 1 : 0;
  llvm::erase_if(G.Inputs, [&](const Edge &E) {
    const auto Probability = BPI.getEdgeProbability(E.From, E.To);
    const double Taken = static_cast<double>(Probability.getNumerator()) /
                         Probability.getDenominator();
    const double Gain = RecipeCost * (1.0 - Taken) - JumpCost * Taken;
    if (Gain < 0)
      return true;
    Score += BFI.getBlockFreqRelativeToEntryBlock(E.From) * Gain;
    return false;
  });
  return Score;
}

void apply(MachineFunction &MF, Group &G, Placement Where,
           MachineBasicBlock &Anchor) {
  const auto &TII = *MF.getSubtarget().getInstrInfo();
  MachineBasicBlock *Shared = G.Destination;
  if (Where != Placement::Destination) {
    Shared = MF.CreateMachineBasicBlock();
    if (Where == Placement::Fallthrough)
      Shared->setSectionID(G.Destination->getSectionID());
    MF.insert(Where == Placement::Fallthrough ? G.Destination->getIterator()
                                              : std::next(Anchor.getIterator()),
              Shared);
    Shared->addSuccessor(G.Destination, BranchProbability::getOne());
  }
  const auto Insert = Shared->begin();
  for (const Restore &R : G.Inputs.front().Recipes) {
    MachineInstr *Copy = MF.CloneMachineInstr(R.Instruction);
    Copy->setDebugLoc(DebugLoc());
    Shared->insert(Insert, Copy);
  }
  if (Where == Placement::Detached)
    TII.insertBranch(*Shared, G.Destination, nullptr, {}, DebugLoc());
  recomputeLiveIns(*Shared);
  for (Edge &E : G.Inputs) {
    if (Where != Placement::Destination) {
      // Keep probabilities and skip debug instructions between terminators.
      for (MachineInstr &MI : E.From->terminators())
        for (MachineOperand &Op : MI.operands())
          if (Op.isMBB() && Op.getMBB() == G.Destination)
            Op.setMBB(Shared);
      E.From->replaceSuccessor(G.Destination, Shared);
    }
    for (Restore &R : E.Recipes)
      R.Instruction->eraseFromParent();
    recomputeLivenessFlags(*E.From);
  }
  recomputeLivenessFlags(*Shared);
}

template <bool Transform> class EdgeRestorePass : public MachineFunctionPass {
public:
  static char ID;
  EdgeRestorePass() : MachineFunctionPass(ID) {}

  StringRef getPassName() const override {
    return Transform ? "RSSN machine edge restores"
                     : "RSSN machine edge restore report";
  }

  void getAnalysisUsage(AnalysisUsage &AU) const override {
    AU.addRequired<MachineBlockFrequencyInfoWrapperPass>();
    AU.addRequired<MachineBranchProbabilityInfoWrapperPass>();
    AU.addRequired<MachineOptimizationRemarkEmitterPass>();
    if (!Transform)
      AU.setPreservesAll();
    MachineFunctionPass::getAnalysisUsage(AU);
  }

  bool runOnMachineFunction(MachineFunction &MF) override {
    const Function &F = MF.getFunction();
    const auto &Properties = MF.getProperties();
    const auto &MFI = MF.getFrameInfo();
    const char *Filter = std::getenv("RSSN_SPILL_FILTER");
    if (MF.getTarget().getTargetTriple().getArch() != Triple::x86_64 ||
        (Filter && !F.getName().contains(Filter)) || F.hasOptNone() ||
        F.hasOptSize() || F.hasMinSize() || MF.hasEHFunclets() ||
        !Properties.hasNoVRegs() || !Properties.hasNoPHIs() ||
        !Properties.hasTracksLiveness() || MFI.getStackSize() ||
        MFI.hasVarSizedObjects() || MFI.hasOpaqueSPAdjustment() || MF.empty() ||
        MF.size() > 4096)
      return false;
    auto Layout = reverse(MF);
    auto Anchor = find_if(Layout, [](MachineBasicBlock &B) {
      const auto Last = B.getLastNonDebugInstr();
      return Last != B.end() && Last->isBarrier() && !Last->isInlineAsm() &&
             !B.isEHPad() && !B.canFallThrough();
    });
    auto Groups = findGroups(MF);
    auto &BFI = getAnalysis<MachineBlockFrequencyInfoWrapperPass>().getMBFI();
    if (!BFI.getEntryFreq().getFrequency())
      return false;
    auto &BPI =
        getAnalysis<MachineBranchProbabilityInfoWrapperPass>().getMBPI();
    auto &ORE = getAnalysis<MachineOptimizationRemarkEmitterPass>().getORE();
    const auto &TII = *MF.getSubtarget().getInstrInfo();
    struct Plan {
      double Gain;
      Group *Candidate;
      Placement Where;
    };
    SmallVector<Plan, 8> Choices;
    for (Group &G : Groups) {
      const Placement Where = choosePlacement(G, TII);
      const double Gain = selectEdges(G, Where, BFI, BPI);
      if (G.Inputs.empty())
        continue;
      const size_t Removed =
          (G.Inputs.size() - 1) * G.Inputs.front().Recipes.size();
      if ((G.Inputs.size() == 1 && Where == Placement::Destination &&
           Gain > 0) ||
          (Removed > (Where == Placement::Detached ? 1u : 0u) &&
           (Where != Placement::Detached || Anchor != Layout.end())))
        Choices.push_back({Gain, &G, Where});
    }
    llvm::stable_sort(
        Choices, [](const auto &L, const auto &R) { return L.Gain > R.Gain; });
    SmallPtrSet<MachineBasicBlock *, 16> Used;
    SmallPtrSet<MachineBasicBlock *, 8> RewrittenDestinations;
    bool Changed = false;
    for (const Plan &Choice : Choices) {
      Group *G = Choice.Candidate;
      if (RewrittenDestinations.contains(G->Destination) ||
          any_of(G->Inputs,
                 [&](const Edge &E) { return Used.contains(E.From); }) ||
          (Choice.Where == Placement::Fallthrough &&
           !canPlaceBefore(*G->Destination, TII)))
        continue;
      RewrittenDestinations.insert(G->Destination);
      for (const Edge &E : G->Inputs)
        Used.insert(E.From);
      ORE.emit([&] {
        using NV = DiagnosticInfoOptimizationBase::Argument;
        return MachineOptimizationRemarkAnalysis(
                   PassName, "SharedRestores",
                   G->Inputs.front().Recipes.front().Instruction)
               << NV("Function", F.getName()) << ", "
               << (G->Inputs.size() == 1
                       ? (Transform ? "sank " : "can sink ")
                       : (Transform ? "shared " : "can share "))
               << NV("Recipes", G->Inputs.front().Recipes.size())
               << " restores on " << NV("Edges", G->Inputs.size())
               << (G->Inputs.size() == 1 ? " incoming edge, "
                                         : " incoming edges, ")
               << (Choice.Where == Placement::Destination   ? "in destination"
                   : Choice.Where == Placement::Fallthrough ? "on fallthrough"
                                                            : "with extra jump")
               << ", estimated gain "
               << NV("Gain", static_cast<float>(Choice.Gain));
      });
      if (Transform) {
        MachineBasicBlock &InsertionAnchor =
            Choice.Where == Placement::Detached ? *Anchor : *G->Destination;
        apply(MF, *G, Choice.Where, InsertionAnchor);
        Changed = true;
      }
    }
    return Changed;
  }
};

template <bool Transform> char EdgeRestorePass<Transform>::ID = 0;

RegisterPass<EdgeRestorePass<true>>
    Transform(PassName, "RSSN machine edge restores", false, false);
RegisterPass<EdgeRestorePass<false>> Report("rssn-edge-restore-report",
                                            "RSSN machine edge restore report",
                                            false, true);

} // namespace
