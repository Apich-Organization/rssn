#include "llvm/ADT/DenseMap.h"
#include "llvm/ADT/STLExtras.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/CodeGen/MachineDominators.h"
#include "llvm/CodeGen/MachineFrameInfo.h"
#include "llvm/CodeGen/MachineFunctionPass.h"
#include "llvm/CodeGen/MachineMemOperand.h"
#include "llvm/CodeGen/MachineOptimizationRemarkEmitter.h"
#include "llvm/CodeGen/TargetInstrInfo.h"
#include "llvm/CodeGen/TargetRegisterInfo.h"
#include "llvm/CodeGen/TargetSubtargetInfo.h"
#include "llvm/IR/Function.h"
#include "llvm/Pass.h"
#include "llvm/Target/TargetMachine.h"
#include <cstdlib>
#include <optional>

using namespace llvm;

namespace {

constexpr char PassName[] = "rssn-machine-stack-traffic";

enum class OriginKind { Unreachable, Entry, Call, Ambiguous };

struct CallOrigin {
  OriginKind Kind = OriginKind::Unreachable;
  const MachineInstr *Call = nullptr;

  bool merge(CallOrigin Other) {
    if (Kind == OriginKind::Unreachable) {
      *this = Other;
      return true;
    }
    if (Kind == OriginKind::Ambiguous ||
        (Kind == Other.Kind && Call == Other.Call))
      return false;
    *this = {OriginKind::Ambiguous};
    return true;
  }
};

CallOrigin originAfter(const MachineInstr &MI) {
  if (MI.isBundled())
    return {OriginKind::Ambiguous};
  return {OriginKind::Call, &MI};
}

DenseMap<const MachineBasicBlock *, CallOrigin>
findBlockEntryOrigins(const MachineFunction &MF) {
  DenseMap<const MachineBasicBlock *, CallOrigin> Entry, LastLocalCall;
  SmallVector<const MachineBasicBlock *, 32> Worklist;
  auto mergeInto = [&](const MachineBasicBlock *Block, CallOrigin Origin) {
    if (Entry[Block].merge(Origin))
      Worklist.push_back(Block);
  };

  mergeInto(&MF.front(), {OriginKind::Entry});
  for (const MachineBasicBlock &Block : MF) {
    if (Block.isEHPad() || Block.isEHFuncletEntry() || Block.isEHContTarget() ||
        Block.hasAddressTaken())
      mergeInto(&Block, {OriginKind::Ambiguous});
    for (const MachineInstr &MI : reverse(Block.instrs())) {
      if (MI.isCall(MachineInstr::IgnoreBundle)) {
        LastLocalCall[&Block] = originAfter(MI);
        break;
      }
    }
  }

  // A block's entry changes at most twice. A local call kills any incoming
  // fact. This keeps loops and reconverging paths linear in the CFG size.
  while (!Worklist.empty()) {
    const auto *Block = Worklist.pop_back_val();
    auto Local = LastLocalCall.find(Block);
    const CallOrigin Out = Local == LastLocalCall.end() ? Entry.lookup(Block)
                                                       : Local->second;
    for (const MachineBasicBlock *Successor : Block->successors())
      mergeInto(Successor, Out);
  }
  return Entry;
}

enum class MemoryKind { Other, Spill, Fixed, Constant };

struct Reads {
  uint64_t Count = 0;
  uint64_t AfterCall = 0;
  uint64_t KnownBytes = 0;
  uint64_t Wide = 0;

  void add(const MachineMemOperand &MMO, bool SawCall) {
    ++Count;
    AfterCall += SawCall;
    const auto Size = MMO.getSize();
    if (Size.hasValue() && !Size.isScalable()) {
      const uint64_t Bytes = Size.getValue().getFixedValue();
      KnownBytes += Bytes;
      Wide += Bytes > 8;
    }
  }
};

struct Traffic {
  Reads Spill, Fixed, Constant;

  void add(MemoryKind Kind, const MachineMemOperand &MMO, bool FollowsCall) {
    switch (Kind) {
    case MemoryKind::Spill:
      Spill.add(MMO, FollowsCall);
      break;
    case MemoryKind::Fixed:
      Fixed.add(MMO, FollowsCall);
      break;
    case MemoryKind::Constant:
      Constant.add(MMO, FollowsCall);
      break;
    case MemoryKind::Other:
      break;
    }
  }
};

bool isSpillSlot(const MachineFrameInfo &Frame, int FI) {
  return Frame.isSpillSlotObjectIndex(FI) &&
         !Frame.isStatepointSpillSlotObjectIndex(FI);
}

const FixedStackPseudoSourceValue *frameSource(const MachineMemOperand &MMO) {
  return dyn_cast_or_null<FixedStackPseudoSourceValue>(MMO.getPseudoValue());
}

MemoryKind classify(const MachineMemOperand &MMO, const MachineFrameInfo &Frame) {
  const auto *PSV = MMO.getPseudoValue();
  if (const auto *Slot = frameSource(MMO)) {
    const int FI = Slot->getFrameIndex();
    if (isSpillSlot(Frame, FI))
      return MemoryKind::Spill;
    if (Frame.isFixedObjectIndex(FI) && Frame.isImmutableObjectIndex(FI))
      return MemoryKind::Fixed;
  } else if (PSV && PSV->isConstantPool()) {
    return MemoryKind::Constant;
  }
  return MemoryKind::Other;
}

struct SlotAccesses {
  const MachineInstr *Store = nullptr;
  SmallVector<const MachineInstr *, 4> Loads;
  bool Unknown = false;
};

const MachineInstr *findConstantRecipe(const MachineInstr &Store, Register Reg,
                                       const TargetInstrInfo &TII,
                                       const TargetRegisterInfo &TRI) {
  unsigned Remaining = 64;
  for (const MachineInstr &MI : reverse(make_range(
           Store.getParent()->instr_begin(), Store.getIterator()))) {
    if (MI.isDebugInstr() || MI.isCFIInstruction())
      continue;
    if (!Remaining-- || MI.isCall() || MI.isBundled() ||
        MI.hasUnmodeledSideEffects())
      return nullptr;
    if (!MI.modifiesRegister(Reg, &TRI))
      continue;
    if (MI.getNumExplicitDefs() != 1 || !MI.getOperand(0).isReg() ||
        MI.getOperand(0).getReg() != Reg || MI.getOperand(0).getSubReg() ||
        MI.mayStore() || !MI.hasOneMemOperand() || MI.hasOrderedMemoryRef())
      return nullptr;
    const auto *PSV = (*MI.memoperands_begin())->getPseudoValue();
    const StringRef Opcode = TII.getName(MI.getOpcode());
    if (MI.mayLoad() && PSV && PSV->isConstantPool() &&
        (TII.isTriviallyReMaterializable(MI) || Opcode.starts_with("VPBROADCAST") ||
         Opcode.starts_with("VBROADCAST") || Opcode.starts_with("VMOVDDUP")) &&
        any_of(MI.uses(), [](const MachineOperand &Op) { return Op.isCPI(); }) &&
        none_of(MI.uses(), [&](const MachineOperand &Op) {
          return Op.isRegMask() ||
                 (Op.isReg() && Op.getReg() &&
                  (!Op.isUse() || StringRef(TRI.getName(Op.getReg())) != "RIP"));
        }))
      return &MI;
    return nullptr;
  }
  return nullptr;
}

void reportConstantSpills(MachineFunction &MF, MachineOptimizationRemarkEmitter &ORE) {
  // Limit the dominance proof to ordinary entries and control-flow edges.
  if (MF.getTarget().getTargetTriple().getArch() != Triple::x86_64 ||
      any_of(MF, [](const MachineBasicBlock &BB) {
        return BB.isEHPad() || BB.isEHFuncletEntry() || BB.isEHContTarget() ||
               BB.hasAddressTaken();
      }))
    return;
  const auto &Frame = MF.getFrameInfo();
  const auto &TII = *MF.getSubtarget().getInstrInfo();
  const auto &TRI = *MF.getSubtarget().getRegisterInfo();
  DenseMap<int, SlotAccesses> Slots;
  for (const MachineBasicBlock &BB : MF) {
    for (const MachineInstr &MI : BB.instrs()) {
      if (MI.isBundle())
        continue;
      SmallVector<int, 2> Indices;
      auto addIndex = [&](int FI) {
        if (isSpillSlot(Frame, FI) && !Frame.isDeadObjectIndex(FI) &&
            !Frame.isAliasedObjectIndex(FI) && !is_contained(Indices, FI))
          Indices.push_back(FI);
      };
      for (const MachineOperand &Op : MI.operands())
        if (Op.isFI())
          addIndex(Op.getIndex());
      for (const auto *MMO : MI.memoperands())
        if (const auto *Slot = frameSource(*MMO))
          addIndex(Slot->getFrameIndex());

      for (int FI : Indices) {
        auto &Access = Slots[FI];
        const auto *MMO = MI.hasOneMemOperand() ? *MI.memoperands_begin() : nullptr;
        const auto *Slot = MMO ? frameSource(*MMO) : nullptr;
        if (Indices.size() != 1 || MI.isBundled() || !Slot ||
            Slot->getFrameIndex() != FI || MI.hasOrderedMemoryRef() ||
            MMO->getOffset() != 0 ||
            MMO->getSize() != LocationSize::precise(Frame.getObjectSize(FI))) {
          Access.Unknown = true;
        } else if (MMO->isStore()) {
          Access.Unknown |= Access.Store != nullptr || MMO->isLoad();
          Access.Store = &MI;
        } else if (MMO->isLoad()) {
          Access.Loads.push_back(&MI);
        } else {
          Access.Unknown = true;
        }
      }
    }
  }

  std::optional<MachineDominatorTree> DomTree;
  for (int FI = 0, End = Frame.getObjectIndexEnd(); FI < End; ++FI) {
    const auto Found = Slots.find(FI);
    if (Found == Slots.end())
      continue;
    const SlotAccesses &Access = Found->second;
    if (Access.Unknown || !Access.Store || Access.Loads.empty())
      continue;
    int StoredIndex = 0;
    const Register Reg = TII.isStoreToStackSlot(*Access.Store, StoredIndex);
    if (!Reg || StoredIndex != FI)
      continue;
    const MachineInstr *Recipe = findConstantRecipe(*Access.Store, Reg, TII, TRI);
    if (!Recipe)
      continue;
    if (!DomTree)
      DomTree.emplace(MF);
    if (!DomTree->isReachableFromEntry(Access.Store->getParent()) ||
        any_of(Access.Loads, [&](const MachineInstr *Load) {
          return !DomTree->isReachableFromEntry(Load->getParent()) ||
                 !DomTree->dominates(Access.Store, Load);
        }))
      continue;
    const unsigned Direct = count_if(Access.Loads, [&](const MachineInstr *Load) {
      int Index = 0;
      return TII.isLoadFromStackSlot(*Load, Index) && Index == FI;
    });
    ORE.emit([&] {
      using NV = DiagnosticInfoOptimizationBase::Argument;
      return MachineOptimizationRemarkAnalysis(PassName, "ConstantSpill", Access.Store)
             << NV("Function", MF.getName()) << " constant spill slot " << NV("FrameIndex", FI)
             << ", " << NV("SlotBytes", Frame.getObjectSize(FI)) << " bytes from "
             << NV("Recipe", TII.getName(Recipe->getOpcode())) << ", "
             << NV("Reads", Access.Loads.size()) << " reads ("
             << NV("DirectReads", Direct) << " direct reloads)";
    });
  }
}

StringRef calleeName(const MachineInstr &MI) {
  for (const MachineOperand &Op : MI.explicit_operands()) {
    if (Op.isGlobal())
      return Op.getGlobal()->getName();
    if (Op.isSymbol())
      return Op.getSymbolName();
  }
  return "<unknown>";
}

class StackTrafficPass : public MachineFunctionPass {
public:
  static char ID;
  StackTrafficPass() : MachineFunctionPass(ID) {}

  void getAnalysisUsage(AnalysisUsage &AU) const override {
    AU.addRequired<MachineOptimizationRemarkEmitterPass>();
    AU.setPreservesAll();
    MachineFunctionPass::getAnalysisUsage(AU);
  }

  bool runOnMachineFunction(MachineFunction &MF) override {
    const Function &F = MF.getFunction();
    const char *Filter = std::getenv("RSSN_SPILL_FILTER");
    if (!MF.getProperties().hasNoVRegs() || MF.empty() ||
        (Filter && !F.getName().contains(Filter)))
      return false;

    const MachineFrameInfo &Frame = MF.getFrameInfo();
    auto &ORE = getAnalysis<MachineOptimizationRemarkEmitterPass>().getORE();
    const char *Detail = std::getenv("RSSN_MACHINE_STACK_DETAIL");
    const bool ReportCalls = Detail && StringRef(Detail) == "1" &&
                             ORE.allowExtraAnalysis(PassName);
    const auto Entry = ReportCalls
                           ? findBlockEntryOrigins(MF)
                           : DenseMap<const MachineBasicBlock *, CallOrigin>();
    DenseMap<const MachineInstr *, Traffic> ByCall;
    Traffic Totals, Ambiguous;
    uint64_t Calls = 0, SpillWrites = 0;
    for (const MachineBasicBlock &Block : MF) {
      bool SawCall = false;
      CallOrigin Origin = Entry.lookup(&Block);
      for (const MachineInstr &MI : Block.instrs()) {
        if (MI.isBundle())
          continue;
        for (const MachineMemOperand *MMO : MI.memoperands()) {
          // Memory operands identify spill-slot accesses even when the target
          // folds the reload into a compare or another consuming operation.
          const auto Kind = classify(*MMO, Frame);
          SpillWrites += Kind == MemoryKind::Spill && MMO->isStore();
          if (!MMO->isLoad() || Kind == MemoryKind::Other)
            continue;
          Totals.add(Kind, *MMO, SawCall);
          if (Origin.Kind == OriginKind::Unreachable)
            continue;
          if (Origin.Kind == OriginKind::Ambiguous || MI.isBundled())
            Ambiguous.add(Kind, *MMO, false);
          else if (Origin.Kind == OriginKind::Call)
            ByCall[Origin.Call].add(Kind, *MMO, !SawCall);
        }
        if (MI.isCall(MachineInstr::IgnoreBundle)) {
          ++Calls;
          SawCall = true;
          if (ReportCalls && Origin.Kind != OriginKind::Unreachable)
            Origin = originAfter(MI);
        }
      }
    }

    using NV = DiagnosticInfoOptimizationBase::Argument;
    ORE.emit([&] {
      const auto &[Spill, Fixed, Constant] = Totals;
      return MachineOptimizationRemarkAnalysis(PassName, "StackTraffic",
                                               F.getSubprogram(), &MF.front())
             << NV("Function", F.getName()) << ", calls " << NV("Calls", Calls)
             << ", spill reads " << NV("SpillReads", Spill.Count)
             << " (" << NV("SpillAfterCall", Spill.AfterCall) << " after calls, "
             << NV("SpillWide", Spill.Wide) << " wider than 8 bytes, "
             << NV("SpillKnownBytes", Spill.KnownBytes) << " known bytes)"
             << ", spill writes " << NV("SpillWrites", SpillWrites)
             << ", immutable fixed-stack reads " << NV("FixedReads", Fixed.Count)
             << " (" << NV("FixedAfterCall", Fixed.AfterCall) << " after calls)"
             << ", constant-pool reads " << NV("ConstantReads", Constant.Count)
             << " (" << NV("ConstantAfterCall", Constant.AfterCall) << " after calls)";
    });
    if (ReportCalls) {
      for (const MachineBasicBlock &Block : MF) {
        unsigned CallIndex = 0;
        for (const MachineInstr &MI : Block.instrs()) {
          if (!MI.isCall(MachineInstr::IgnoreBundle))
            continue;
          ++CallIndex;
          const auto Found = ByCall.find(&MI);
          if (Found == ByCall.end())
            continue;
          ORE.emit([&] {
            const auto &[Spill, Fixed, Constant] = Found->second;
            return MachineOptimizationRemarkAnalysis(PassName, "CallTraffic", &MI)
                   << NV("Function", MF.getName()) << " bb." << NV("Block", Block.getNumber())
                   << " call " << NV("CallIndex", CallIndex) << " " << NV("Callee", calleeName(MI))
                   << ", spill reads " << NV("SpillReads", Spill.Count)
                   << " (" << NV("SpillAcrossBlocks", Spill.AfterCall) << " across blocks, "
                   << NV("SpillWide", Spill.Wide) << " wider than 8 bytes, "
                   << NV("SpillKnownBytes", Spill.KnownBytes) << " known bytes)"
                   << ", immutable fixed-stack reads " << NV("FixedReads", Fixed.Count)
                   << " (" << NV("FixedAcrossBlocks", Fixed.AfterCall) << " across blocks)"
                   << ", constant-pool reads " << NV("ConstantReads", Constant.Count)
                   << " (" << NV("ConstantAcrossBlocks", Constant.AfterCall) << " across blocks)";
          });
        }
      }
      ORE.emit([&] {
        return MachineOptimizationRemarkAnalysis(PassName, "AmbiguousCallTraffic",
                                                 F.getSubprogram(), &MF.front())
               << NV("Function", MF.getName()) << " ambiguous call origins, spill reads "
               << NV("SpillReads", Ambiguous.Spill.Count) << ", immutable fixed-stack reads "
               << NV("FixedReads", Ambiguous.Fixed.Count) << ", constant-pool reads "
               << NV("ConstantReads", Ambiguous.Constant.Count);
      });
      reportConstantSpills(MF, ORE);
    }
    return false;
  }
};

char StackTrafficPass::ID = 0;
RegisterPass<StackTrafficPass> StackTraffic(
    PassName, "RSSN machine stack traffic report", false, true);

} // namespace
