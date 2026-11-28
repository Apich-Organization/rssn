#include "llvm/CodeGen/MachineFrameInfo.h"
#include "llvm/CodeGen/MachineFunctionPass.h"
#include "llvm/CodeGen/MachineMemOperand.h"
#include "llvm/CodeGen/MachineOptimizationRemarkEmitter.h"
#include "llvm/IR/Function.h"
#include "llvm/Pass.h"
#include <cstdlib>

using namespace llvm;

namespace {

constexpr char PassName[] = "rssn-machine-stack-traffic";

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
    Reads Spill, Fixed, Constant;
    uint64_t Calls = 0, SpillWrites = 0;
    for (const MachineBasicBlock &Block : MF) {
      bool SawCall = false;
      for (const MachineInstr &MI : Block.instrs()) {
        if (MI.isBundle())
          continue;
        for (const MachineMemOperand *MMO : MI.memoperands()) {
          const auto *PSV = MMO->getPseudoValue();
          // Memory operands identify spill-slot accesses even when the target
          // folds the reload into a compare or another consuming operation.
          if (const auto *Slot = dyn_cast_or_null<FixedStackPseudoSourceValue>(PSV)) {
            const int FI = Slot->getFrameIndex();
            if (Frame.isSpillSlotObjectIndex(FI) &&
                !Frame.isStatepointSpillSlotObjectIndex(FI)) {
              SpillWrites += MMO->isStore();
              if (MMO->isLoad())
                Spill.add(*MMO, SawCall);
            } else if (MMO->isLoad() && Frame.isFixedObjectIndex(FI) &&
                       Frame.isImmutableObjectIndex(FI)) {
              Fixed.add(*MMO, SawCall);
            }
          } else if (MMO->isLoad() && PSV && PSV->isConstantPool()) {
            Constant.add(*MMO, SawCall);
          }
        }
        if (MI.isCall()) {
          ++Calls;
          SawCall = true;
        }
      }
    }

    auto &ORE = getAnalysis<MachineOptimizationRemarkEmitterPass>().getORE();
    ORE.emit([&] {
      using NV = DiagnosticInfoOptimizationBase::Argument;
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
    return false;
  }
};

char StackTrafficPass::ID = 0;
RegisterPass<StackTrafficPass> StackTraffic(
    PassName, "RSSN machine stack traffic report", false, true);

} // namespace
