#ifndef RSSN_LLVM_SPILL_PRESSURE_H
#define RSSN_LLVM_SPILL_PRESSURE_H

#include "llvm/ADT/ArrayRef.h"
#include "llvm/ADT/BitVector.h"
#include "llvm/ADT/DenseMap.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/IR/BasicBlock.h"
#include <optional>

namespace llvm {
class AttributeList;
class DataLayout;
class Function;
class FunctionType;
class TargetTransformInfo;
} // namespace llvm

namespace rssn {

struct Pressure {
  int GP = 0;
  int SIMD = 0;
  int SIMDSpill = 0;
  int PreservableSIMD = 0;
  int VolatileSIMD = 0;
  int VolatileSIMDSpill = 0;

  Pressure &operator+=(Pressure Other);
  Pressure &operator-=(Pressure Other);
  int excess(Pressure Limit) const;
  int congestion(Pressure Limit) const;
  int callCost(Pressure Limit) const;
  static constexpr int Win64PreservedGP = 8;
  static constexpr int Win64PreservedXMM = 10;
  int win64CallExcess() const;
  int win64CallCost() const;
  int aapcs64CallCost() const;
  bool worsens(Pressure Before, Pressure Limit) const;
};

struct Summary {
  Pressure Peak;
  unsigned ExcessPoints = 0;
  unsigned CallSites = 0;
  unsigned GPThroughCalls = 0;
  unsigned GPCallExcess = 0;
  unsigned SIMDThroughCalls = 0;
  unsigned SIMDSpillThroughCalls = 0;
  uint64_t Congestion = 0;
};

class Liveness {
  struct State {
    llvm::BitVector Gen, Kill, PhiUses, In, Out;

    explicit State(unsigned Size)
        : Gen(Size), Kill(Size), PhiUses(Size), In(Size), Out(Size) {}
  };

  llvm::DenseMap<const llvm::Value *, unsigned> Values;
  llvm::SmallVector<Pressure, 128> Weights;
  llvm::DenseMap<const llvm::BasicBlock *, unsigned> Blocks;
  llvm::SmallVector<State, 16> States;
  const llvm::TargetTransformInfo &TTI;
  const llvm::DataLayout &DL;
  const bool AArch64;
  llvm::DenseMap<llvm::Type *, Pressure> TypeWeights;

  void track(const llvm::Value &V);

public:
  static constexpr unsigned MaxBlockSize = 4096;
  static bool supports(const llvm::Function &F);
  // A trial can change instruction placement while preserving the CFG and PHI
  // edges. Each original instruction must occur once in the resulting schedule.
  struct Schedule {
    const llvm::BasicBlock *Block;
    llvm::ArrayRef<llvm::Instruction *> Code;
  };

  Pressure Limit;
  static constexpr Pressure CallLimit{5, 0};

  Liveness(const llvm::TargetTransformInfo &TTI, const llvm::Function &F);
  bool compute(llvm::Function &F, llvm::ArrayRef<Schedule> Trial = {});
  void rebind(const llvm::Value &From, const llvm::Value &To);
  Pressure typeWeight(llvm::Type *Ty) const;
  unsigned id(const llvm::Value *V) const;
  Pressure weight(const llvm::Value *V) const;
  Pressure weight(unsigned ID) const { return Weights[ID]; }
  const llvm::BitVector &liveOut(const llvm::BasicBlock &B) const;
};

// Ordinary Win64 scalar slots. Complex signatures require target lowering.
std::optional<unsigned> win64ScalarStackSlots(const llvm::FunctionType &Ty,
                                              const llvm::AttributeList &Attrs);

bool isCallPressurePoint(const llvm::Instruction &I);
bool isWin64CallPressurePoint(const llvm::Instruction &I);

enum class CallPressureABI { SysV, Win64, AAPCS64 };

struct BlockPressure {
  llvm::SmallVector<llvm::Instruction *, 64> Code;
  llvm::DenseMap<const llvm::Instruction *, unsigned> Positions;
  llvm::DenseMap<const llvm::Value *, unsigned> FirstUse, LastUse;
  llvm::SmallVector<Pressure, 64> Before;
  llvm::SmallVector<Pressure, 64> Through;
  llvm::BitVector Calls;
  bool Pressured = false;

  BlockPressure(llvm::BasicBlock &B, const Liveness &Live,
                bool IncludeScalarLibm = false,
                llvm::ArrayRef<llvm::Instruction *> Trial = {},
                CallPressureABI ABI = CallPressureABI::SysV);
  void summarize(Summary &S, Pressure Limit) const;
  bool isLiveBefore(const llvm::Value *V, unsigned Point) const;
};

Summary summarize(llvm::Function &F, const Liveness &Live);

} // namespace rssn
#endif
