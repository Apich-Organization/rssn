#ifndef RSSN_LLVM_SPILL_PLANNING_H
#define RSSN_LLVM_SPILL_PLANNING_H

#include "MemoryModel.h"
#include "Pressure.h"
#include "llvm/ADT/DenseSet.h"
#include "llvm/IR/Instructions.h"
#include <optional>

namespace rssn {

constexpr unsigned MaxDistance = 192;
constexpr unsigned CallWeight = 8;
constexpr unsigned MaxPlanWidth = 4;

struct Changes {
  unsigned Sunk = 0;
  unsigned LoadsSunk = 0;
  unsigned LoadsMerged = 0;
  unsigned PacksFused = 0;
  unsigned LaneGroupsFused = 0;
  unsigned ReadsFused = 0;
  unsigned MotionGroups = 0;
  unsigned GroupMembers = 0;
  unsigned CrossBlock = 0;
  unsigned Copies = 0;

  bool any() const {
    return Sunk || LoadsSunk || LoadsMerged || PacksFused || LaneGroupsFused ||
           CrossBlock;
  }
};

using AddressKey = std::pair<const llvm::Value *, int64_t>;
std::optional<AddressKey> accessAddress(const llvm::LoadInst &Load);
bool simpleLoad(const llvm::LoadInst &L, const Liveness &Live);
bool cheapPure(const llvm::Instruction &I, const Liveness &Live,
               const llvm::TargetTransformInfo &TTI);
bool trySink(BlockPressure &P, const Liveness &Live,
             const llvm::TargetTransformInfo &TTI, MemoryModel &Memory,
             llvm::DenseSet<const llvm::Instruction *> &Moved, Changes &Result,
             unsigned RewriteBudget, unsigned PlanWidth);
bool tryFuseLaneLoads(llvm::BasicBlock &B, Liveness &Live,
                      const llvm::TargetTransformInfo &TTI, MemoryModel &Memory,
                      Changes &Result);
bool tryFusePackedLoads(llvm::BasicBlock &B, Liveness &Live,
                        const llvm::TargetTransformInfo &TTI,
                        MemoryModel &Memory, Changes &Result);

} // namespace rssn
#endif
