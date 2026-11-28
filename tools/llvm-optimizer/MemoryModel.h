#ifndef RSSN_LLVM_SPILL_MEMORY_MODEL_H
#define RSSN_LLVM_SPILL_MEMORY_MODEL_H

#include "llvm/ADT/ArrayRef.h"
#include "llvm/ADT/DenseMap.h"
#include "llvm/ADT/DenseSet.h"

namespace llvm {
class AAResults;
class AllocaInst;
class CallBase;
class Instruction;
class LoadInst;
} // namespace llvm

namespace rssn {

class MemoryModel {
  struct FrameObject {
    bool Private = false;
    llvm::DenseSet<const llvm::CallBase *> ReceivesAddress;
  };

  llvm::AAResults &AA;
  llvm::DenseMap<const llvm::AllocaInst *, FrameObject> Frames;
  llvm::DenseSet<std::pair<const llvm::AllocaInst *, const llvm::CallBase *>>
      ProvenBoundaries;

  static FrameObject inspect(const llvm::AllocaInst &Root);
  bool inaccessibleAt(const llvm::LoadInst &L, const llvm::CallBase &Call);

public:
  explicit MemoryModel(llvm::AAResults &AA) : AA(AA) {}

  bool unchanged(const llvm::LoadInst &Load,
                 llvm::ArrayRef<llvm::Instruction *> Code, unsigned From,
                 unsigned To);

  unsigned privateBoundaries() const { return ProvenBoundaries.size(); }
};

} // namespace rssn
#endif
