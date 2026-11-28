#ifndef RSSN_LLVM_LOOP_EXIT_COMPRESSION_H
#define RSSN_LLVM_LOOP_EXIT_COMPRESSION_H

#include "llvm/IR/PassManager.h"
#include <string>

namespace rssn {
class LoopExitCompressionPass
    : public llvm::PassInfoMixin<LoopExitCompressionPass> {
  bool Transform;
  std::string Filter;

public:
  LoopExitCompressionPass(bool Transform, llvm::StringRef Filter)
      : Transform(Transform), Filter(Filter.str()) {}

  llvm::PreservedAnalyses run(llvm::Function &F,
                              llvm::FunctionAnalysisManager &AM);
};
} // namespace rssn
#endif
