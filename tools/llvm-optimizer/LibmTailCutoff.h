#ifndef RSSN_LLVM_LIBM_TAIL_CUTOFF_H
#define RSSN_LLVM_LIBM_TAIL_CUTOFF_H

#include "llvm/IR/PassManager.h"
#include <string>

namespace rssn {

class LibmTailCutoffPass : public llvm::PassInfoMixin<LibmTailCutoffPass> {
  bool Transform;
  std::string Filter;

public:
  LibmTailCutoffPass(bool Transform, llvm::StringRef Filter)
      : Transform(Transform), Filter(Filter.str()) {}

  llvm::PreservedAnalyses run(llvm::Function &F,
                              llvm::FunctionAnalysisManager &AM);
};

} // namespace rssn
#endif
