#ifndef RSSN_LLVM_LIBM_REUSE_H
#define RSSN_LLVM_LIBM_REUSE_H

#include "llvm/IR/PassManager.h"
#include <string>

namespace rssn {

class LibmReusePass : public llvm::PassInfoMixin<LibmReusePass> {
  bool Transform;
  std::string Filter;

public:
  LibmReusePass(bool Transform, llvm::StringRef Filter)
      : Transform(Transform), Filter(Filter.str()) {}

  llvm::PreservedAnalyses run(llvm::Function &F,
                              llvm::FunctionAnalysisManager &AM);
};

} // namespace rssn
#endif
