#ifndef RSSN_LLVM_SCALAR_LOAD_PACK_H
#define RSSN_LLVM_SCALAR_LOAD_PACK_H

#include "llvm/IR/PassManager.h"
#include <string>

namespace rssn {

class ScalarLoadPackPass : public llvm::PassInfoMixin<ScalarLoadPackPass> {
  bool Transform;
  std::string Filter;

public:
  ScalarLoadPackPass(bool Transform, llvm::StringRef Filter)
      : Transform(Transform), Filter(Filter.str()) {}

  llvm::PreservedAnalyses run(llvm::Function &F,
                              llvm::FunctionAnalysisManager &AM);
};

} // namespace rssn

#endif
