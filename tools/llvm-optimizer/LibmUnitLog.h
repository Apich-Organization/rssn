#ifndef RSSN_LLVM_LIBM_UNIT_LOG_H
#define RSSN_LLVM_LIBM_UNIT_LOG_H

#include "llvm/IR/PassManager.h"
#include <string>

namespace rssn {

class LibmUnitLogPass : public llvm::PassInfoMixin<LibmUnitLogPass> {
  bool Transform;
  std::string Filter;

public:
  LibmUnitLogPass(bool Transform, llvm::StringRef Filter)
      : Transform(Transform), Filter(Filter.str()) {}
  llvm::PreservedAnalyses run(llvm::Function &F, llvm::FunctionAnalysisManager &AM);
};

} // namespace rssn
#endif
