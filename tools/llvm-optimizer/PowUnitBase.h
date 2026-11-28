#ifndef RSSN_LLVM_POW_UNIT_BASE_H
#define RSSN_LLVM_POW_UNIT_BASE_H

#include "llvm/IR/PassManager.h"
#include <string>

namespace rssn {

class PowUnitBasePass : public llvm::PassInfoMixin<PowUnitBasePass> {
  bool Transform;
  std::string Filter;

public:
  PowUnitBasePass(bool Transform, llvm::StringRef Filter)
      : Transform(Transform), Filter(Filter.str()) {}

  llvm::PreservedAnalyses run(llvm::Function &F,
                              llvm::FunctionAnalysisManager &AM);
};

} // namespace rssn
#endif
