#ifndef RSSN_LLVM_FP_ENDPOINT_GUARDS_H
#define RSSN_LLVM_FP_ENDPOINT_GUARDS_H

#include "llvm/IR/PassManager.h"
#include <string>

namespace rssn {
class FPEndpointGuardsPass : public llvm::PassInfoMixin<FPEndpointGuardsPass> {
  bool Transform;
  std::string Filter;

public:
  FPEndpointGuardsPass(bool Transform, llvm::StringRef Filter)
      : Transform(Transform), Filter(Filter.str()) {}

  llvm::PreservedAnalyses run(llvm::Function &F,
                             llvm::FunctionAnalysisManager &AM);
};
} // namespace rssn
#endif
