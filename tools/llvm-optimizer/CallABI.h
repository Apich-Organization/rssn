#ifndef RSSN_LLVM_SPILL_CALL_ABI_H
#define RSSN_LLVM_SPILL_CALL_ABI_H

#include "llvm/IR/PassManager.h"
#include <string>

namespace rssn {

class CallABIPass : public llvm::PassInfoMixin<CallABIPass> {
  bool Transform;
  std::string Filter;

public:
  CallABIPass(bool Transform, llvm::StringRef Filter)
      : Transform(Transform), Filter(Filter.str()) {}

  llvm::PreservedAnalyses run(llvm::Module &M, llvm::ModuleAnalysisManager &AM);
};

} // namespace rssn
#endif
