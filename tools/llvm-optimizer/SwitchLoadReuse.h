#ifndef RSSN_LLVM_SWITCH_LOAD_REUSE_H
#define RSSN_LLVM_SWITCH_LOAD_REUSE_H

#include "llvm/ADT/StringRef.h"
#include "llvm/IR/PassManager.h"
#include <string>

namespace rssn {

class SwitchLoadReusePass : public llvm::PassInfoMixin<SwitchLoadReusePass> {
  bool Transform;
  std::string Filter;

public:
  SwitchLoadReusePass(bool Transform, llvm::StringRef Filter)
      : Transform(Transform), Filter(Filter.str()) {}
  llvm::PreservedAnalyses run(llvm::Function &F,
                              llvm::FunctionAnalysisManager &AM);
};

} // namespace rssn
#endif
