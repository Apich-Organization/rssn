#ifndef RSSN_LLVM_FP_CLASS_INTERVAL_H
#define RSSN_LLVM_FP_CLASS_INTERVAL_H

#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/PassManager.h"
#include <string>

namespace rssn {

llvm::Value *buildFiniteNonzeroTest(llvm::IRBuilder<> &B, llvm::Value *Input,
                                   bool Invalid);
llvm::ICmpInst *findNonfiniteOrZeroTest(llvm::BitCastInst &Bits);

class FPClassIntervalPass : public llvm::PassInfoMixin<FPClassIntervalPass> {
  bool Transform;
  std::string Filter;

public:
  FPClassIntervalPass(bool Transform, llvm::StringRef Filter)
      : Transform(Transform), Filter(Filter.str()) {}

  llvm::PreservedAnalyses run(llvm::Function &F,
                             llvm::FunctionAnalysisManager &AM);
};

} // namespace rssn
#endif
