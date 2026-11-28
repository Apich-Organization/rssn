#ifndef RSSN_LLVM_FP_REGISTER_BANK_H
#define RSSN_LLVM_FP_REGISTER_BANK_H

#include "llvm/IR/PassManager.h"
#include <string>

namespace rssn {

class FPRegisterBankPass : public llvm::PassInfoMixin<FPRegisterBankPass> {
  bool Transform;
  std::string Filter;

public:
  FPRegisterBankPass(bool Transform, llvm::StringRef Filter)
      : Transform(Transform), Filter(Filter.str()) {}

  llvm::PreservedAnalyses run(llvm::Function &F,
                              llvm::FunctionAnalysisManager &AM);
};

} // namespace rssn
#endif
