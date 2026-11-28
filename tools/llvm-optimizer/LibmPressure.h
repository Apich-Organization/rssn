#ifndef RSSN_LLVM_LIBM_PRESSURE_H
#define RSSN_LLVM_LIBM_PRESSURE_H

#include "llvm/IR/PassManager.h"
#include <string>

namespace llvm {
class Instruction;
}

namespace rssn {

bool isScalarLibmCut(const llvm::Instruction &I);

class LibmPressurePass : public llvm::PassInfoMixin<LibmPressurePass> {
  bool Transform;
  bool TrialCones;
  std::string Filter;

public:
  LibmPressurePass(bool Transform, llvm::StringRef Filter,
                   bool TrialCones = false)
      : Transform(Transform), TrialCones(TrialCones), Filter(Filter.str()) {}

  llvm::PreservedAnalyses run(llvm::Function &F,
                              llvm::FunctionAnalysisManager &AM);
};

} // namespace rssn
#endif
