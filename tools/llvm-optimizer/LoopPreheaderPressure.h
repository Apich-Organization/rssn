#ifndef RSSN_LLVM_LOOP_PREHEADER_PRESSURE_H
#define RSSN_LLVM_LOOP_PREHEADER_PRESSURE_H

#include "llvm/ADT/StringRef.h"
#include "llvm/IR/PassManager.h"
#include <string>

namespace rssn {

class LoopPreheaderPressurePass
    : public llvm::PassInfoMixin<LoopPreheaderPressurePass> {
  bool Transform;
  std::string Filter;

public:
  LoopPreheaderPressurePass(bool Transform, llvm::StringRef Filter)
      : Transform(Transform), Filter(Filter.str()) {}
  llvm::PreservedAnalyses run(llvm::Function &F,
                              llvm::FunctionAnalysisManager &AM);
};

} // namespace rssn
#endif
