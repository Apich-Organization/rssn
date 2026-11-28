#include "LibmUnitLog.h"
#include "llvm/ADT/FloatingPointMode.h"
#include "llvm/ADT/MapVector.h"
#include "llvm/Analysis/BlockFrequencyInfo.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/OptimizationRemarkEmitter.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Operator.h"
#include "llvm/IR/PatternMatch.h"
#include "llvm/Transforms/Utils/BasicBlockUtils.h"
#include "llvm/Transforms/Utils/Local.h"
#include <optional>

using namespace llvm;
using namespace llvm::PatternMatch;

namespace rssn {
namespace {

constexpr char PassName[] = "rssn-libm-unit-log";
constexpr unsigned MaxChains = 16;
constexpr uint64_t OneBits = UINT64_C(0x3ff0000000000000);

bool exactFP(const Instruction &I) {
  return isa<FPMathOperator>(&I) && I.getFastMathFlags().none() &&
         !I.getMetadata(LLVMContext::MD_fpmath) &&
         !I.getMetadata(LLVMContext::MD_nofpclass);
}

IntrinsicInst *pureMath(Value *V, Intrinsic::ID ID) {
  auto *I = dyn_cast_or_null<IntrinsicInst>(V);
  if (!I || I->getIntrinsicID() != ID || !I->getType()->isDoubleTy() ||
      !I->doesNotAccessMemory() || !I->hasFnAttr(Attribute::Speculatable) ||
      I->hasFnAttr(Attribute::StrictFP) || I->isConvergent() ||
      I->hasOperandBundles() || I->isMustTailCall() ||
      I->hasRetAttr(Attribute::NoFPClass) || I->paramHasAttr(0, Attribute::NoFPClass) ||
      !exactFP(*I))
    return nullptr;
  return I;
}

Instruction *uniqueNext(Instruction &I) {
  if (!I.hasOneUse())
    return nullptr;
  Instruction *Next = I.getNextNode();
  while (isa_and_nonnull<DbgInfoIntrinsic>(Next))
    Next = Next->getNextNode();
  return Next == *I.user_begin() ? Next : nullptr;
}

struct Chain {
  IntrinsicInst *Log;
  Instruction *Multiply;
  Instruction *Subtract;
  IntrinsicInst *Exp;
  Value *Node;
  Value *Coefficient;
};

std::optional<Chain> matchChain(Instruction &Start, const Loop &L) {
  auto *Log = pureMath(&Start, Intrinsic::log);
  if (!Log)
    return std::nullopt;
  Value *Node = Log->getArgOperand(0);
  Value *Offset = nullptr;
  const auto *Add = dyn_cast<Instruction>(Node);
  if (!Add || !exactFP(*Add) || !match(Add, m_c_FAdd(m_SpecificFP(1.0), m_Value(Offset))) ||
      !pureMath(Offset, Intrinsic::exp))
    return std::nullopt;

  Instruction *Mul = uniqueNext(*Log);
  Value *C = nullptr;
  if (!Mul || !exactFP(*Mul) || !match(Mul, m_c_FMul(m_Specific(Log), m_Value(C))) ||
      !L.isLoopInvariant(C) || L.isLoopInvariant(Node))
    return std::nullopt;
  Instruction *Sub = uniqueNext(*Mul);
  if (!Sub || !exactFP(*Sub) || !match(Sub, m_FSub(m_Specific(Mul), m_Specific(Node))))
    return std::nullopt;
  IntrinsicInst *Exp = pureMath(uniqueNext(*Sub), Intrinsic::exp);
  if (!Exp || Exp->use_empty())
    return std::nullopt;
  return Chain{Log, Mul, Sub, Exp, Node, C};
}

bool excludesUnit(const Chain &C, const DataLayout &DL) {
  const auto *One = ConstantFP::get(C.Node->getType(), 1.0);
  if (isImpliedByDomCondition(CmpInst::FCMP_OEQ, C.Node, One, C.Log, DL) == false)
    return true;
  const auto *Pred = C.Log->getParent()->getSinglePredecessor();
  const auto *Branch = Pred ? dyn_cast<CondBrInst>(Pred->getTerminator()) : nullptr;
  const auto *Cmp = Branch ? dyn_cast<ICmpInst>(Branch->getCondition()) : nullptr;
  if (!Cmp)
    return false;
  for (const Value *V : Cmp->operand_values()) {
    const auto *Bits = dyn_cast<BitCastInst>(V);
    if (Bits && Bits->getOperand(0) == C.Node && Bits->getType()->isIntegerTy(64)) {
      const auto *Key = ConstantInt::get(Bits->getType(), OneBits);
      if (isImpliedByDomCondition(CmpInst::ICMP_EQ, Bits, Key, C.Log, DL) == false)
        return true;
    }
  }
  return false;
}

struct Group {
  BasicBlock *Setup;
  SmallVector<Chain, 4> Chains;
};

struct Cache {
  Value *Coefficient;
  Value *Finite;
  Value *Key;
  Value *Result;
};

Cache materialize(const Group &G, const DominatorTree &DT) {
  Instruction *Insert = G.Setup->getTerminator();
  IRBuilder<> B(Insert);
  Value *Coefficient = G.Chains.front().Coefficient;
  if (!isGuaranteedNotToBeUndefOrPoison(Coefficient, nullptr, Insert, &DT))
    Coefficient = B.CreateFreeze(Coefficient, "unit.log.coefficient");
  Value *Finite = B.createIsFPClass(Coefficient, fcFinite);
  Value *Zero = ConstantFP::get(B.getDoubleTy(), 0.0);
  Value *Safe = B.CreateSelect(Finite, Coefficient, Zero);
  // Retain the math intrinsic. The cloned argument has the same rounding
  // as the original expression at log(1) and uses the runtime's exp result.
  Value *AtUnit = B.CreateFSub(B.CreateFMul(Safe, Zero),
                              ConstantFP::get(B.getDoubleTy(), 1.0));
  Value *Result = B.CreateUnaryIntrinsic(Intrinsic::exp, AtUnit, nullptr,
                                         "unit.log.cached");
  Value *Key = B.CreateSelect(Finite, B.getInt64(OneBits), B.getInt64(0),
                              "unit.log.key");
  return Cache{Coefficient, Finite, Key, Result};
}

void specialize(const Chain &C, const Cache &K, Value *Skip, Value *Node) {
  BasicBlock *Head = C.Log->getParent();
  IRBuilder<> B(C.Log);
  Instruction *Jump = SplitBlockAndInsertIfThen(B.CreateNot(Skip), C.Log, false);
  Jump->getParent()->setName("unit.log.evaluate");
  C.Log->setArgOperand(0, Node);
  C.Subtract->setOperand(1, Node);
  C.Multiply->replaceUsesOfWith(C.Coefficient, K.Coefficient);
  auto *Result = PHINode::Create(B.getDoubleTy(), 2, "unit.log.result",
                                 C.Log->getIterator());
  Result->setDebugLoc(C.Exp->getDebugLoc());
  C.Exp->replaceAllUsesWith(Result);
  for (Instruction *I : {cast<Instruction>(C.Log), C.Multiply, C.Subtract,
                         cast<Instruction>(C.Exp)}) {
    salvageDebugInfo(*I);
    I->moveBefore(Jump->getIterator());
  }
  Result->addIncoming(K.Result, Head);
  Result->addIncoming(C.Exp, Jump->getParent());
}

struct Rewrite {
  const Chain *Site;
  Cache Values;
  Value *Skip;
  Value *Node;
};

} // namespace

PreservedAnalyses LibmUnitLogPass::run(Function &F, FunctionAnalysisManager &AM) {
  if (!F.getName().contains(Filter) || F.isDeclaration() || F.hasOptNone() ||
      F.hasOptSize() || F.hasMinSize() || F.hasFnAttribute(Attribute::StrictFP) ||
      F.hasGC() ||
      F.getDenormalMode(APFloat::IEEEdouble()) != DenormalMode::getIEEE() ||
      F.size() > 1024 || F.getInstructionCount() > 24000)
    return PreservedAnalyses::all();
  const auto &LI = AM.getResult<LoopAnalysis>(F);
  if (LI.empty())
    return PreservedAnalyses::all();
  const auto &DT = AM.getResult<DominatorTreeAnalysis>(F);
  const DataLayout &DL = F.getDataLayout();
  MapVector<Value *, Group> Groups;
  unsigned Count = 0;
  for (BasicBlock &B : F) {
    const Loop *L = LI.getLoopFor(&B);
    BasicBlock *Preheader = L ? L->getLoopPreheader() : nullptr;
    if (!Preheader || B.isEHPad() || !DT.isReachableFromEntry(&B))
      continue;
    for (Instruction &I : B) {
      auto C = matchChain(I, *L);
      if (!C)
        continue;
      // A source-level or earlier IR guard can already handle the unit node.
      // Its remaining log/exp chain cannot benefit from another cache.
      if (excludesUnit(*C, DL))
        continue;
      auto [G, New] = Groups.try_emplace(C->Coefficient, Group{Preheader, {}});
      if (!New)
        G->second.Setup = DT.findNearestCommonDominator(G->second.Setup, Preheader);
      G->second.Chains.push_back(*C);
      ++Count;
      break;
    }
    if (Count == MaxChains)
      break;
  }
  if (Groups.empty())
    return PreservedAnalyses::all();
  const auto &BFI = AM.getResult<BlockFrequencyAnalysis>(F);
  auto &ORE = AM.getResult<OptimizationRemarkEmitterAnalysis>(F);
  SmallVector<Rewrite, MaxChains> Rewrites;
  for (const auto &Entry : Groups) {
    Value *Coefficient = Entry.first;
    const auto &G = Entry.second;
    if (G.Chains.size() < 2 || !G.Setup || G.Setup->isEHPad() ||
        !DT.dominates(Coefficient, G.Setup->getTerminator()))
      continue;
    const double SetupFreq = BFI.getBlockFreq(G.Setup).getFrequency();
    double Visits = 0;
    for (const Chain &C : G.Chains)
      Visits += BFI.getBlockFreq(C.Log->getParent()).getFrequency();
    if (Visits < 8 * SetupFreq)
      continue;
    ORE.emit([&] {
      using NV = ore::NV;
      return OptimizationRemarkAnalysis(PassName, "UnitLogChain", G.Setup->getTerminator())
             << NV("Function", F.getName()) << " "
             << NV("Action", Transform ? "cached" : "can cache") << " a unit-node result for "
             << NV("Chains", G.Chains.size()) << " log/exp chains";
    });
    if (!Transform)
      continue;
    Cache K = materialize(G, DT);
    for (const Chain &C : G.Chains) {
      IRBuilder<> B(C.Log);
      Value *Node = C.Node;
      Value *Skip;
      if (isGuaranteedNotToBeUndefOrPoison(Node, nullptr, C.Log, &DT)) {
        // The matched 1 + exp() is >= 1 or a NaN. Zero cannot match it
        // when the coefficient is nonfinite and the cache is disabled.
        Skip = B.CreateICmpEQ(B.CreateBitCast(Node, B.getInt64Ty()), K.Key,
                              "unit.log.match");
      } else {
        // A frozen undef can have any bits, so this branch cannot rely
        // on the impossible sentinel from the defined case above.
        Node = B.CreateFreeze(Node, "unit.log.node");
        Skip = B.CreateAnd(K.Finite, B.CreateFCmpOEQ(
                                        Node, ConstantFP::get(B.getDoubleTy(), 1.0)),
                           "unit.log.match");
      }
      Rewrites.push_back({&C, K, Skip, Node});
    }
  }
  for (const Rewrite &R : Rewrites)
    specialize(*R.Site, R.Values, R.Skip, R.Node);
  return Rewrites.empty() ? PreservedAnalyses::all() : PreservedAnalyses::none();
}

} // namespace rssn
