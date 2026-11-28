#include "LibmTailCutoff.h"
#include "FPClassInterval.h"
#include "llvm/ADT/FloatingPointMode.h"
#include "llvm/Analysis/BlockFrequencyInfo.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/OptimizationRemarkEmitter.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/PatternMatch.h"
#include "llvm/Transforms/Utils/BasicBlockUtils.h"
#include "llvm/Transforms/Utils/Local.h"
#include <optional>

using namespace llvm;
using namespace llvm::PatternMatch;

namespace rssn {
namespace {

constexpr char PassName[] = "rssn-libm-tail-cutoff";
constexpr unsigned MaxSites = 16;
constexpr double ExpZeroCutoff = 746;

bool exactFP(const Instruction &I) {
  return isa<FPMathOperator>(&I) && I.getFastMathFlags().none() &&
         !I.getMetadata(LLVMContext::MD_fpmath) &&
         !I.getMetadata(LLVMContext::MD_nofpclass);
}

IntrinsicInst *pureMath(Value *V, Intrinsic::ID ID) {
  auto *I = dyn_cast_or_null<IntrinsicInst>(V);
  if (!I || I->getIntrinsicID() != ID || !I->getType()->isDoubleTy() ||
      !I->doesNotAccessMemory() || !I->hasFnAttr(Attribute::Speculatable) ||
      I->hasFnAttr(Attribute::StrictFP) || I->isConvergent() || I->isMustTailCall() ||
      I->hasOperandBundles() || I->hasRetAttr(Attribute::NoFPClass) ||
      I->paramHasAttr(0, Attribute::NoFPClass) || !exactFP(*I))
    return nullptr;
  return I;
}

Instruction *nextDebugFree(Instruction &I) {
  Instruction *Next = I.getNextNode();
  while (isa_and_nonnull<DbgInfoIntrinsic>(Next))
    Next = Next->getNextNode();
  return Next;
}

struct Shape {
  Value *Node;
  Value *Coefficient;
  Value *Scale; // Null when the subtraction uses the node directly.
  Value *Offset; // Optional second subtraction after the linear decay.
  IntrinsicInst *Log;
  BinaryOperator *Product;
  BinaryOperator *Decay;
  BinaryOperator *Subtract;
  BinaryOperator *Shift;
};

std::optional<Shape> parseScore(Value *V, const Loop &L) {
  auto *Sub = dyn_cast<BinaryOperator>(V);
  if (!Sub || Sub->getOpcode() != Instruction::FSub || !exactFP(*Sub))
    return std::nullopt;
  Value *Offset = nullptr;
  BinaryOperator *Shift = nullptr;
  if (auto *Inner = dyn_cast<BinaryOperator>(Sub->getOperand(0));
      Inner && Inner->getOpcode() == Instruction::FSub &&
      L.isLoopInvariant(Sub->getOperand(1))) {
    if (!exactFP(*Inner) || Inner->getParent() != Sub->getParent())
      return std::nullopt;
    Shift = Sub;
    Offset = Sub->getOperand(1);
    Sub = Inner;
  }
  auto *Product = dyn_cast<BinaryOperator>(Sub->getOperand(0));
  if (!Product || !exactFP(*Product) || Product->getOpcode() != Instruction::FMul)
    return std::nullopt;
  auto *Log = pureMath(Product->getOperand(0), Intrinsic::log);
  Value *Coefficient = Product->getOperand(1);
  if (!Log) {
    Log = pureMath(Product->getOperand(1), Intrinsic::log);
    Coefficient = Product->getOperand(0);
  }
  if (!Log || !L.isLoopInvariant(Coefficient))
    return std::nullopt;
  Value *Node = Log->getArgOperand(0), *Scale = nullptr;
  auto *Decay = dyn_cast<BinaryOperator>(Sub->getOperand(1));
  if (Sub->getOperand(1) != Node) {
    if (!Decay || !exactFP(*Decay) ||
        !match(Decay, m_c_FMul(m_Specific(Node), m_Value(Scale))) ||
        !L.isLoopInvariant(Scale))
      return std::nullopt;
  } else {
    Decay = nullptr;
  }
  BasicBlock *B = Sub->getParent();
  if (Log->getParent() != B || Product->getParent() != B ||
      (Decay && Decay->getParent() != B) || L.isLoopInvariant(Node))
    return std::nullopt;
  return Shape{Node, Coefficient, Scale, Offset, Log, Product, Decay, Sub, Shift};
}

// A valid weight proves the same, defined raw exp value was finite and
// strictly positive. Neither zero, infinity nor NaN times any FP operand
// can round to a finite nonzero weight. Adding one preserves finiteness
// here, even at the largest positive finite binary64 value.
bool validNode(const Shape &S, BasicBlock &At, const DominatorTree &DT) {
  Value *Raw = S.Node;
  if (auto *Add = dyn_cast<Instruction>(S.Node);
      Add && exactFP(*Add) && match(Add, m_c_FAdd(m_SpecificFP(1.0), m_Value(Raw)))) {
    if (!pureMath(Raw, Intrinsic::exp))
      return false;
  } else if (!pureMath(Raw, Intrinsic::exp)) {
    return false;
  }
  if (!isGuaranteedNotToBeUndefOrPoison(Raw, nullptr, At.getTerminator(), &DT) ||
      Raw->getNumUses() > 16)
    return false;
  constexpr unsigned Invalid = fcNan | fcInf | fcZero;
  for (User *U : Raw->users()) {
    auto *Weight = dyn_cast<BinaryOperator>(U);
    if (!Weight || Weight->getOpcode() != Instruction::FMul ||
        !exactFP(*Weight) || Weight->getNumUses() > 16)
      continue;
    for (User *TestUse : Weight->users()) {
      Instruction *Test = nullptr;
      if (auto *Class = dyn_cast<IntrinsicInst>(TestUse);
          Class && Class->getIntrinsicID() == Intrinsic::is_fpclass &&
          Class->getArgOperand(0) == Weight &&
          (cast<ConstantInt>(Class->getArgOperand(1))->getZExtValue() & Invalid) == Invalid)
        Test = Class;
      else if (auto *Bits = dyn_cast<BitCastInst>(TestUse))
        Test = findNonfiniteOrZeroTest(*Bits);
      if (!Test || !Test->hasOneUse())
        continue;
      auto *Branch = dyn_cast<CondBrInst>(*Test->user_begin());
      if (Branch && Branch->getCondition() == Test &&
          Branch->getSuccessor(0) != Branch->getSuccessor(1) &&
          DT.dominates(BasicBlockEdge(Branch->getParent(), Branch->getSuccessor(1)), &At))
        return true;
    }
  }
  return false;
}

struct Site {
  Shape Score;
  double Limit;
  BasicBlock *Preheader;
  CondBrInst *NodeGuard;
  CondBrInst *ScoreGuard;
  IntrinsicInst *ZeroableExp;
};

bool isPositiveTest(Value *V, Value *Node) {
  auto *Cmp = dyn_cast<FCmpInst>(V);
  if (!Cmp || !exactFP(*Cmp))
    return false;
  const bool Swapped = Cmp->getOperand(1) == Node;
  const auto Predicate = Swapped ? Cmp->getSwappedPredicate() : Cmp->getPredicate();
  return (Predicate == FCmpInst::FCMP_OGT || Predicate == FCmpInst::FCMP_UGT) &&
         Cmp->getOperand(Swapped ? 1 : 0) == Node &&
         match(Cmp->getOperand(Swapped ? 0 : 1), m_SpecificFP(0.0));
}

std::optional<Site> existingCutoff(CondBrInst &Branch, const Loop &L,
                                   const DominatorTree &DT) {
  auto *Cmp = dyn_cast<FCmpInst>(Branch.getCondition());
  if (!Cmp || !exactFP(*Cmp))
    return std::nullopt;
  const bool Swapped = Cmp->getPredicate() == FCmpInst::FCMP_OGT;
  if (!Swapped && Cmp->getPredicate() != FCmpInst::FCMP_OLT)
    return std::nullopt;
  auto *Bound = dyn_cast<ConstantFP>(Cmp->getOperand(Swapped ? 0 : 1));
  if (!Bound || !Bound->getType()->isDoubleTy())
    return std::nullopt;
  const double Limit = -Bound->getValueAPF().convertToDouble();
  if (!(Limit >= 8 && Limit <= 1024))
    return std::nullopt;
  auto S = parseScore(Cmp->getOperand(Swapped ? 1 : 0), L);
  auto *B = Branch.getParent();
  auto *Prev = B->getSinglePredecessor();
  auto *Gate = Prev ? dyn_cast<CondBrInst>(Prev->getTerminator()) : nullptr;
  auto *Failure = Branch.getSuccessor(0);
  if (!S || !S->Scale || !Gate || Branch.getSuccessor(1) == Failure ||
      Gate->getSuccessor(0) != B || Gate->getSuccessor(1) != Failure ||
      !isPositiveTest(Gate->getCondition(), S->Node) || !validNode(*S, *Prev, DT))
    return std::nullopt;
  if (!all_of(Failure->phis(), [&](const PHINode &Phi) {
        return Phi.getIncomingValueForBlock(Prev) == Phi.getIncomingValueForBlock(B);
      }) || !all_of(*B, [](const Instruction &I) {
        return I.isTerminator() || isa<DbgInfoIntrinsic>(I) ||
               (!I.mayReadOrWriteMemory() && isSafeToSpeculativelyExecute(&I));
      }))
    return std::nullopt;
  return Site{*S, Limit, L.getLoopPreheader(), Gate, &Branch, nullptr};
}

std::optional<Site> zeroNumerator(Instruction &I, const Loop &L,
                                  const DominatorTree &DT) {
  IntrinsicInst *Exp = pureMath(&I, Intrinsic::exp);
  if (!Exp || Exp->use_empty())
    return std::nullopt;
  auto S = parseScore(Exp->getArgOperand(0), L);
  if (!S || S->Subtract->getParent() != Exp->getParent() ||
      !validNode(*S, *Exp->getParent(), DT))
    return std::nullopt;
  SmallVector<Instruction *, 6> Corridor{S->Log, S->Product};
  if (S->Decay)
    Corridor.push_back(S->Decay);
  Corridor.push_back(S->Subtract);
  if (S->Shift)
    Corridor.push_back(S->Shift);
  Corridor.push_back(Exp);
  for (unsigned J = 0; J + 1 < Corridor.size(); ++J)
    if (!Corridor[J]->hasOneUse() || nextDebugFree(*Corridor[J]) != Corridor[J + 1])
      return std::nullopt;
  return Site{*S, ExpZeroCutoff, L.getLoopPreheader(), nullptr, nullptr, Exp};
}

struct Group {
  Value *Coefficient;
  Value *Scale;
  Value *Offset;
  double Limit;
  BasicBlock *Setup;
  SmallVector<Site, 4> Sites;
};

Value *defined(IRBuilder<> &B, Value *V, const DominatorTree &DT) {
  if (isGuaranteedNotToBeUndefOrPoison(V, nullptr, &*B.GetInsertPoint(), &DT))
    return V;
  return B.CreateFreeze(V, "exp.tail.defined");
}

Value *buildCutoff(const Group &G, const DominatorTree &DT) {
  IRBuilder<> B(G.Setup->getTerminator());
  Type *Double = B.getDoubleTy();
  Value *Zero = ConstantFP::get(Double, 0.0);
  Value *One = ConstantFP::get(Double, 1.0);
  Value *C = defined(B, G.Coefficient, DT);
  Value *A = G.Scale ? defined(B, G.Scale, DT) : One;
  Value *Offset = G.Offset ? defined(B, G.Offset, DT) : Zero;
  Value *Valid = B.CreateAnd(B.createIsFPClass(C, fcFinite),
                             B.createIsFPClass(A, fcPosNormal | fcPosSubnormal));
  Value *PositiveC = B.CreateSelect(B.CreateFCmpOGT(C, Zero), C, Zero);
  // For finite x >= 1, log(x) is in [0, 1024]. Its rounded product with C
  // is at most 1024*max(C, 0). Subtracting a negative offset can add at most
  // -min(offset, 0) to the score. Double both bounds and add a unit of
  // margin to cover rounding at setup, decay and both scored subtractions.
  Value *Margin = B.CreateFMul(PositiveC, ConstantFP::get(Double, 2048.0));
  if (G.Offset) {
    Valid = B.CreateAnd(Valid, B.createIsFPClass(Offset, fcFinite));
    Value *OffsetIncrease =
        B.CreateSelect(B.CreateFCmpOLT(Offset, Zero), B.CreateFNeg(Offset), Zero);
    Margin = B.CreateFAdd(Margin, B.CreateFMul(OffsetIncrease,
                                                ConstantFP::get(Double, 2.0)));
  }
  Value *Target = B.CreateFAdd(Margin, ConstantFP::get(Double, G.Limit + 1.0));
  Value *Cut = G.Scale ? B.CreateFDiv(Target, A) : Target;
  Cut = B.CreateSelect(B.CreateFCmpOLT(Cut, One), One, Cut);
  for (const Site &S : G.Sites) {
    S.Score.Product->replaceUsesOfWith(G.Coefficient, C);
    if (S.Score.Decay)
      S.Score.Decay->replaceUsesOfWith(G.Scale, A);
    if (S.Score.Shift)
      S.Score.Shift->replaceUsesOfWith(G.Offset, Offset);
  }
  // A disabled or overflowing cutoff is unreachable for a proven finite x.
  return B.CreateSelect(Valid, Cut, ConstantFP::getInfinity(Double), "exp.tail.cutoff");
}

void applyCutoff(const Site &S, Value *Cutoff) {
  if (S.NodeGuard) {
    IRBuilder<> B(S.NodeGuard);
    auto *Old = cast<Instruction>(S.NodeGuard->getCondition());
    S.NodeGuard->setCondition(B.CreateFCmpOLT(S.Score.Node, Cutoff, "exp.tail.keep"));
    for (Instruction *I : {S.NodeGuard, S.ScoreGuard}) {
      I->setMetadata(LLVMContext::MD_prof, nullptr);
      I->setMetadata(LLVMContext::MD_unpredictable, nullptr);
    }
    RecursivelyDeleteTriviallyDeadInstructions(Old);
    return;
  }
  auto &C = S.Score;
  IRBuilder<> B(C.Log);
  BasicBlock *Head = C.Log->getParent();
  Value *Keep = B.CreateFCmpOLT(C.Node, Cutoff, "exp.tail.keep");
  Instruction *Jump = SplitBlockAndInsertIfThen(Keep, C.Log, false);
  Jump->getParent()->setName("exp.tail.evaluate");
  auto *Result = PHINode::Create(B.getDoubleTy(), 2, "exp.tail.result",
                                 C.Log->getIterator());
  Result->setDebugLoc(S.ZeroableExp->getDebugLoc());
  S.ZeroableExp->replaceAllUsesWith(Result);
  Instruction *Slow[] = {C.Log, C.Product, C.Decay, C.Subtract, C.Shift, S.ZeroableExp};
  for (Instruction *I : Slow) {
    if (!I)
      continue;
    salvageDebugInfo(*I);
    I->moveBefore(Jump->getIterator());
  }
  Result->addIncoming(ConstantFP::get(B.getDoubleTy(), 0.0), Head);
  Result->addIncoming(S.ZeroableExp, Jump->getParent());
}

} // namespace

PreservedAnalyses LibmTailCutoffPass::run(Function &F, FunctionAnalysisManager &AM) {
  const Triple &T = F.getParent()->getTargetTriple();
  if (T.getArch() != Triple::x86_64 ||
      (!T.isWindowsMSVCEnvironment() && (!T.isOSLinux() || T.isX32())) ||
      !F.getName().contains(Filter) || F.isDeclaration() || F.hasOptNone() ||
      F.hasOptSize() || F.hasMinSize() || F.hasFnAttribute(Attribute::StrictFP) ||
      F.hasGC() || F.getDenormalMode(APFloat::IEEEdouble()) != DenormalMode::getIEEE() ||
      F.size() > 1024 || F.getInstructionCount() > 24000)
    return PreservedAnalyses::all();
  const auto &LI = AM.getResult<LoopAnalysis>(F);
  if (LI.empty())
    return PreservedAnalyses::all();
  const auto &DT = AM.getResult<DominatorTreeAnalysis>(F);
  SmallVector<Group, 4> Groups;
  auto Add = [&](Site S) {
    auto G = find_if(Groups, [&](const Group &G) {
      return G.Coefficient == S.Score.Coefficient && G.Scale == S.Score.Scale &&
             G.Offset == S.Score.Offset && G.Limit == S.Limit;
    });
    if (G == Groups.end()) {
      Groups.push_back({S.Score.Coefficient, S.Score.Scale, S.Score.Offset,
                        S.Limit, S.Preheader, {}});
      G = std::prev(Groups.end());
    }
    if (G->Setup)
      G->Setup = DT.findNearestCommonDominator(G->Setup, S.Preheader);
    G->Sites.push_back(S);
  };
  unsigned Count = 0;
  for (BasicBlock &B : F) {
    const auto *L = LI.getLoopFor(&B);
    if (!L || !L->getLoopPreheader() || B.isEHPad() || !DT.isReachableFromEntry(&B))
      continue;
    bool Found = false;
    if (auto *Term = dyn_cast<CondBrInst>(B.getTerminator()))
      if (auto S = existingCutoff(*Term, *L, DT)) {
        Add(*S);
        Found = true;
      }
    if (!Found)
      for (Instruction &I : B)
        if (auto S = zeroNumerator(I, *L, DT)) {
          Add(*S);
          Found = true;
          break;
        }
    if (Found && ++Count == MaxSites)
      break;
  }
  if (Groups.empty())
    return PreservedAnalyses::all();
  const auto &BFI = AM.getResult<BlockFrequencyAnalysis>(F);
  auto &ORE = AM.getResult<OptimizationRemarkEmitterAnalysis>(F);
  SmallVector<std::pair<const Group *, Value *>, 4> Work;
  for (Group &G : Groups) {
    if (G.Sites.size() < 2 || !G.Setup || G.Setup->isEHPad())
      continue;
    Instruction *Insert = G.Setup->getTerminator();
    if (!DT.dominates(G.Coefficient, Insert) || (G.Scale && !DT.dominates(G.Scale, Insert)) ||
        (G.Offset && !DT.dominates(G.Offset, Insert)))
      continue;
    double Frequency = 0;
    for (const Site &S : G.Sites)
      Frequency += BFI.getBlockFreq(S.Score.Log->getParent()).getFrequency();
    if (Frequency < 8.0 * BFI.getBlockFreq(G.Setup).getFrequency())
      continue;
    ORE.emit([&] {
      using NV = ore::NV;
      return OptimizationRemarkAnalysis(PassName, "EarlyTailCutoff", Insert)
             << NV("Function", F.getName()) << " "
             << NV("Action", Transform ? "placed" : "can place") << " a shared cutoff before "
             << NV("Sites", G.Sites.size()) << " log chains";
    });
    if (Transform)
      Work.emplace_back(&G, buildCutoff(G, DT));
  }
  for (auto [G, Cutoff] : Work)
    for (const Site &S : G->Sites)
      applyCutoff(S, Cutoff);
  return Work.empty() ? PreservedAnalyses::all() : PreservedAnalyses::none();
}

} // namespace rssn
