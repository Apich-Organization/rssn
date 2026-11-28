#include "PowUnitBase.h"
#include "llvm/ADT/FloatingPointMode.h"
#include "llvm/ADT/MapVector.h"
#include "llvm/ADT/SmallPtrSet.h"
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

constexpr char PassName[] = "rssn-pow-unit-base";
constexpr unsigned MinPowers = 3;
constexpr unsigned MaxPowers = 4;
constexpr unsigned MaxChains = 16;

bool exactFP(const Instruction &I) {
  return isa<FPMathOperator>(&I) && I.getFastMathFlags().none() &&
         !I.getMetadata(LLVMContext::MD_fpmath) &&
         !I.getMetadata(LLVMContext::MD_nofpclass);
}

IntrinsicInst *purePow(Value *V) {
  auto *I = dyn_cast<IntrinsicInst>(V);
  if (!I || I->getIntrinsicID() != Intrinsic::pow ||
      !I->getType()->isDoubleTy() || !I->hasOneUse() ||
      !I->doesNotAccessMemory() || !I->hasFnAttr(Attribute::Speculatable) ||
      I->hasFnAttr(Attribute::StrictFP) || I->hasOperandBundles() ||
      I->isMustTailCall() || I->isConvergent() ||
      I->hasRetAttr(Attribute::NoFPClass) ||
      I->paramHasAttr(0, Attribute::NoFPClass) ||
      I->paramHasAttr(1, Attribute::NoFPClass) || !exactFP(*I))
    return nullptr;
  return I;
}

struct Factor {
  IntrinsicInst *Pow;
  BinaryOperator *Base;
  BinaryOperator *Scaled;
  Value *Coefficient;
  BinaryOperator *Multiply;
};

struct Chain {
  Value *Node;
  Value *Initial;
  BasicBlock *Preheader;
  SmallVector<Factor, MaxPowers> Factors;
  SmallVector<Instruction *, 4 * MaxPowers> Region;
  IntrinsicInst *PrefixPow = nullptr;
  Value *GuardedNode = nullptr;
};

BinaryOperator *unitMinus(IntrinsicInst &Pow, Value *&RHS) {
  auto *Base = dyn_cast<BinaryOperator>(Pow.getArgOperand(0));
  if (!Base || !exactFP(*Base) ||
      !match(Base, m_FSub(m_SpecificFP(1.0), m_Value(RHS))))
    return nullptr;
  return Base;
}

std::optional<Chain> matchChain(IntrinsicInst &First, const Loop &L,
                                const DominatorTree &DT) {
  BasicBlock *B = First.getParent();
  BasicBlock *Preheader = L.getLoopPreheader();
  Value *Node = nullptr;
  auto *FirstBase = unitMinus(First, Node);
  if (!Preheader || !FirstBase || L.isLoopInvariant(Node))
    return std::nullopt;
  // Early CSE can share 1-t between both halves of a quadrature loop or
  // with an endpoint check. Leave such a definition in place.
  const bool MoveFirstBase =
      FirstBase->getParent() == B && FirstBase->hasOneUse();
  Instruction *Start = MoveFirstBase ? cast<Instruction>(FirstBase) : &First;

  auto *Multiply = dyn_cast<BinaryOperator>(*First.user_begin());
  Value *Initial = nullptr;
  if (!Multiply || Multiply->getParent() != B || !exactFP(*Multiply) ||
      !match(Multiply, m_c_FMul(m_Specific(&First), m_Value(Initial))) ||
      !DT.dominates(Initial, Start) ||
      !L.isLoopInvariant(First.getArgOperand(1)))
    return std::nullopt;

  Chain C{Node,
          Initial,
          Preheader,
          {{&First, MoveFirstBase ? FirstBase : nullptr, nullptr, nullptr,
            Multiply}},
          {}};
  while (C.Factors.size() < MaxPowers && Multiply->hasOneUse()) {
    auto *Next = dyn_cast<BinaryOperator>(*Multiply->user_begin());
    Value *Other = nullptr;
    if (!Next || Next->getParent() != B || !exactFP(*Next) ||
        !match(Next, m_c_FMul(m_Specific(Multiply), m_Value(Other))))
      break;
    auto *Pow = purePow(Other);
    if (!Pow || Pow->getParent() != B || !Multiply->comesBefore(Pow) ||
        !L.isLoopInvariant(Pow->getArgOperand(1)))
      break;
    Value *RHS = nullptr, *Coefficient = nullptr;
    auto *Base = unitMinus(*Pow, RHS);
    auto *Scaled = dyn_cast_or_null<BinaryOperator>(RHS);
    if (!Base || Base->getParent() != B || !Base->hasOneUse() || !Scaled ||
        Scaled->getParent() != B || !Scaled->hasOneUse() || !exactFP(*Scaled) ||
        !match(Scaled, m_c_FMul(m_Specific(Node), m_Value(Coefficient))) ||
        !L.isLoopInvariant(Coefficient))
      break;
    C.Factors.push_back({Pow, Base, Scaled, Coefficient, Next});
    Multiply = Next;
  }
  if (C.Factors.size() < MinPowers)
    return std::nullopt;

  // The slow branch must be a contiguous, side-effect-free region and only
  // the final accumulator may escape it.
  SmallPtrSet<Instruction *, 16> Expected;
  for (const Factor &F : C.Factors) {
    if (F.Base)
      Expected.insert(F.Base);
    Expected.insert(F.Pow);
    Expected.insert(F.Multiply);
    if (F.Scaled)
      Expected.insert(F.Scaled);
  }
  for (Instruction *I = Start; I; I = I->getNextNode()) {
    if (Expected.contains(I))
      C.Region.push_back(I);
    else if (!isa<DbgInfoIntrinsic>(I))
      return std::nullopt;
    if (I == Multiply)
      break;
  }
  if (C.Region.size() != Expected.size())
    return std::nullopt;
  if (auto *Prefix = purePow(Initial);
      Prefix && Prefix->getArgOperand(0) == Node &&
      L.isLoopInvariant(Prefix->getArgOperand(1)))
    C.PrefixPow = Prefix;
  return C;
}

using Parameters = SmallVector<Value *, 2 * MaxPowers>;

Parameters setupInputs(const Chain &C) {
  Parameters Inputs;
  for (const Factor &F : C.Factors) {
    Inputs.push_back(F.Pow->getArgOperand(1));
    if (F.Coefficient)
      Inputs.push_back(F.Coefficient);
  }
  if (C.PrefixPow)
    Inputs.push_back(C.PrefixPow->getArgOperand(1));
  return Inputs;
}

struct Group {
  Parameters Inputs;
  SmallVector<Chain, 4> Chains;
  BasicBlock *Setup;
};

Value *defined(IRBuilder<> &B, Value *V, const DominatorTree &DT) {
  if (isGuaranteedNotToBeUndefOrPoison(V, nullptr, &*B.GetInsertPoint(), &DT))
    return V;
  return B.CreateFreeze(V, "unit.defined");
}

Value *buildCutoff(const Group &G, MapVector<Value *, Value *> &Frozen,
                   const DominatorTree &DT) {
  IRBuilder<> B(G.Setup->getTerminator());
  Value *Largest = B.getInt64(1023);
  Value *Ordered = B.getTrue();
  SmallPtrSet<Value *, MaxPowers> Scaled, Exponents;
  for (Value *V : G.Inputs)
    if (!Frozen.count(V))
      Frozen[V] = defined(B, V, DT);
  auto RequireNumber = [&](Value *Exp) {
    if (Exponents.insert(Exp).second) {
      Value *V = Frozen[Exp];
      // Keep NaN exponents on the original libm path, including sNaNs on
      // runtimes that return a NaN for pow(1, sNaN).
      Ordered = B.CreateAnd(Ordered, B.CreateFCmpORD(V, V));
    }
  };
  for (const Factor &F : G.Chains.front().Factors) {
    if (F.Coefficient && Scaled.insert(F.Coefficient).second) {
      Value *Bits = B.CreateBitCast(Frozen[F.Coefficient], B.getInt64Ty());
      Value *Exponent = B.CreateAnd(B.CreateLShr(Bits, 52), B.getInt64(2047));
      Largest = B.CreateBinaryIntrinsic(Intrinsic::umax, Largest, Exponent);
    }
    RequireNumber(F.Pow->getArgOperand(1));
  }
  if (auto *Prefix = G.Chains.front().PrefixPow)
    RequireNumber(Prefix->getArgOperand(1));

  // If the greatest biased coefficient exponent is E, every |coefficient|
  // is less than 2^(E-1022), including the implicit coefficient of one.
  // For 0 <= node <= 2^(968-E), |node * coefficient| < 2^-54. The
  // rounded product can reach 2^-54, but 1 - 2^-54 still rounds to 1
  // (ties to even). The negative product has an even wider rounding bin.
  // E >= 1991 would make the encoded bound non-normal and also includes
  // infinities and NaNs. Disable the fast path in this case.
  Value *Exponent = B.CreateSub(B.getInt64(1991), Largest);
  Value *Bits = B.CreateShl(Exponent, 52);
  Value *Exclusive = B.CreateAdd(Bits, B.getInt64(1));
  Value *Valid =
      B.CreateAnd(Ordered, B.CreateICmpULT(Largest, B.getInt64(1991)));
  return B.CreateSelect(Valid, Exclusive, B.getInt64(0), "pow.unit.cutoff");
}

void specialize(const Chain &C, Value *Cutoff,
                const MapVector<Value *, Value *> &Frozen) {
  Instruction *Start = C.Region.front();
  IRBuilder<> Check(Start);
  Value *Node = C.GuardedNode;
  Value *Bits = Check.CreateBitCast(Node, Check.getInt64Ty());
  Value *Skip = Check.CreateICmpULT(Bits, Cutoff, "pow.unit");
  // Multiplication by one may disappear during lowering. Keep a NaN
  // accumulator on the original path so its quieting/payload behavior stays
  // with the original operations. A prefix pow with a nonnegative finite base
  // and a non-NaN exponent is non-NaN here, without a test in the loop.
  Value *Initial = C.Initial;
  if (auto *Prefix = C.PrefixPow) {
    Prefix->setArgOperand(0, Node);
    Prefix->setArgOperand(1, Frozen.lookup(Prefix->getArgOperand(1)));
  } else {
    Initial = Check.CreateFreeze(Initial, "pow.unit.initial");
    Skip = Check.CreateAnd(Skip, Check.CreateFCmpORD(Initial, Initial));
    C.Factors.front().Multiply->replaceUsesOfWith(C.Initial, Initial);
  }
  Instruction *FastEnd = nullptr, *SlowEnd = nullptr;
  SplitBlockAndInsertIfThenElse(Skip, Start, &FastEnd, &SlowEnd);
  FastEnd->getParent()->setName("pow.unit.fast");
  SlowEnd->getParent()->setName("pow.unit.evaluate");

  IRBuilder<> Fast(FastEnd);
  Value *One = ConstantFP::get(Check.getDoubleTy(), 1.0);
  Value *Product = Initial;
  for (const Factor &F : C.Factors) {
    // Keep the original multiplication order and operand positions.
    if (F.Multiply->getOperand(0) == F.Pow)
      Product = Fast.CreateFMul(One, Product);
    else
      Product = Fast.CreateFMul(Product, One);

    if (F.Scaled)
      F.Scaled->replaceUsesOfWith(F.Coefficient, Frozen.lookup(F.Coefficient));
    Instruction *User = F.Scaled ? F.Scaled : F.Base;
    if (User)
      User->replaceUsesOfWith(C.Node, Node);
    F.Pow->setArgOperand(1, Frozen.lookup(F.Pow->getArgOperand(1)));
  }

  auto *Last = C.Factors.back().Multiply;
  auto *Result = PHINode::Create(Last->getType(), 2, "pow.unit.result",
                                 Start->getIterator());
  Result->setDebugLoc(Last->getDebugLoc());
  Last->replaceAllUsesWith(Result);
  for (Instruction *I : C.Region) {
    salvageDebugInfo(*I);
    I->moveBefore(SlowEnd->getIterator());
  }
  const Factor &First = C.Factors.front();
  if (!First.Base && Node != C.Node) {
    IRBuilder<> AtPow(First.Pow);
    First.Pow->setArgOperand(0, AtPow.CreateFSub(One, Node));
  }
  Result->addIncoming(Product, FastEnd->getParent());
  Result->addIncoming(Last, SlowEnd->getParent());
}

} // namespace

PreservedAnalyses PowUnitBasePass::run(Function &F,
                                       FunctionAnalysisManager &AM) {
  if (!F.getName().contains(Filter) || F.isDeclaration() || F.hasOptNone() ||
      F.hasOptSize() || F.hasMinSize() ||
      F.hasFnAttribute(Attribute::StrictFP) ||
      F.getDenormalMode(APFloat::IEEEdouble()) != DenormalMode::getIEEE() ||
      F.size() > 1024 || F.getInstructionCount() > 24000)
    return PreservedAnalyses::all();
  const auto &LI = AM.getResult<LoopAnalysis>(F);
  if (LI.empty())
    return PreservedAnalyses::all();
  const auto &DT = AM.getResult<DominatorTreeAnalysis>(F);
  SmallVector<Group, 4> Groups;
  unsigned Count = 0;
  for (BasicBlock &B : F) {
    const Loop *L = LI.getLoopFor(&B);
    if (!L || B.isEHPad() || !DT.isReachableFromEntry(&B))
      continue;
    for (Instruction &I : B) {
      auto *Pow = purePow(&I);
      if (!Pow)
        continue;
      auto Candidate = matchChain(*Pow, *L, DT);
      if (!Candidate)
        continue;
      auto Inputs = setupInputs(*Candidate);
      auto G = find_if(
          Groups, [&](const Group &Group) { return Group.Inputs == Inputs; });
      if (G == Groups.end()) {
        BasicBlock *Preheader = Candidate->Preheader;
        Groups.push_back({std::move(Inputs), {}, Preheader});
        G = std::prev(Groups.end());
      }
      G->Setup = DT.findNearestCommonDominator(G->Setup, Candidate->Preheader);
      G->Chains.push_back(std::move(*Candidate));
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
  SmallVector<Group *, 4> Selected;
  for (Group &G : Groups) {
    if (!G.Setup || G.Setup->isEHPad() || G.Chains.size() < 2 ||
        !all_of(G.Inputs, [&](Value *V) {
          return DT.dominates(V, G.Setup->getTerminator());
        }))
      continue;
    const double SetupFreq = BFI.getBlockFreq(G.Setup).getFrequency();
    double Visits = 0;
    for (const Chain &C : G.Chains)
      Visits += BFI.getBlockFreq(C.Region.front()->getParent()).getFrequency();
    if (Visits < 8 * SetupFreq)
      continue;
    ORE.emit([&] {
      using NV = ore::NV;
      return OptimizationRemarkAnalysis(PassName, "UnitPowerBases",
                                        G.Setup->getTerminator())
             << NV("Function", F.getName()) << " "
             << NV("Action", Transform ? "guarded" : "can guard") << " "
             << NV("Powers", G.Chains.front().Factors.size())
             << " powers per node in " << NV("Chains", G.Chains.size())
             << " blocks with one shared unit-base cutoff";
    });
    Selected.push_back(&G);
  }

  if (!Transform || Selected.empty())
    return PreservedAnalyses::all();
  // Compute all value-tracking facts while the analyzed CFG still exists.
  struct Rewrite {
    Group *G;
    MapVector<Value *, Value *> Inputs;
    Value *Cutoff;
  };
  SmallVector<Rewrite, 4> Rewrites;
  for (Group *G : Selected) {
    Rewrites.push_back({G, {}, nullptr});
    Rewrite &R = Rewrites.back();
    R.Cutoff = buildCutoff(*G, R.Inputs, DT);
    for (Chain &C : G->Chains) {
      IRBuilder<> B(C.PrefixPow ? cast<Instruction>(C.PrefixPow)
                                : C.Region.front());
      C.GuardedNode = defined(B, C.Node, DT);
    }
  }
  for (Rewrite &R : Rewrites)
    for (const Chain &C : R.G->Chains)
      specialize(C, R.Cutoff, R.Inputs);
  return PreservedAnalyses::none();
}

} // namespace rssn
