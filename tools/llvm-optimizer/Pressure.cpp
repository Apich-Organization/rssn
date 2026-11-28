#include "Pressure.h"
#include "LibmPressure.h"
#include "llvm/ADT/STLExtras.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/Analysis/TargetTransformInfo.h"
#include "llvm/IR/CFG.h"
#include "llvm/IR/DerivedTypes.h"
#include "llvm/IR/Function.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/Module.h"

using namespace llvm;

namespace rssn {

Pressure &Pressure::operator+=(Pressure Other) {
  GP += Other.GP;
  SIMD += Other.SIMD;
  SIMDSpill += Other.SIMDSpill;
  PreservableSIMD += Other.PreservableSIMD;
  VolatileSIMD += Other.VolatileSIMD;
  VolatileSIMDSpill += Other.VolatileSIMDSpill;
  return *this;
}

Pressure &Pressure::operator-=(Pressure Other) {
  GP -= Other.GP;
  SIMD -= Other.SIMD;
  SIMDSpill -= Other.SIMDSpill;
  PreservableSIMD -= Other.PreservableSIMD;
  VolatileSIMD -= Other.VolatileSIMD;
  VolatileSIMDSpill -= Other.VolatileSIMDSpill;
  return *this;
}

int Pressure::excess(Pressure Limit) const {
  return std::max(0, GP - Limit.GP) + std::max(0, SIMD - Limit.SIMD);
}

// Distinguish sparse pressure near the budget from long dense overlap. The
// quadratic term saturates so large legalization estimates cannot dominate.
static int bankCongestion(int Used, int Limit) {
  const int Excess = std::max(0, Used - Limit);
  const int Dense = std::min(Excess, 6);
  return std::clamp(Used - (Limit - 2), 0, 2) + 3 * Excess + Dense * Dense;
}

int Pressure::congestion(Pressure Limit) const {
  return bankCongestion(GP, Limit.GP) + bankCongestion(SIMD, Limit.SIMD);
}

int Pressure::win64CallExcess() const {
  // Upper bounds on preserved registers under the ordinary Win64 ABI.
  return std::max(0, GP - Win64PreservedGP) +
         std::max(0, PreservableSIMD - Win64PreservedXMM) + VolatileSIMD;
}

int Pressure::win64CallCost() const {
  // As in callCost, account for a store/reload pair and a width term.
  // Win64 preserves the low 128 bits of XMM6–15, never a whole YMM or ZMM.
  return 5 * (std::max(0, GP - Win64PreservedGP) +
              std::max(0, PreservableSIMD - Win64PreservedXMM)) +
         VolatileSIMDSpill;
}

int Pressure::aapcs64CallCost() const {
  // X19–X28 are allocatable nonvolatile GPRs. Only D8–D15 survive as
  // ordinary SIMD values, so a live Q register is charged its full spill.
  return 5 * (std::max(0, GP - 10) + std::max(0, PreservableSIMD - 8)) +
         VolatileSIMDSpill;
}

int Pressure::callCost(Pressure Limit) const {
  // SysV C ABI preserves no SIMD registers. A pair of L1 accesses has both
  // an instruction cost and a width-dependent transfer cost.
  return 5 * std::max(0, GP - Limit.GP) + SIMDSpill;
}

bool Pressure::worsens(Pressure Before, Pressure Limit) const {
  return GP > std::max(Before.GP, Limit.GP) ||
         SIMD > std::max(Before.SIMD, Limit.SIMD);
}

bool Liveness::supports(const Function &F) {
  if (F.size() > 1024)
    return false;
  unsigned Count = 0;
  for (const BasicBlock &B : F) {
    if (B.size() > MaxBlockSize)
      return false;
    Count += B.size();
    if (Count > 24000)
      return false;
  }
  return Count >= 4;
}

Liveness::Liveness(const TargetTransformInfo &Target, const Function &F)
    : TTI(Target), DL(F.getParent()->getDataLayout()),
      AArch64(F.getParent()->getTargetTriple().getArch() == Triple::aarch64) {
  LLVMContext &C = F.getContext();
  auto *Vector = FixedVectorType::get(Type::getDoubleTy(C), 2);
  const unsigned NumSIMD =
      TTI.getNumberOfRegisters(TTI.getRegisterClassForType(true, Vector));
  const unsigned NumGP = TTI.getNumberOfRegisters(
      TTI.getRegisterClassForType(false, Type::getInt64Ty(C)));
  Limit = {NumGP ? std::clamp(static_cast<int>(NumGP) - 4, 1, 28) : 12,
           NumSIMD ? std::clamp(static_cast<int>(NumSIMD) - 1, 1, 31) : 15};
}

Pressure Liveness::typeWeight(Type *Ty) const {
  // X86 lowers raw byte values using the corresponding integer machine
  // types. Account for the carrier even when the IR keeps its bits opaque.
  if (Ty->isByteOrByteVectorTy())
    Ty = Type::getIntFromByteType(Ty);
  if (Ty->isIntegerTy() || Ty->isPointerTy()) {
    const unsigned Minimum =
        Ty->isIntegerTy() ? divideCeil(Ty->getIntegerBitWidth(), 64) : 1;
    const unsigned Parts = std::max(TTI.getRegUsageForType(Ty), Minimum);
    return {static_cast<int>(std::min(Parts, 128u)), 0, 0};
  }

  const auto *Vector = dyn_cast<FixedVectorType>(Ty);
  if (!Vector && !Ty->isHalfTy() && !Ty->isBFloatTy() && !Ty->isFloatTy() &&
      !Ty->isDoubleTy())
    return {};
  // X86 TTI exposes vector and scalar FP classes, but gives no reliable
  // context-free mapping from an IR i1 vector to a SIMD or mask register.
  if (Vector && Vector->getElementType()->isIntegerTy(1))
    return {};
  const unsigned Parts = std::clamp(TTI.getRegUsageForType(Ty), 1u, 128u);
  if (Vector && !Vector->getElementType()->isFloatingPointTy() &&
      Parts >= Vector->getNumElements() && !TTI.isTypeLegal(Ty))
    return {static_cast<int>(Parts), 0, 0};

  const auto Load = TTI.getMemoryOpCost(Instruction::Load, Ty, Align(16), 0);
  const auto Store = TTI.getMemoryOpCost(Instruction::Store, Ty, Align(16), 0);
  const auto Pair = Load + Store;
  const int Operations =
      Pair.isValid() ? std::clamp<int64_t>(Pair.getValue(), 2, 256) : 2 * Parts;
  const uint64_t Bytes = DL.getTypeStoreSize(Ty).getFixedValue();
  const int Transfer = std::min<uint64_t>(128, divideCeil(Bytes, 16));
  // A full-width YMM/ZMM value cannot survive a Win64 call in a single
  // preserved XMM. Keep it in the general pressure estimate.
  const int Spill = 2 * Operations + Transfer;
  const auto Registers = static_cast<int>(Parts);
  // Only a legal low-64-bit vector or a scalar FP value fits into D8–D15.
  // Wider or legalized AArch64 values have no dependable preserved carrier.
  if (AArch64) {
    const bool Fits = Bytes <= 8 && (!Vector || TTI.isTypeLegal(Ty));
    return {0,
            Registers,
            Spill,
            Fits ? Registers : 0,
            Fits ? 0 : Registers,
            Fits ? 0 : Spill};
  }
  const bool NativeWide = Bytes > 16 && TTI.isTypeLegal(Ty);
  return {0,
          Registers,
          Spill,
          Bytes <= 16 ? Registers : 0,
          NativeWide ? Registers : 0,
          NativeWide ? Spill : 0};
}

void Liveness::track(const Value &V) {
  auto [It, New] = TypeWeights.try_emplace(V.getType());
  if (New)
    It->second = typeWeight(V.getType());
  const Pressure W = It->second;
  if (W.GP || W.SIMD) {
    Values[&V] = Weights.size();
    Weights.push_back(W);
  }
}

void Liveness::rebind(const Value &From, const Value &To) {
  auto Old = Values.find(&From);
  if (Old == Values.end())
    return;
  const unsigned Slot = Old->second;
  Values.erase(Old);
  Values[&To] = Slot;
}

unsigned Liveness::id(const Value *V) const {
  auto It = Values.find(V);
  return It == Values.end() ? Weights.size() : It->second;
}

Pressure Liveness::weight(const Value *V) const {
  const unsigned ID = id(V);
  return ID == Weights.size() ? Pressure{} : Weights[ID];
}

const BitVector &Liveness::liveOut(const BasicBlock &B) const {
  return States[Blocks.lookup(&B)].Out;
}

bool Liveness::compute(Function &F, ArrayRef<Schedule> Trial) {
  for (const Argument &Arg : F.args())
    track(Arg);
  for (BasicBlock &B : F) {
    Blocks[&B] = Blocks.size();
    for (Instruction &I : B)
      track(I);
  }
  const unsigned N = Weights.size();
  if (N > 12000)
    return false;
  States.reserve(F.size());
  for (size_t I = 0; I < F.size(); ++I)
    States.emplace_back(N);

  for (BasicBlock &B : F) {
    State &S = States[Blocks.lookup(&B)];
    auto Account = [&](Instruction &I) {
      if (auto *Phi = dyn_cast<PHINode>(&I)) {
        for (unsigned K = 0; K < Phi->getNumIncomingValues(); ++K) {
          const unsigned V = id(Phi->getIncomingValue(K));
          if (V != N)
            States[Blocks.lookup(Phi->getIncomingBlock(K))].PhiUses.set(V);
        }
      } else {
        for (const Use &U : I.operands()) {
          const unsigned V = id(U.get());
          if (V != N && !S.Kill.test(V))
            S.Gen.set(V);
        }
      }
      const unsigned Def = id(&I);
      if (Def != N)
        S.Kill.set(Def);
    };
    auto Override = find_if(
        Trial, [&](const Schedule &Order) { return Order.Block == &B; });
    if (Override == Trial.end()) {
      for (Instruction &I : B)
        Account(I);
    } else {
      for (Instruction *I : Override->Code)
        Account(*I);
    }
  }

  SmallVector<BasicBlock *, 32> Work;
  SmallPtrSet<BasicBlock *, 32> Queued;
  for (BasicBlock &B : F) {
    Work.push_back(&B);
    Queued.insert(&B);
  }
  BitVector In(N), Out(N);
  unsigned Steps = 0;
  const unsigned WorkLimit = 32 * F.size() + 512;
  while (!Work.empty()) {
    if (++Steps > WorkLimit)
      return false;
    BasicBlock *B = Work.pop_back_val();
    Queued.erase(B);
    State &S = States[Blocks.lookup(B)];
    Out = S.PhiUses;
    for (BasicBlock *Succ : successors(B))
      Out |= States[Blocks.lookup(Succ)].In;
    S.Out = Out;
    In = Out;
    In.reset(S.Kill);
    In |= S.Gen;
    if (In == S.In)
      continue;
    S.In = In;
    for (BasicBlock *Pred : predecessors(B))
      if (Queued.insert(Pred).second)
        Work.push_back(Pred);
  }
  return true;
}

static bool win64SlotScalar(Type *Ty) {
  return (Ty->isPointerTy() && Ty->getPointerAddressSpace() == 0) ||
         (Ty->isIntegerTy() && Ty->getIntegerBitWidth() <= 64) ||
         Ty->isFloatTy() || Ty->isDoubleTy();
}

std::optional<unsigned> win64ScalarStackSlots(const FunctionType &Ty,
                                              const AttributeList &Attrs) {
  if (Ty.isVarArg() || Attrs.hasRetAttr(Attribute::InReg) ||
      (!Ty.getReturnType()->isVoidTy() && !win64SlotScalar(Ty.getReturnType())))
    return std::nullopt;
  constexpr Attribute::AttrKind Special[] = {
      Attribute::InReg,     Attribute::Nest,      Attribute::ByVal,
      Attribute::ByRef,     Attribute::InAlloca,  Attribute::Preallocated,
      Attribute::StructRet, Attribute::SwiftSelf, Attribute::SwiftError,
      Attribute::SwiftAsync};
  for (auto [I, Param] : enumerate(Ty.params())) {
    if (!win64SlotScalar(Param))
      return std::nullopt;
    const AttributeSet Set = Attrs.getParamAttrs(I);
    if (any_of(Special,
               [&](Attribute::AttrKind K) { return Set.hasAttribute(K); }))
      return std::nullopt;
  }
  const unsigned Count = Ty.getNumParams();
  return Count > 4 ? Count - 4 : 0;
}

bool isCallPressurePoint(const Instruction &I) {
  const auto *Call = dyn_cast<CallInst>(&I);
  return Call && !Call->isInlineAsm() && !isa<IntrinsicInst>(I) &&
         !Call->isMustTailCall() &&
         (Call->getCallingConv() == CallingConv::C ||
          Call->getCallingConv() == CallingConv::Fast ||
          Call->getCallingConv() == CallingConv::X86_64_SysV ||
          Call->getCallingConv() == CallingConv::X86_RegCall);
}

bool isWin64CallPressurePoint(const Instruction &I) {
  // On Windows x86-64, LLVM 23 uses csr_win64 for ordinary fastcc calls
  // as well as C and explicit win64cc calls.
  const auto *Call = dyn_cast<CallInst>(&I);
  return Call && !Call->isInlineAsm() && !isa<IntrinsicInst>(I) &&
         !Call->isMustTailCall() && !Call->hasOperandBundles() &&
         !Call->hasFnAttr(Attribute::ReturnsTwice) &&
         !Call->hasFnAttr("no_caller_saved_registers") &&
         !Call->hasFnAttr("no_callee_saved_registers") &&
         (Call->getCallingConv() == CallingConv::C ||
          Call->getCallingConv() == CallingConv::Fast ||
          Call->getCallingConv() == CallingConv::Win64);
}

static bool isAAPCS64CallPressurePoint(const Instruction &I) {
  const auto *Call = dyn_cast<CallInst>(&I);
  if (!Call || Call->isInlineAsm() || isa<IntrinsicInst>(I) ||
      Call->isMustTailCall() || Call->hasOperandBundles() ||
      Call->hasFnAttr(Attribute::ReturnsTwice) ||
      Call->hasFnAttr("no_caller_saved_registers") ||
      Call->hasFnAttr("no_callee_saved_registers") ||
      (Call->getCallingConv() != CallingConv::C &&
       Call->getCallingConv() != CallingConv::Fast))
    return false;
  // A scalable argument or return can select the SVE PCS implicitly.
  return !Call->getType()->isScalableTy() &&
         none_of(Call->args(),
                 [](const Use &U) { return U->getType()->isScalableTy(); });
}

static bool isExplicitCallPressurePoint(const Instruction &I,
                                        CallPressureABI ABI) {
  switch (ABI) {
  case CallPressureABI::SysV:
    return isCallPressurePoint(I);
  case CallPressureABI::Win64:
    return isWin64CallPressurePoint(I);
  case CallPressureABI::AAPCS64:
    return isAAPCS64CallPressurePoint(I);
  }
  llvm_unreachable("unknown call pressure ABI");
}

BlockPressure::BlockPressure(BasicBlock &B, const Liveness &Live,
                             bool IncludeScalarLibm,
                             ArrayRef<Instruction *> Trial,
                             CallPressureABI ABI) {
  if (Trial.empty())
    for (Instruction &I : B)
      Code.push_back(&I);
  else
    append_range(Code, Trial);
  for (auto [Point, I] : enumerate(Code))
    Positions[I] = Point;
  Before.resize(Code.size());
  Through.resize(Code.size());
  Calls.resize(Code.size());
  BitVector Set = Live.liveOut(B);
  Pressure P;
  for (unsigned ID : Set.set_bits())
    P += Live.weight(ID);
  const unsigned N = Set.size();

  for (unsigned Point = Code.size(); Point-- > 0;) {
    Instruction *I = Code[Point];
    const unsigned Def = Live.id(I);
    if (Def != N && Set.test(Def)) {
      if (!LastUse.count(I))
        LastUse[I] = Code.size();
      Set.reset(Def);
      P -= Live.weight(Def);
    }
    const bool ExplicitCall = isExplicitCallPressurePoint(*I, ABI);
    const bool MayLowerToLibm =
        IncludeScalarLibm &&
        (ABI == CallPressureABI::SysV || isa<IntrinsicInst>(I)) &&
        isScalarLibmCut(*I);
    if (ExplicitCall || MayLowerToLibm) {
      Calls.set(Point);
      Through[Point] = P;
      if (ABI == CallPressureABI::SysV)
        Pressured |= P.excess(Liveness::CallLimit) != 0;
    }
    if (!isa<PHINode>(I)) {
      for (const Use &U : I->operands()) {
        Value *V = U.get();
        const unsigned ID = Live.id(V);
        if (ID == N)
          continue;
        FirstUse[V] = Point;
        if (!Set.test(ID)) {
          LastUse[V] = Point;
          P += Live.weight(ID);
          Set.set(ID);
        } else if (!LastUse.count(V)) {
          LastUse[V] = Code.size();
        }
      }
    }
    Before[Point] = P;
    Pressured |= P.excess(Live.Limit) != 0;
  }
}

bool BlockPressure::isLiveBefore(const Value *V, unsigned Point) const {
  const auto Last = LastUse.find(V);
  return Last != LastUse.end() && Last->second >= Point;
}

void BlockPressure::summarize(Summary &S, Pressure Limit) const {
  for (unsigned I = 0; I < Code.size(); ++I) {
    S.Peak.GP = std::max(S.Peak.GP, Before[I].GP);
    S.Peak.SIMD = std::max(S.Peak.SIMD, Before[I].SIMD);
    S.ExcessPoints += Before[I].excess(Limit) != 0;
    S.Congestion += Before[I].congestion(Limit);
    if (Calls.test(I)) {
      ++S.CallSites;
      S.GPThroughCalls += Through[I].GP;
      S.GPCallExcess += std::max(0, Through[I].GP - Liveness::CallLimit.GP);
      S.SIMDThroughCalls += Through[I].SIMD;
      S.SIMDSpillThroughCalls += Through[I].SIMDSpill;
    }
  }
}

Summary summarize(Function &F, const Liveness &Live) {
  Summary S;
  for (BasicBlock &B : F)
    BlockPressure(B, Live).summarize(S, Live.Limit);
  return S;
}

} // namespace rssn
