#include "Planning.h"
#include "llvm/ADT/STLExtras.h"
#include "llvm/Analysis/TargetTransformInfo.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/Transforms/Utils/Local.h"

using namespace llvm;

namespace rssn {
namespace {

constexpr unsigned MaxPackWindow = 48;

struct Pack {
  SmallVector<LoadInst *, 8> Reads;
  SmallVector<Instruction *, 8> Chain;
};

bool scalarLadder(InsertElementInst &Final, Pack &P) {
  auto *Ty = cast<FixedVectorType>(Final.getType());
  P.Reads.resize(Ty->getNumElements());
  Value *Cursor = &Final;
  while (auto *Insert = dyn_cast<InsertElementInst>(Cursor)) {
    if (Insert != &Final && !Insert->hasOneUse())
      return false;
    auto *Lane = dyn_cast<ConstantInt>(Insert->getOperand(2));
    auto *Read = dyn_cast<LoadInst>(Insert->getOperand(1));
    if (!Lane || Lane->getValue().uge(P.Reads.size()) || !Read ||
        P.Reads[Lane->getZExtValue()])
      return false;
    P.Reads[Lane->getZExtValue()] = Read;
    P.Chain.push_back(Insert);
    Cursor = Insert->getOperand(0);
  }
  return isa<UndefValue>(Cursor) && all_of(P.Reads, [](auto *L) { return L; });
}

bool vectorConcat(ShuffleVectorInst &Final, Pack &P) {
  auto *Lo = dyn_cast<LoadInst>(Final.getOperand(0));
  auto *Hi = dyn_cast<LoadInst>(Final.getOperand(1));
  if (!Lo || !Hi || Lo == Hi)
    return false;
  auto *Half = dyn_cast<FixedVectorType>(Lo->getType());
  if (!Half || Final.getShuffleMask().size() != 2 * Half->getNumElements())
    return false;
  for (auto [Lane, Mask] : enumerate(Final.getShuffleMask()))
    if (Mask != static_cast<int>(Lane))
      return false;
  P.Reads = {Lo, Hi};
  P.Chain = {&Final};
  return true;
}

bool harmlessMetadata(const LoadInst &L) {
  SmallVector<std::pair<unsigned, MDNode *>, 4> MD;
  L.getAllMetadataOtherThanDebugLoc(MD);
  return all_of(MD, [](const auto &Entry) {
    switch (Entry.first) {
    case LLVMContext::MD_tbaa:
    case LLVMContext::MD_alias_scope:
    case LLVMContext::MD_noalias:
    case LLVMContext::MD_noundef:
      return true;
    default:
      return false;
    }
  });
}

bool contiguous(const Pack &P, Type *Ty) {
  const auto Start = accessAddress(*P.Reads.front());
  if (!Start)
    return false;
  const DataLayout &DL = P.Reads.front()->getDataLayout();
  const uint64_t Stride = DL.getTypeStoreSize(P.Reads.front()->getType());
  if (Stride * P.Reads.size() != DL.getTypeStoreSize(Ty))
    return false;

  for (auto [Lane, L] : enumerate(P.Reads)) {
    const auto Address = accessAddress(*L);
    int64_t Expected;
    if (AddOverflow(Start->second, static_cast<int64_t>(Lane * Stride),
                    Expected) ||
        !Address || Address->first != Start->first ||
        Address->second != Expected)
      return false;
  }
  return true;
}

bool fitsPressure(const Pack &Pack, BasicBlock &B, const Liveness &Live,
                  const Value *Pointer, unsigned Insert) {
  const BlockPressure P(B, Live);
  const auto LastPointer = P.LastUse.find(Pointer);
  if (LastPointer == P.LastUse.end() || LastPointer->second >= Insert)
    return true;
  const Pressure Added = Live.weight(Pointer);
  SmallVector<const Instruction *, 16> Removed(Pack.Reads.begin(),
                                                Pack.Reads.end());
  append_range(Removed, drop_begin(Pack.Chain));
  for (unsigned K = LastPointer->second; K <= Insert; ++K) {
    Pressure New = P.Before[K];
    Pressure NewCall = P.Through[K];
    if (K > LastPointer->second)
      New += Added;
    if (K < Insert)
      NewCall += Added;
    for (const Instruction *Old : Removed) {
      if (P.Positions.lookup(Old) >= K)
        continue;
      const unsigned End = P.LastUse.lookup(Old);
      if (End >= K)
        New -= Live.weight(Old);
      if (End > K)
        NewCall -= Live.weight(Old);
    }
    if (New.worsens(P.Before[K], Live.Limit) ||
        (P.Calls.test(K) && NewCall.worsens(P.Through[K], Liveness::CallLimit)))
      return false;
  }
  return true;
}

bool profitable(const Pack &P, Type *Ty, const TargetTransformInfo &TTI) {
  if (!TTI.isTypeLegal(Ty) || TTI.getRegUsageForType(Ty) != 1)
    return false;
  const Align Alignment = P.Reads.front()->getAlign();
  const auto WideRate = TTI.getMemoryOpCost(Instruction::Load, Ty, Alignment, 0);
  const auto WideSize = TTI.getMemoryOpCost(
      Instruction::Load, Ty, Alignment, 0, TargetTransformInfo::TCK_CodeSize);
  const auto WideLatency = TTI.getMemoryOpCost(
      Instruction::Load, Ty, Alignment, 0, TargetTransformInfo::TCK_Latency);
  InstructionCost NarrowRate, NarrowSize, AssemblyRate, AssemblyLatency;
  InstructionCost MaxLoadLatency;
  for (const LoadInst *L : P.Reads) {
    NarrowRate += TTI.getInstructionCost(L, TargetTransformInfo::TCK_RecipThroughput);
    NarrowSize += TTI.getInstructionCost(L, TargetTransformInfo::TCK_CodeSize);
    const auto Latency = TTI.getInstructionCost(L, TargetTransformInfo::TCK_Latency);
    if (!Latency.isValid())
      return false;
    MaxLoadLatency = std::max(MaxLoadLatency, Latency);
  }
  for (const Instruction *I : P.Chain) {
    AssemblyRate += TTI.getInstructionCost(I, TargetTransformInfo::TCK_RecipThroughput);
    AssemblyLatency += TTI.getInstructionCost(I, TargetTransformInfo::TCK_Latency);
  }
  const auto OldRate = NarrowRate + AssemblyRate;
  const auto OldPath = MaxLoadLatency + AssemblyLatency;
  return WideRate.isValid() && OldRate.isValid() && WideRate < OldRate &&
         WideSize.isValid() && NarrowSize.isValid() && WideSize <= NarrowSize &&
         WideLatency.isValid() && OldPath.isValid() && WideLatency <= OldPath;
}

LoadInst *buildWideLoad(ArrayRef<LoadInst *> Loads, Type *Ty,
                        Instruction &Anchor, const Twine &Name) {
  LoadInst *First = Loads.front();
  auto *Wide = new LoadInst(Ty, First->getPointerOperand(), Name, false,
                            First->getAlign(), Anchor.getIterator());
  Wide->setDebugLoc(Anchor.getDebugLoc());
  if (all_of(Loads, [](const LoadInst *L) {
        return L->hasMetadata(LLVMContext::MD_noundef);
      }))
    Wide->setMetadata(LLVMContext::MD_noundef,
                      First->getMetadata(LLVMContext::MD_noundef));
  return Wide;
}

void eraseReads(ArrayRef<LoadInst *> Loads, BasicBlock &B) {
  SmallPtrSet<Instruction *, 8> LocalAddresses;
  for (LoadInst *L : Loads) {
    if (auto *Address = dyn_cast<GetElementPtrInst>(L->getPointerOperand()))
      if (Address->getParent() == &B)
        LocalAddresses.insert(Address);
    salvageDebugInfo(*L);
    L->eraseFromParent();
  }
  // Leave cross-block address cleanup to the ordinary dead-code passes.
  for (Instruction *Address : LocalAddresses)
    if (Address->use_empty()) {
      salvageDebugInfo(*Address);
      Address->eraseFromParent();
    }
}

struct Lane {
  LoadInst *Read;
  BinaryOperator *Use;
  unsigned DefPoint;
  unsigned UsePoint;
};

bool laneOpcode(const Instruction &I) {
  return I.getOpcode() == Instruction::FAdd ||
         I.getOpcode() == Instruction::FSub ||
         I.getOpcode() == Instruction::FMul;
}

struct PackedOperation {
  BinaryOperator *First;
  Value *Common;
  unsigned Side;
};

std::optional<PackedOperation> compatibleOperation(ArrayRef<Lane> Lanes,
                                                    const BlockPressure &P) {
  const Lane &First = *min_element(Lanes, [](const Lane &A, const Lane &B) {
    return A.UsePoint < B.UsePoint;
  });
  BinaryOperator *Prototype = First.Use;
  const unsigned Side = Prototype->getOperand(0) == First.Read ? 0 : 1;
  Value *Common = Prototype->getOperand(1 - Side);
  unsigned Last = First.UsePoint;
  for (const Lane &L : Lanes) {
    Last = std::max(Last, L.UsePoint);
    if (L.Use->getOpcode() != Prototype->getOpcode() ||
        L.Use->getFastMathFlags() != Prototype->getFastMathFlags() ||
        L.Use->hasMetadata(LLVMContext::MD_fpmath) ||
        L.Use->getOperand(Side) != L.Read ||
        L.Use->getOperand(1 - Side) != Common)
      return std::nullopt;
  }
  if (Last - First.UsePoint > 8)
    return std::nullopt;
  for (unsigned K = First.UsePoint; K < Last; ++K)
    if (isa<CallBase>(P.Code[K]))
      return std::nullopt;
  return PackedOperation{Prototype, Common, Side};
}

InstructionCost packedOperationCost(const PackedOperation &Op, FixedVectorType *Ty,
                                    const TargetTransformInfo &TTI,
                                    TargetTransformInfo::TargetCostKind Kind) {
  InstructionCost Cost = TTI.getArithmeticInstrCost(Op.First->getOpcode(), Ty, Kind);
  if (!isa<Constant>(Op.Common)) {
    Cost += TTI.getVectorInstrCost(Instruction::InsertElement, Ty, Kind, 0);
    Cost += TTI.getShuffleCost(TargetTransformInfo::SK_Broadcast, Ty, Ty, {}, Kind);
  }
  return Cost;
}

// Score carrier and projections together. At an original user, the scalar
// projection briefly coexists with the carrier if some later lane still needs
// it. At a call, independent scalar spill pairs are traded for one vector pair.
int laneGain(ArrayRef<Lane> Lanes, FixedVectorType *Ty, unsigned Anchor,
             const PackedOperation &Op, const BlockPressure &P, const Liveness &Live,
             const TargetTransformInfo &TTI) {
  const Pressure Carrier = Live.typeWeight(Ty);
  if (Carrier.SIMD != 1 || !TTI.isTypeLegal(Ty))
    return 0;
  const unsigned Start = min_element(Lanes, [](const Lane &A, const Lane &B) {
                           return A.DefPoint < B.DefPoint;
                         })->DefPoint;
  const unsigned End = max_element(Lanes, [](const Lane &A, const Lane &B) {
                         return A.UsePoint < B.UsePoint;
                       })->UsePoint;
  const unsigned VectorUse = P.Positions.lookup(Op.First);
  const Value *Pointer = Lanes.front().Read->getPointerOperand();
  const unsigned PointerEnd = P.LastUse.lookup(Pointer);
  const Pressure PW = Live.weight(Pointer);
  int Gain = 0;
  for (unsigned K = Start + 1; K <= End; ++K) {
    Pressure New = P.Before[K];
    Pressure Through = P.Through[K];
    if (K <= Anchor) {
      if (K > PointerEnd)
        New += PW;
      if (K >= PointerEnd && K < Anchor)
        Through += PW;
    } else {
      New += Carrier;
      if (K < End)
        Through += Carrier;
    }

    Pressure Projections;
    for (const Lane &L : Lanes) {
      if (L.DefPoint >= K)
        continue;
      const Pressure W = Live.weight(L.Read);
      if (L.UsePoint >= K)
        New -= W;
      if (L.UsePoint > K)
        Through -= W;
      if (L.UsePoint == K)
        Projections += W;
    }
    if (New.worsens(P.Before[K], Live.Limit) ||
        (P.Calls.test(K) && Through.worsens(P.Through[K], Liveness::CallLimit)))
      return 0;
    Gain += P.Before[K].congestion(Live.Limit) - New.congestion(Live.Limit);
    if (P.Calls.test(K))
      Gain += CallWeight * (P.Through[K].callCost(Liveness::CallLimit) -
                           Through.callCost(Liveness::CallLimit));

    if (K == VectorUse && !isa<Constant>(Op.Common)) {
      // Account for the transient broadcast register before the vector op.
      Pressure WithSplat = New;
      WithSplat += Carrier;
      if (WithSplat.worsens(P.Before[K], Live.Limit))
        return 0;
    }
    if (K == End)
      New -= Carrier;
    New += Projections;
    if (New.worsens(P.Before[K], Live.Limit))
      return 0;
  }

  const Align Alignment = Lanes.front().Read->getAlign();
  auto NewRate = TTI.getMemoryOpCost(Instruction::Load, Ty, Alignment, 0);
  auto NewSize = TTI.getMemoryOpCost(Instruction::Load, Ty, Alignment, 0,
                                    TargetTransformInfo::TCK_CodeSize);
  NewRate += packedOperationCost(Op, Ty, TTI, TargetTransformInfo::TCK_RecipThroughput);
  NewSize += packedOperationCost(Op, Ty, TTI, TargetTransformInfo::TCK_CodeSize);
  const auto PackedLatency =
      packedOperationCost(Op, Ty, TTI, TargetTransformInfo::TCK_Latency);
  const auto ScalarLatency =
      TTI.getInstructionCost(Op.First, TargetTransformInfo::TCK_Latency);
  if (!PackedLatency.isValid() || !ScalarLatency.isValid())
    return 0;
  const int AddedLatency =
      std::clamp<int64_t>((PackedLatency - ScalarLatency).getValue(), 0, 16);
  InstructionCost OldRate, OldSize;
  int Scheduling = 2 * AddedLatency;
  for (auto [Index, L] : enumerate(Lanes)) {
    const unsigned LaneNo = Index;
    OldRate += TTI.getInstructionCost(L.Read, TargetTransformInfo::TCK_RecipThroughput) +
               TTI.getInstructionCost(L.Use, TargetTransformInfo::TCK_RecipThroughput);
    OldSize += TTI.getInstructionCost(L.Read, TargetTransformInfo::TCK_CodeSize) +
               TTI.getInstructionCost(L.Use, TargetTransformInfo::TCK_CodeSize);
    SmallVector<std::tuple<Value *, User *, int>, 4> Users;
    for (User *User : L.Use->users())
      Users.emplace_back(L.Use, User, LaneNo);
    NewRate += TTI.getVectorInstrCost(Instruction::ExtractElement, Ty,
                                      TargetTransformInfo::TCK_RecipThroughput,
                                      LaneNo, L.Use, Users);
    NewSize += TTI.getVectorInstrCost(Instruction::ExtractElement, Ty,
                                      TargetTransformInfo::TCK_CodeSize, LaneNo);
    const auto Latency = TTI.getVectorInstrCost(Instruction::ExtractElement, Ty,
                                                TargetTransformInfo::TCK_Latency,
                                                LaneNo);
    if (!Latency.isValid())
      return 0;
    Scheduling += 2 * std::clamp<int64_t>(Latency.getValue(), 0, 16);
    Scheduling += 3; // The old scalar operation could fold its stack reload.
  }
  if (!NewRate.isValid() || !NewSize.isValid() || !OldRate.isValid() ||
      !OldSize.isValid() || NewRate > OldRate + 2 || NewSize > OldSize + 4)
    return 0;
  const auto RateDelta = NewRate - OldRate;
  const auto SizeDelta = NewSize - OldSize;
  return Gain - Scheduling - 8 * RateDelta.getValue() - 2 * SizeDelta.getValue();
}

} // namespace

bool tryFusePackedLoads(BasicBlock &B, Liveness &Live,
                        const TargetTransformInfo &TTI,
                        MemoryModel &Memory, Changes &Result) {
  SmallVector<Instruction *, 64> Code;
  DenseMap<const Instruction *, unsigned> Positions;
  for (Instruction &I : B) {
    Positions[&I] = Code.size();
    Code.push_back(&I);
  }
  for (auto [Index, Final] : enumerate(Code)) {
    const unsigned Point = Index;
    Type *Ty = Final->getType();
    auto *VT = dyn_cast<FixedVectorType>(Ty);
    if (!VT || Final->use_empty() || Ty->getPrimitiveSizeInBits() < 128 ||
        Ty->getPrimitiveSizeInBits() > 256 || VT->getNumElements() > 8)
      continue;
    const Type *Element = VT->getElementType();
    if (!Element->isFloatTy() && !Element->isDoubleTy() &&
        !Element->isIntegerTy(32) && !Element->isIntegerTy(64))
      continue;
    // A vector that is immediately scalarized offers no packed consumer.
    if (all_of(Final->users(), [](const User *U) {
          return isa<ExtractElementInst>(U);
        }))
      continue;

    Pack P;
    if (auto *Insert = dyn_cast<InsertElementInst>(Final)) {
      if (!scalarLadder(*Insert, P))
        continue;
    } else if (auto *Shuffle = dyn_cast<ShuffleVectorInst>(Final)) {
      if (!vectorConcat(*Shuffle, P))
        continue;
    } else {
      continue;
    }

    if (any_of(P.Chain, [&](const Instruction *I) { return I->getParent() != &B; }) ||
        any_of(P.Reads, [&](const LoadInst *L) {
          return L->getParent() != &B || !L->hasOneUse() || !L->isSimple() ||
                 L->getPointerAddressSpace() != 0 || !harmlessMetadata(*L) ||
                 Positions.lookup(L) >= Point ||
                 Point - Positions.lookup(L) > MaxPackWindow;
        }) ||
        !contiguous(P, Ty) || !profitable(P, Ty, TTI) ||
        !fitsPressure(P, B, Live, P.Reads.front()->getPointerOperand(), Point))
      continue;

    if (any_of(P.Reads, [&](const LoadInst *L) {
          return !Memory.unchanged(*L, Code, Positions.lookup(L) + 1, Point);
        }))
      continue;

    // Put the wide load at the old packed value's definition. It has exactly
    // the same result live range, and all old lanes and partial packs die.
    auto *Wide = buildWideLoad(P.Reads, Ty, *Final, Final->getName() + ".rssn.wide");
    // Successors see a different SSA identity with the same type and
    // edge-level lifetime until liveness is rebuilt after this prepass.
    Live.rebind(*Final, *Wide);
    Final->replaceAllUsesWith(Wide);
    for (Instruction *I : P.Chain) {
      salvageDebugInfo(*I);
      I->eraseFromParent();
    }
    eraseReads(P.Reads, B);
    ++Result.PacksFused;
    Result.ReadsFused += P.Reads.size();
    return true;
  }
  return false;
}

bool tryFuseLaneLoads(BasicBlock &B, Liveness &Live,
                     const TargetTransformInfo &TTI, MemoryModel &Memory,
                     Changes &Result) {
  if (!any_of(B, isCallPressurePoint))
    return false;
  BlockPressure P(B, Live);
  SmallVector<Lane, 4> Best;
  FixedVectorType *BestType = nullptr;
  unsigned BestAnchor = 0;
  std::optional<PackedOperation> BestOp;
  int BestGain = 7;
  for (unsigned Call : P.Calls.set_bits()) {
    if (P.Through[Call].SIMD < 2)
      continue;
    using LoadKey = std::pair<AddressKey, Type *>;
    SmallVector<std::pair<LoadKey, Lane>, 16> Available;
    DenseMap<LoadKey, unsigned> ByAddress;
    const unsigned Start = Call > MaxPackWindow ? Call - MaxPackWindow : 0;
    for (unsigned K = Start; K < Call; ++K) {
      auto *L = dyn_cast<LoadInst>(P.Code[K]);
      if (!L || !L->hasOneUse() || !simpleLoad(*L, Live) || !harmlessMetadata(*L))
        continue;
      Type *Scalar = L->getType();
      if (!Scalar->isDoubleTy() && !Scalar->isFloatTy())
        continue;
      auto *Use = dyn_cast<BinaryOperator>(*L->user_begin());
      if (!Use || !laneOpcode(*Use) || Use->getParent() != &B ||
          Use->hasNUsesOrMore(9))
        continue;
      const unsigned End = P.Positions.lookup(Use);
      if (End <= Call || End - Call > MaxPackWindow)
        continue;
      if (const auto Address = accessAddress(*L)) {
        const LoadKey Key{*Address, Scalar};
        if (ByAddress.try_emplace(Key, Available.size()).second)
          Available.push_back({Key, Lane{L, Use, K, End}});
      }
    }

    for (const auto &[Key, First] : Available) {
      Type *Scalar = First.Read->getType();
      const unsigned Count = Scalar->isDoubleTy() ? 2 : 4;
      const unsigned Stride = 16 / Count;
      SmallVector<Lane, 4> Lanes{First};
      unsigned Anchor = First.DefPoint;
      for (unsigned K = 1; K < Count; ++K) {
        int64_t Offset;
        if (AddOverflow(Key.first.second, static_cast<int64_t>(K * Stride), Offset))
          break;
        const auto It = ByAddress.find({{Key.first.first, Offset}, Scalar});
        if (It == ByAddress.end())
          break;
        const Lane &L = Available[It->second].second;
        Anchor = std::max(Anchor, L.DefPoint);
        Lanes.push_back(L);
      }
      if (Lanes.size() != Count)
        continue;
      auto Op = compatibleOperation(Lanes, P);
      if (!Op)
        continue;

      unsigned Blocked = 0;
      for (const Lane &L : Lanes) {
        if (!Memory.unchanged(*L.Read, P.Code, L.DefPoint + 1, Anchor)) {
          Blocked = 0;
          break;
        }
        // Prefer moving original loads when they can cross the call. Keep
        // the packed carrier for lanes whose original reads must stay early.
        Blocked += !Memory.unchanged(*L.Read, P.Code, Anchor + 1, Call + 1);
      }
      if (Blocked < 2)
        continue;
      auto *Vector = FixedVectorType::get(Scalar, Count);
      const int Gain = laneGain(Lanes, Vector, Anchor, *Op, P, Live, TTI);
      if (Gain > BestGain) {
        BestGain = Gain;
        Best = Lanes;
        BestType = Vector;
        BestAnchor = Anchor;
        BestOp = Op;
      }
    }
  }

  if (Best.empty())
    return false;
  SmallVector<LoadInst *, 4> Reads;
  for (const Lane &L : Best)
    Reads.push_back(L.Read);
  auto *Wide = buildWideLoad(Reads, BestType, *P.Code[BestAnchor], "rssn.carrier");
  IRBuilder<> Packed(BestOp->First);
  Packed.setFastMathFlags(BestOp->First->getFastMathFlags());
  Value *Common = Packed.CreateVectorSplat(Best.size(), BestOp->Common, "rssn.splat");
  Value *Left = BestOp->Side ? Common : Wide;
  Value *Right = BestOp->Side ? Wide : Common;
  Value *Vector = Packed.CreateBinOp(
      static_cast<Instruction::BinaryOps>(BestOp->First->getOpcode()),
      Left, Right, "rssn.packed");

  for (auto [I, L] : enumerate(Best)) {
    IRBuilder<> Builder(L.Use);
    Value *Projection = Builder.CreateExtractElement(
        Vector, Builder.getInt32(I), L.Use->getName() + ".rssn.lane");
    Live.rebind(*L.Use, *Projection);
    L.Use->replaceAllUsesWith(Projection);
    L.Use->eraseFromParent();
  }
  eraseReads(Reads, B);
  ++Result.LaneGroupsFused;
  Result.ReadsFused += Reads.size();
  return true;
}

} // namespace rssn
