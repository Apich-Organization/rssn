#include "MemoryModel.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/Analysis/AliasAnalysis.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/IR/IntrinsicInst.h"

using namespace llvm;

namespace rssn {

MemoryModel::FrameObject MemoryModel::inspect(const AllocaInst &Root) {
  FrameObject Object;
  if (!Root.isStaticAlloca() || Root.getAddressSpace() != 0 ||
      Root.isUsedWithInAlloca() || Root.isSwiftError() ||
      Root.getFunction()->hasGC() ||
      Root.getFunction()->hasFnAttribute(Attribute::PresplitCoroutine))
    return Object;

  // A syntactic ownership proof. Every derived pointer is followed, and
  // every call that receives provenance is retained even if noncapturing.
  SmallVector<const Value *, 16> Work{&Root};
  SmallPtrSet<const Value *, 16> Seen;
  unsigned Uses = 0;
  while (!Work.empty()) {
    const Value *Ptr = Work.pop_back_val();
    if (!Seen.insert(Ptr).second)
      continue;
    for (const Use &U : Ptr->uses()) {
      if (++Uses > 96)
        return {};
      const User *User = U.getUser();
      if (const auto *GEP = dyn_cast<GetElementPtrInst>(User)) {
        if (U.getOperandNo() != GEP->getPointerOperandIndex())
          return {};
        Work.push_back(GEP);
        continue;
      }
      if (const auto *Load = dyn_cast<LoadInst>(User)) {
        if (Load->getPointerOperand() != Ptr)
          return {};
        continue;
      }
      if (const auto *Store = dyn_cast<StoreInst>(User)) {
        if (U.getOperandNo() != Store->getPointerOperandIndex())
          return {};
        continue;
      }
      if (const auto *CB = dyn_cast<CallBase>(User)) {
        if (CB->isInlineAsm() || CB->hasOperandBundles() ||
            !CB->isArgOperand(&U))
          return {};
        if (const auto *II = dyn_cast<IntrinsicInst>(CB)) {
          if (II->isLifetimeStartOrEnd())
            continue;
          // localescape and other implicit frame access mechanisms must
          // never be treated as an ordinary noncapturing call.
          if (!isa<MemIntrinsic>(II))
            return {};
        }
        if (!CB->doesNotCapture(CB->getArgOperandNo(&U)))
          return {};
        Object.ReceivesAddress.insert(CB);
        continue;
      }
      // Includes storing the pointer, ptrtoint/ptrtoaddr, integer provenance
      // recovery, pointer comparison, ambiguous PHIs and escaping returns.
      return {};
    }
  }
  Object.Private = true;
  return Object;
}

bool MemoryModel::inaccessibleAt(const LoadInst &L, const CallBase &Call) {
  const auto *Root =
      dyn_cast<AllocaInst>(getUnderlyingObject(L.getPointerOperand()));
  if (!Root)
    return false;
  auto [It, New] = Frames.try_emplace(Root);
  if (New)
    It->second = inspect(*Root);
  if (!It->second.Private || It->second.ReceivesAddress.contains(&Call))
    return false;
  ProvenBoundaries.insert({Root, &Call});
  return true;
}

static bool volatileAccess(const Instruction &I) {
  if (auto *Load = dyn_cast<LoadInst>(&I))
    return Load->isVolatile();
  if (auto *Store = dyn_cast<StoreInst>(&I))
    return Store->isVolatile();
  return false;
}

bool MemoryModel::unchanged(const LoadInst &L, ArrayRef<Instruction *> Code,
                            unsigned From, unsigned To) {
  const MemoryLocation Loc = MemoryLocation::get(&L);
  for (unsigned K = From; K < To; ++K) {
    const Instruction *I = Code[K];
    if (I->isAtomic() || volatileAccess(*I))
      return false;
    if (auto *Call = dyn_cast<CallBase>(I)) {
      if (Call->isInlineAsm() || Call->isConvergent() ||
          Call->hasOperandBundles() || Call->hasFnAttr(Attribute::ReturnsTwice))
        return false;
      if (isa<IntrinsicInst>(Call)) {
        // In particular do not bypass stackrestore, lifetime or coroutine
        // intrinsics with the frame proof.
        if (!Call->doesNotAccessMemory() ||
            !Call->hasFnAttr(Attribute::Speculatable))
          return false;
      } else {
        const bool Regular = Call->doesNotThrow() &&
                             Call->hasFnAttr(Attribute::NoFree) &&
                             Call->hasFnAttr(Attribute::NoSync);
        if ((!Regular ||
             (I->mayWriteToMemory() && isModSet(AA.getModRefInfo(I, Loc)))) &&
            !inaccessibleAt(L, *Call))
          return false;
        continue;
      }
    }
    if (I->mayWriteToMemory() && isModSet(AA.getModRefInfo(I, Loc)))
      return false;
  }
  return true;
}

} // namespace rssn
