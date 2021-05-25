#include "SimTracker/TrackTriggerAssociation/interface/StubAssociation.h"

#include <map>
#include <vector>
#include <utility>
#include <numeric>

using namespace std;

namespace tt {

  // insert a TPPtr and its associated collection of TTstubRefs into the underlayering maps
  void StubAssociation::insert(const TPPtr& tpPtr, const vector<TTStubRef>& ttSTubRefs) {
    mapTPPtrsTTStubRefs_.insert({tpPtr, ttSTubRefs});
    for (const TTStubRef& ttSTubRef : ttSTubRefs)
      mapTTStubRefsTPPtrs_[ttSTubRef].push_back(tpPtr);
  }

  // returns collection of TPPtrs associated to given TTstubRef
  const vector<TPPtr>& StubAssociation::findTrackingParticlePtrs(const TTStubRef& ttStubRef) const {
    const auto it = mapTTStubRefsTPPtrs_.find(ttStubRef);
    return it != mapTTStubRefsTPPtrs_.end() ? it->second : emptyTPPtrs_;
  }

  // returns collection of TTStubRefs associated to given TPPtr
  const vector<TTStubRef>& StubAssociation::findTTStubRefs(const TPPtr& tpPtr) const {
    const auto it = mapTPPtrsTTStubRefs_.find(tpPtr);
    return it != mapTPPtrsTTStubRefs_.end() ? it->second : emptyTTStubRefs_;
  }

  // Get all TPs that are matched to these stubs in at least 'tpMinLayers' layers 
  vector<TPPtr> StubAssociation::associate(const vector<TTStubRef>& ttStubRefs) const {
    // count associated layer for each TP
    map<TPPtr, set<int>> m;
    for (const TTStubRef& ttStubRef : ttStubRefs)
      for (const TPPtr& tpPtr : findTrackingParticlePtrs(ttStubRef))
        m[tpPtr].insert(setup_->layerId(ttStubRef));
    // count matched TPs
    auto sum = [this](int& sum, const pair<TPPtr, set<int>>& p) {
      return sum += ((int)p.second.size() < setup_->tpMinLayers() ? 0 : 1);
    };
    const int nTPs = accumulate(m.begin(), m.end(), 0, sum);
    vector<TPPtr> tpPtrs;
    tpPtrs.reserve(nTPs);
    // fill and return matched TPs
    for (const auto& p : m)
      if ((int)p.second.size() >= setup_->tpMinLayers())
        tpPtrs.push_back(p.first);
    return tpPtrs;
  }

} // namespace tt