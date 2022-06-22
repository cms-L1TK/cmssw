#include "L1Trigger/TrackerTFP/interface/MiniHoughTransform.h"

#include <numeric>
#include <algorithm>
#include <iterator>
#include <deque>
#include <vector>
#include <set>
#include <utility>
#include <cmath>

using namespace std;
using namespace edm;
using namespace tt;

namespace trackerTFP {

  MiniHoughTransform::MiniHoughTransform(const ParameterSet& iConfig,
                                         const Setup* setup,
                                         const DataFormats* dataFormats,
                                         int region)
      : enableTruncation_(iConfig.getParameter<bool>("EnableTruncation")),
        setup_(setup),
        dataFormats_(dataFormats),
        inv2R_(dataFormats_->format(Variable::inv2R, Process::ht)),
        phiT_(dataFormats_->format(Variable::phiT, Process::ht)),
        region_(region),
        numBinsInv2R_(setup_->htNumBinsInv2R()),
        numStages_(setup->mhtNumStages()),
        numCells_(setup_->mhtNumCells()),
        input_(numBinsInv2R_) {}

  // read in and organize input product (fill vector input_)
  void MiniHoughTransform::consume(const StreamsStub& streams) {
    auto valid = [](int sum, const FrameStub& frame) { return sum + (frame.first.isNonnull() ? 1 : 0); };
    int nStubsHT(0);
    for (int binInv2R = 0; binInv2R < numBinsInv2R_; binInv2R++) {
      const StreamStub& stream = streams[region_ * numBinsInv2R_ + binInv2R];
      nStubsHT += accumulate(stream.begin(), stream.end(), 0, valid);
    }
    stubs_.reserve(nStubsHT);
    states_.reserve(nStubsHT);
    for (int binInv2R = 0; binInv2R < numBinsInv2R_; binInv2R++) {
      const int inv2R = inv2R_.toSigned(binInv2R);
      const StreamStub& stream = streams[region_ * numBinsInv2R_ + binInv2R];
      vector<State*>& states = input_[binInv2R];
      states.reserve(stream.size());
      // Store input stubs in vector, so rest of MHT algo can work with pointers to them (saves CPU)
      for (const FrameStub& frame : stream) {
        State* state = nullptr;
        if (frame.first.isNonnull()) {
          stubs_.emplace_back(frame, dataFormats_, inv2R);
          states_.emplace_back(stubs_.back(), setup_);
          state = &states_.back();
        }
        states.push_back(state);
      }
    }
  }

  // fill output products
  void MiniHoughTransform::produce(StreamsStub& accepted, StreamsStub& lost) {
    // fill MHT cells
    for (int channel = 0; channel < numBinsInv2R_; channel++) {
      vector<State*>& stream = input_[channel];
      /*for (State* stub : stream) {
        if (!stub)
          continue;
        cout << setw(12)<< stub->r_ << " "
             << setw(12)<< stub->phi_ << " "
             << setw(12)<< stub->dPhi_ << " "
             << stub->trackId_ << " "
             << endl;
      }
      if (!stream.empty()) {
        cout << dataFormats_->base(Variable::inv2R, Process::mht) * 2. << " " << dataFormats_->base(Variable::phiT, Process::mht) * 2. << endl;
        cout << dataFormats_->base(Variable::inv2R, Process::ht) << " " << dataFormats_->base(Variable::phiT, Process::ht) << endl;
        cout << endl;
      }*/
      for (int iter = 0; iter < numStages_; iter++)
        stage(iter, stream);
      // fill output productd
      StreamStub& out = accepted[region_ * numBinsInv2R_ + channel];
      out.reserve(stream.size());
      for (State* s : stream) {
        if (!s) {
          out.emplace_back(FrameStub());
          continue;
        }
        StubMHT stub(*s->stubHT_, s->phi_, s->inv2R_, s->phiT_);
        out.emplace_back(stub.frame());
      }
    }
  }

  void MiniHoughTransform::stage(int iter, vector<State*>& stream) {
    const double bp = pow(2., -iter);
    const double basePhiT = phiT_.base() * bp;
    const double baseInv2R = inv2R_.base() * bp;
    int id;
    auto different = [&id](State* stub) { return stub && (id != stub->trackId_); };
    auto moreLayer = [](const TTBV& lhs, const TTBV& rhs) { return lhs.count() < rhs.count(); };
    for (auto it = stream.begin(); it != stream.end();) {
      auto start = it;
      id = (*it) ? (*it)->trackId_ : -1;
      auto end = find_if(it, stream.end(), different);
      it = end;
      if (id == -1)
        continue;
      const int size = accumulate(start, end, 0, [](int& sum, const auto& stub){ return sum += (stub ? 1 : 0); });
      // create finer track candidates stub container
      vector<vector<State*>> mhtCells(numCells_);
      for (vector<State*>& mhtCell : mhtCells)
        mhtCell.reserve(size);
      // first loop over track: fill finer track candidates stub container
      for (auto stub = start; stub != end; stub++) {
        if (!*stub)
          continue;
        const double r = (*stub)->r_;
        const double chi = (*stub)->phi_;
        const double dChi = (*stub)->dPhi_;
        // identify finer track candidates for this stub
        // 0 and 1 belong to the MHT cells with larger inv2R; 0 and 2 belong to those with smaller track PhiT
        TTBV cells(0, numCells_);
        const bool compA = 2. * abs(chi) < basePhiT + abs(r * baseInv2R) + dChi;
        const bool compB = 2. * abs(chi) < basePhiT + dChi;
        const bool compC = 2. * abs(chi) < abs(r * baseInv2R) + dChi;
        const bool compD = 2. * abs(chi) < dChi;
        if (chi >= 0. && r >= 0.) { if (compA) cells.set(3); if (compB) cells.set(1); if (compC) cells.set(2); if (compD) cells.set(0); }
        if (chi >= 0. && r <  0.) { if (compA) cells.set(1); if (compB) cells.set(3); if (compC) cells.set(0); if (compD) cells.set(2); }
        if (chi <  0. && r >= 0.) { if (compA) cells.set(0); if (compB) cells.set(2); if (compC) cells.set(1); if (compD) cells.set(3); }
        if (chi <  0. && r <  0.) { if (compA) cells.set(2); if (compB) cells.set(0); if (compC) cells.set(3); if (compD) cells.set(1); }
        /*cout << setw(12) << r << " "
             << setw(12) << chi << " "
             << setw(12) << dChi << " "
             << "| "
             << setw(12) << baseInv2R << " "
             << setw(12) << basePhiT << " "
             << "| ";
        for (int cell : cells.ids())
          cout << cell << " ";
        cout << "| "
             << 2. * abs(chi) << " " << basePhiT + dChi << " " << abs(r * baseInv2R) + dChi << " " << dChi << " ";
        cout << endl;*/
        // organise stubs in finer track candidates
        for (int cell : cells.ids())
          mhtCells[cell].push_back(*stub);
      }
      // perform finer pattern recognition and keep only cells with max number of layer
      vector<TTBV> hitPatterns(numCells_, TTBV(0, setup_->numLayers()));
      for (int cell = 0; cell < numCells_; cell++)
        for (State* stub : mhtCells[cell])
          hitPatterns[cell].set(stub->layer_);
      const int maxLayer = max_element(hitPatterns.begin(), hitPatterns.end(), moreLayer)->count();
      const int limit = max(maxLayer, setup_->mhtMinLayers());
      TTBV cells(0, numCells_);
      for (int cell = 0; cell < numCells_; cell++)
        if (hitPatterns[cell].count() == limit)
          cells.set(cell);
      // identify merged finer track
      set<State*> track;
      pair<int, int> inv2Rs(setup_->mhtNumBinsInv2R(), -1);
      pair<int, int> phiTs(setup_->mhtNumBinsPhiT(), -1);
      for (int cell : cells.ids()) {
        const int inv2R = cell / setup_->mhtNumBinsPhiT();
        const int phiT = cell % setup_->mhtNumBinsPhiT();
        inv2Rs = {min(inv2Rs.first, inv2R), max(inv2Rs.second, inv2R)};
        phiTs = {min(phiTs.first, phiT), max(phiTs.second, phiT)};
        for (State* stub : mhtCells[cell])
          track.insert(stub);
      }
      const int p = pow(2, numStages_ - iter);
      const int inv2R = (inv2Rs.first + inv2Rs.second - setup_->mhtNumBinsInv2R() + 1) * p / 2;
      const int phiT = (phiTs.first + phiTs.second - setup_->mhtNumBinsPhiT() + 1) * p / 2;
      /*if (cells.any()) {
        cout << setw(2) << inv2R << " " << setw(2) << phiT << " | ";
        for (int cell : cells.ids())
          cout << cell << " ";
        cout << endl;
      }*/
      // second loop over track, kill bad stubs and update track parameter and stub residual
      for (auto& stub = start; stub != end; stub++) {
        if (!*stub)
          continue;
        const auto s = find(track.begin(), track.end(), *stub);
        if (s == track.end())
          *stub = nullptr;
        else
          (*stub)->update(inv2R, phiT, dataFormats_);
      }
    }
    stream.erase(remove(stream.begin(), stream.end(), nullptr), stream.end());
    /*if (!stream.empty())
      cout << endl;*/
  }

}  // namespace trackerTFP
