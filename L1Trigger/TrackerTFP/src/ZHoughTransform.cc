#include "L1Trigger/TrackerTFP/interface/ZHoughTransform.h"

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

  ZHoughTransform::ZHoughTransform(const ParameterSet& iConfig,
                                   const Setup* setup,
                                   const DataFormats* dataFormats,
                                   int region)
      : enableTruncation_(iConfig.getParameter<bool>("EnableTruncation")),
        setup_(setup),
        dataFormats_(dataFormats),
        region_(region),
        input_(dataFormats->numChannel(Process::mht))  {}

  // read in and organize input product (fill vector input_)
  void ZHoughTransform::consume(const StreamsStub& streams) {
    auto valid = [](int sum, const FrameStub& frame) { return sum + (frame.first.isNonnull() ? 1 : 0); };
    const int offset = region_ * dataFormats_->numChannel(Process::mht);
    int nStubsMHT(0);
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      const StreamStub& stream = streams[offset + channel];
      nStubsMHT += accumulate(stream.begin(), stream.end(), 0, valid);
    }
    stubs_.reserve(nStubsMHT);
    states_.reserve(nStubsMHT);
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      const StreamStub& stream = streams[offset + channel];
      vector<State*>& states = input_[channel];
      states.reserve(stream.size());
      // Store input stubs in vector, so rest of ZHT algo can work with pointers to them (saves CPU)
      for (const FrameStub& frame : stream) {
        State* state = nullptr;
        if (frame.first.isNonnull()) {
          StubMHT stubMHT(frame, dataFormats_);
          stubs_.emplace_back(frame, dataFormats_);
          states_.emplace_back(stubs_.back(), setup_);
          state = &states_.back();
        }
        states.push_back(state);
      }
    }
  }

  // fill output products
  void ZHoughTransform::produce(StreamsStub& accepted, StreamsStub& lost) {
    // fill MHT cells
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      vector<State*>& stream = input_[channel];
      /*for (State* stub : stream) {
        if (!stub)
          continue;
        cout << setw(12)<< stub->r_ << " "
             << setw(12)<< stub->z_ << " "
             << setw(12)<< stub->dZ_ << " "
             << stub->trackId_ << " "
             << endl;
      }
      if (!stream.empty()) {
        cout << dataFormats_->base(Variable::cot, Process::zht) * 2. << " " << dataFormats_->base(Variable::zT, Process::zht) * 2. << endl;
        cout << endl;
      }*/
      for (int iter = 0; iter < setup_->zhtNumStages(); iter++)
        stage(iter, stream);
      // fill output productd
      StreamStub& out = accepted[region_ * dataFormats_->numChannel(Process::mht) + channel];
      out.reserve(stream.size());
      for (State* state : stream) {
        if (!state) {
          out.emplace_back(FrameStub());
          continue;
        }
        StubZHT stub(*state->stub_, state->z_, state->cot_, state->zT_);
        out.emplace_back(stub.frame());
      }
    }
  }

  // perform finer pattern recognition per track
  void ZHoughTransform::stage(int iter, vector<State*>& stream) {
    const double bp = pow(2., setup_->zhtNumStages() - iter + 1);
    const double baseZT = dataFormats_->format(Variable::zT, Process::zht).base() * bp;
    const double baseCot = dataFormats_->format(Variable::cot, Process::zht).base() * bp;
    int id;
    auto different = [&id](State* state) { return state && (id != state->trackId_); };
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
      vector<vector<State*>> mhtCells(setup_->zhtNumCells());
      for (vector<State*>& mhtCell : mhtCells)
        mhtCell.reserve(size);
      // fill finer track candidates stub container
      for (auto stub = start; stub != end; stub++) {
        if (!*stub)
          continue;
        const double r = (*stub)->r_;
        const double chi = (*stub)->z_;
        const double dChi = (*stub)->dZ_;
        // identify finer track candidates for this stub
        // 0 and 1 belong to the ZHT cells with smaller cot; 0 and 2 belong to those with smaller zT
        TTBV cells(0, setup_->zhtNumCells());
        const bool compA = 2. * abs(chi) < baseZT + abs(r) * baseCot + dChi;
        const bool compB = 2. * abs(chi) < baseZT + dChi;
        const bool compC = 2. * abs(chi) < abs(r) * baseCot + dChi;
        const bool compD = 2. * abs(chi) < dChi;
        if (chi >= 0. && r >= 0.) { if (compA) cells.set(3); if (compB) cells.set(1); if (compC) cells.set(2); if (compD) cells.set(0); }
        if (chi >= 0. && r <  0.) { if (compA) cells.set(1); if (compB) cells.set(3); if (compC) cells.set(0); if (compD) cells.set(2); }
        if (chi <  0. && r >= 0.) { if (compA) cells.set(0); if (compB) cells.set(2); if (compC) cells.set(1); if (compD) cells.set(3); }
        if (chi <  0. && r <  0.) { if (compA) cells.set(2); if (compB) cells.set(0); if (compC) cells.set(3); if (compD) cells.set(1); }
        /*cout << setw(12) << r << " "
             << setw(12) << chi << " "
             << setw(12) << dChi << " "
             << "| "
             << setw(12) << baseCot << " "
             << setw(12) << baseZT << " "
             << "| "
             << setw(2) << (*stub)->cot_ << " "
             << setw(2) << (*stub)->zT_ << " "
             << "| ";
        for (int cell : cells.ids())
          cout << cell << " ";
        cout << endl;*/
        for (int cell : cells.ids())
          mhtCells[cell].push_back(*stub);
      }
      // perform pattern recognition
      TTBV cells(0, setup_->zhtNumCells());
      vector<TTBV> hitPatternsPS(setup_->zhtNumCells(), TTBV(0, setup_->numLayers()));
      for (int cell = 0; cell < setup_->zhtNumCells(); cell++) {
        const vector<State*>& mhtCell = mhtCells[cell];
        TTBV hitPattern(0, setup_->numLayers());
        for (State* stub : mhtCell)
          hitPattern.set(stub->layer_);
        if (hitPattern.count() < setup_->zhtMinLayers())
          continue;
        cells.set(cell);
        TTBV& hitPatternPS = hitPatternsPS[cell];
        for (State* stub : mhtCell)
          if (setup_->psModule(stub->stub_->ttStubRef()))
            hitPatternPS.set(stub->layer_);
      }
      // excludes cells outisde of eta sector or z0 fiducial range
      /*for (int cell : cells.ids()) {
        const int dcot = (cell / setup_->zhtNumBinsZT() - setup_->zhtNumBinsCot() + 1) * p;
        const int dzT = (cell % setup_->zhtNumBinsZT() - setup_->zhtNumBinsZT() + 1) * p;
        StubZHT* stub = mhtCells[cell].front();
        const int icot = stub->cot() + dcot;
        const int izT = stub->zT() + dzT;
        const double cot = dataFormats_->format(Variable::cot, Process::zht).floating(icot);
        const double zT = dataFormats_->format(Variable::zT, Process::zht).floating(izT);
        const double sectorzT = (sinh(setup_->boundarieEta(stub->sectorEta() + 1)) - sinh(setup_->boundarieEta(stub->sectorEta()))) * setup_->chosenRofZ();
        const double z0 = zT - setup_->chosenRofZ() * cot;
        cout << sectorzT << " " << zT << " " << z0 << endl;
        if (abs(zT) > sectorzT / 2. || abs(z0) > setup_->beamWindowZ())
          cells.reset(cell);
      }*/
      // only take cells with max number of ps layer into account
      vector<TTBV> hitPatterns(setup_->zhtNumCells(), TTBV(0, setup_->numLayers()));
      const int maxLayerPS = max_element(hitPatternsPS.begin(), hitPatternsPS.end(), moreLayer)->count();
      for (int cell : cells.ids()) {
        if (hitPatternsPS[cell].count() < maxLayerPS)
          cells.reset(cell);
        else
          for (State* stub : mhtCells[cell])
            hitPatterns[cell].set(stub->layer_);
      }
      const int maxLayer = max_element(hitPatterns.begin(), hitPatterns.end(), moreLayer)->count();
      for (int cell : cells.ids())
        if (hitPatterns[cell].count() < maxLayer)
          cells.reset(cell);
      // identify merged finer track
      set<State*> track;
      pair<int, int> cots(setup_->zhtNumBinsCot(), -1);
      pair<int, int> zTs(setup_->zhtNumBinsZT(), -1);
      for (int cell : cells.ids()) {
        const int cot = cell / setup_->zhtNumBinsZT();
        const int zT = cell % setup_->zhtNumBinsZT();
        cots = {min(cots.first, cot), max(cots.second, cot)};
        zTs = {min(zTs.first, zT), max(zTs.second, zT)};
        for (State* stub : mhtCells[cell])
          track.insert(stub);
      }
      const int p = pow(2, setup_->zhtNumStages() - iter);
      const int cot = (cots.first + cots.second - setup_->zhtNumBinsCot() + 1) * p / 2;
      const int zT = (zTs.first + zTs.second - setup_->zhtNumBinsZT() + 1) * p / 2;
      /*if (cells.any())
        cout << setw(2) << cot << " " << setw(2) << zT << endl;*/
      // second loop over track, kill bad stubs and update track parameter and stub residual
      for (auto& stub = start; stub != end; stub++) {
        if (!*stub)
          continue;
        const auto s = find(track.begin(), track.end(), *stub);
        if (s == track.end())
          *stub = nullptr;
        else
          (*stub)->update(cot, zT, dataFormats_);
      }
    }
    stream.erase(remove(stream.begin(), stream.end(), nullptr), stream.end());
    /*if (!stream.empty())
      cout << endl;*/
  }

}  // namespace trackerTFP
