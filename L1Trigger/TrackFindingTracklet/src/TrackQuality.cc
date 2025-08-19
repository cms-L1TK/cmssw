#include "L1Trigger/TrackFindingTracklet/interface/TrackQuality.h"

#include <vector>
#include <string>
#include <numeric>

namespace trklet {

  // read in and organize input tracks and stubs
  void TrackQuality::consume(const tt::StreamsTrack& tracks, const tt::StreamsStub& stubs) {
    auto validT = [](int sum, const tt::FrameTrack& f) { return sum + (f.first.isNull() ? 0 : 1); };
    auto validS = [](int sum, const tt::FrameStub& f) { return sum + (f.first.isNull() ? 0 : 1); };
    const int offset = region_ * setup_->numLayers();
    const tt::StreamTrack& input = tracks[region_];
    // count tracks
    const int nTracks = std::accumulate(input.begin(), input.end(), 0, validT);
    tracks_.reserve(nTracks);
    // count stubs
    int nStubs(0);
    for (int iLayer = 0; iLayer < setup_->numLayers(); iLayer++) {
      const tt::StreamStub& layer = stubs[offset + iLayer];
      nStubs += std::accumulate(layer.begin(), layer.end(), 0, validS);
    }
    stubs_.reserve(nStubs);
    // store input
    input_.reserve(input.size());
    for (int iFrame = 0; iFrame < static_cast<int>(input.size()); iFrame++) {
      const tt::FrameTrack& frameTrack = input[iFrame];
      if (frameTrack.first.isNull()) {
        input_.emplace_back(Frame());
        continue;
      }
      tracks_.emplace_back(frameTrack, dataFormats_);
      input_.emplace_back(&tracks_.back(), setup_->numLayers());
      for (int iLayer = 0; iLayer < setup_->numLayers(); iLayer++) {
        const tt::FrameStub& frameStub = stubs[offset + iLayer][iFrame];
        if (frameStub.first.isNull())
          continue;
        stubs_.emplace_back(frameStub, dataFormats_);
        input_.back().stubs_[iLayer] = &stubs_.back();
      }
    }
  }

  // fills output products
  void TrackQuality::produce(tt::StreamTrack& output) const {
    const DataFormat& dfChi20 = dataFormats_->format(Variable::chi20, Process::tq);
    const DataFormat& dfChi21 = dataFormats_->format(Variable::chi21, Process::tq);
    output.reserve(input_.size());
    for (const Frame& frame : input_) {
      if (!frame.track_) {
        output.emplace_back(tt::FrameTrack());
        continue;
      }
      // analyze stubs
      TTBV hitPattern(0, setup_->numLayers() + 1);
      double chi20(0.);
      double chi21(0.);
      for (int layer = 0; layer < setup_->numLayers(); layer++) {
        StubKF* stub = frame.stubs_[layer];
        if (!stub)
          continue;
        hitPattern.set(layer);
        const double m02 = internalFormats_->m02_.digi(std::pow(stub->phi(), 2));
        const double m12 = internalFormats_->m12_.digi(std::pow(stub->z(), 2));
        const double invV0 = internalFormats_->invV0_.digi(1. / std::pow(2. * stub->dPhi(), 2));
        const double invV1 = internalFormats_->invV1_.digi(1. / std::pow(2. * stub->dZ(), 2));
        chi20 += dfChi20.digi(m02 * invV0);
        chi21 += dfChi21.digi(m12 * invV1);
      }
      chi20 = dfChi20.limit(chi20);
      chi21 = dfChi21.limit(chi21);
      // emulate bdt
      //// to be done ////
      // build output Track
      TrackTQ trackTQ(*frame.track_, chi20, chi21, 0, hitPattern);
      // store result
      output.push_back(trackTQ.frame());
    }
  }

}  // namespace trklet
