#include "L1Trigger/TrackFindingTracklet/interface/TrackQuality.h"

#include <vector>
#include <string>
#include <numeric>

#include "conifer.h"
#include "ap_fixed.h"

namespace trklet {

  // read in and organize input tracks and stubs
  void TrackQuality::consume(const tt::StreamsTrack& tracks, const tt::StreamsStub& stubs) {
    streams_ = tracks;
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
  void TrackQuality::produce(tt::StreamsTrack& outputs) const {
    const int offset = setup_->tqNumChannel() * region_;
    outputs[offset + 0] = streams_[region_];
    tt::StreamTrack& output = outputs[offset + 1];
    const DataFormat& dfChi20 = dataFormats_->format(Variable::chi20, Process::tq);
    const DataFormat& dfChi21 = dataFormats_->format(Variable::chi21, Process::tq);
    output.reserve(input_.size());
    for (const Frame& frame : input_) {
      if (!frame.track_) {
        output.emplace_back(tt::FrameTrack());
        continue;
      }
      // analyze track and stubs
      double chi20(0.);
      double chi21(0.);
      TTBV hitPattern(0, setup_->numLayers());
      for (int layer = 0; layer < setup_->numLayers(); layer++) {
        StubKF* stub = frame.stubs_[layer];
        if (!stub)
          continue;
        hitPattern.set(layer);
        const double m02 = internalFormats_->m02_.digi(std::pow(stub->phi(), 2));
        const double m12 = internalFormats_->m12_.digi(std::pow(stub->z(), 2));
        const double invV0 = internalFormats_->invV0_.digi(1. / std::pow(2. * stub->dPhi(), 2));
        const double invV1 = internalFormats_->invV1_.digi(1. / std::pow(2. * stub->dZ(), 2));
        chi20 += dfChi20.limit(dfChi20.digi(m02 * invV0));
        chi21 += dfChi21.limit(dfChi21.digi(m12 * invV1));
      }

      // Accumulating all BDT Attributes //
      chi20 = dfChi20.limit(chi20);
      chi21 = dfChi21.limit(chi21);
      double zT = (*frame.track_).zT();
      double cot = (*frame.track_).cot();
      int nstub = hitPattern.count();
      int n_missint = get_ninterior(hitPattern);
      // BDT Inference //
      const AP_FIXED_BDT f_nstub    = (AP_FIXED_BDT)nstub;
      const AP_FIXED_BDT f_zT       = transform_zT(zT);
      const AP_FIXED_BDT f_cot      = transform_cot(cot);
      const AP_FIXED_BDT f_chi20    = (AP_FIXED_BDT)chi20;
      const AP_FIXED_BDT f_chi21    = (AP_FIXED_BDT)chi21;
      const AP_FIXED_BDT f_n_miss   = (AP_FIXED_BDT)n_missint;
      const std::vector<AP_FIXED_BDT> features = 
      {
       f_nstub, f_zT, f_cot, f_chi20, f_chi21, f_n_miss
      };
      const AP_FIXED_BDT& mva_raw = bdt_->decision_function(features).at(0);
      const AP_INT_BDT& mva_ = mva_raw.range(mva_raw.width - 1, 0);
      const int mva = packMVA(mva_);
      // build output Track
      TrackTQ trackTQ(*frame.track_, chi20, chi21, mva, hitPattern);
      // store result
      output.push_back(trackTQ.frame());
    }
  }

  const int TrackQuality::packMVA(const AP_INT_BDT& mva) const {
    int out = 0;
    for (int i = 0; i < 8; ++i) {
      if ((mva > kMVABinEdges[i]) && (mva <= kMVABinEdges[i + 1])) {
        out = i;
        break;
      }
    }
    return out;
  }

  const  TrackQuality::AP_FIXED_BDT TrackQuality::transform_zT(const float& zT) const {
      int zT_vivado_view_int = std::floor(187.5351440760983 * zT);
      AP_FIXED_BDT zT_fixed; 
      const AP_INT_BDT zT_int (zT_vivado_view_int);
      zT_fixed.range(19, 0) = zT_int.range(19, 0);
      return zT_fixed;
  }

  const TrackQuality::AP_FIXED_BDT TrackQuality::transform_cot(const float& cot) const {
    int cot_vivado_view_int = std::floor(4869.970502268804 * cot);
    AP_FIXED_BDT cot_fixed;
    const AP_INT_BDT cot_int(cot_vivado_view_int);
    cot_fixed.range(19, 0) = cot_int.range(19, 0);
    return cot_fixed;
  }

  const int TrackQuality::get_ninterior(const TTBV& hitPattern) const {
    std::string s = hitPattern.str();
    std::reverse(s.begin(), s.end());
    if (!s.empty())
      s.pop_back();
    const size_t first_one = s.find('1');
    const size_t last_one  = s.rfind('1');
    if (first_one == std::string::npos || first_one == last_one)
      return 0;
    int count = 0;
    for (size_t i = first_one + 1; i < last_one; ++i) {
      if (s[i] == '0')
        ++count;
    }
  return count;
  }

}  // namespace trklet
