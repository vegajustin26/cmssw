#include "L1Trigger/TrackerTFP/interface/SeedFilter.h"

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
using namespace trackerDTC;

namespace trackerTFP {

  SeedFilter::SeedFilter(const ParameterSet& iConfig, const Setup* setup, const DataFormats* dataFormats, const LayerEncoding* layerEncoding, int region) :
    enableTruncation_(iConfig.getParameter<bool>("EnableTruncation")),
    setup_(setup),
    dataFormats_(dataFormats),
    layerEncoding_(layerEncoding),
    region_(region),
    cot_(dataFormats_->format(Variable::cot, Process::sf)),
    zT_(dataFormats_->format(Variable::zT, Process::sf)),
    z_(dataFormats_->format(Variable::z, Process::sf)),
    input_(dataFormats_->numChannel(Process::mht)) {}

  void SeedFilter::consume(const TTDTC::Streams& streams) {
    auto valid = [](int& sum, const TTDTC::Frame& frame){ return sum += (frame.first.isNonnull() ? 1 : 0); };
    int nStubs(0);
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      const TTDTC::Stream& stream = streams[region_ * dataFormats_->numChannel(Process::mht) + channel];
      nStubs += accumulate(stream.begin(), stream.end(), 0, valid);
    }
    stubsMHT_.reserve(nStubs);
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      const TTDTC::Stream& stream = streams[region_ * dataFormats_->numChannel(Process::mht) + channel];
      vector<StubMHT*>& stubs = input_[channel];
      stubs.reserve(stream.size());
      for (const TTDTC::Frame& frame : stream) {
        StubMHT* stub = nullptr;
        if (frame.first.isNonnull()) {
          stubsMHT_.emplace_back(frame, dataFormats_);
          stub = &stubsMHT_.back();
        }
        stubs.push_back(stub);
      }
    }
  }

  void SeedFilter::produce(TTDTC::Streams& accepted, TTDTC::Streams& lost) {
    auto chi2 = [this](const vector<StubSF*>& stubs) {
      double d(0.);
      TTBV hitPattern(0, setup_->numLayers());
      StubSF* stub = stubs.front();
      const vector<int>& layerEncoding = layerEncoding_->layerEncoding(stub->sectorEta(), zT_.toUnsigned(stub->zT()), cot_.toUnsigned(stub->cot()));
      for (StubSF* stub : stubs) {
        hitPattern.set(distance(layerEncoding.begin(), find(layerEncoding.begin(), layerEncoding.end(), stub->layer())));
        d += stub->chi2();
      }
      return d / (hitPattern.count() - 2) * pow(2, hitPattern.count(0, hitPattern.pmEncode(), false));
    };
    auto smallerChi2 = [this, chi2](const vector<StubSF*>& lhs, const vector<StubSF*>& rhs) {
      return chi2(lhs) < chi2(rhs);
    };
    /*auto smallerChi = [](const vector<StubSF*>& lhs, const vector<StubSF*>& rhs) {
      auto chi = [](const vector<StubSF*>& stubs) {
        double c(0.);
        for (StubSF* stub : stubs)
          c += abs(stub->z());
        return c;
      };
      return chi(lhs) < chi(rhs);
    };
    auto moreLayer = [](const vector<StubSF*>& lhs, const vector<StubSF*>& rhs) {
      auto nLayer = [](const vector<StubSF*>& stubs) {
        set<int> layer;
        for (StubSF* stub : stubs)
          layer.insert(stub->layer());
        return layer.size();
      };
      return nLayer(lhs) > nLayer(rhs);
    };*/
    const int numLayers = 16;
    static const vector<pair<int, int>> seedingLayers = {{1, 2}, {1, 3}, {2, 3}, {1, 11}, {2, 11}, {1, 12}, {2, 12}, {11, 12}};
    const double dT = setup_->chosenRofPhi() - setup_->chosenRofZ();
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      TTDTC::Stream& acceptedStream = accepted[region_ * dataFormats_->numChannel(Process::mht) + channel];
      TTDTC::Stream& lostStream = lost[region_ * dataFormats_->numChannel(Process::mht) + channel];
      deque<StubSF*> output;
      vector<StubMHT*>& stubs = input_[channel];
      stubs.erase(remove(stubs.begin(), stubs.end(), nullptr), stubs.end());
      set<int> ids;
      for (StubMHT* stub : stubs)
        if (stub)
          ids.insert(stub->trackId());
      vector<vector<StubMHT*>> tracks(ids.size());
      int i(0);
      for (auto it = stubs.begin(); it != stubs.end();) {
        const auto start = it;
        const int id = (*it)->trackId();
        auto different = [id](StubMHT* stub){ return id != stub->trackId(); };
        it = find_if(it, stubs.end(), different);
        vector<StubMHT*>& track = tracks[i++];
        track.reserve(distance(start, it));
        copy(start, it, back_inserter(track));
      }
      int iTrack(0);
      const int sectorPhi = 5;
      const int trackId = 192;
      const bool debug = true && region_ == sectorPhi / 2 && channel == trackId / setup_->sfMaxTracks() && iTrack++ == trackId % setup_->sfMaxTracks();
      for (const vector<StubMHT*>& track : tracks) {
        if (debug) {
          cout << "SF" << endl;
          for (StubMHT* stub : track) {
            const GlobalPoint& gp = setup_->stubPos(stub->ttStubRef());
            cout << gp.perp() << " " << gp.phi() << " " << gp.z() << endl;
          }
          cout << endl;
        }
        const int eta = track.front()->sectorEta();
        const double dZT = (sinh(setup_->boundarieEta(eta + 1)) - sinh(setup_->boundarieEta(eta))) * setup_->chosenRofZ();
        vector<vector<StubMHT*>> layerStubs(numLayers);
        for (vector<StubMHT*>& layer : layerStubs)
          layer.reserve(track.size());
        for (StubMHT* stub : track)
          layerStubs[layerId(stub)].push_back(stub);
        deque<pair<StubMHT*, StubMHT*>> seeds;
        for (const pair<int, int>& seedingLayer : seedingLayers)
          for (StubMHT* inner : layerStubs[seedingLayer.first])
            if (setup_->psModule(inner->ttStubRef()))
              for (StubMHT* outer : layerStubs[seedingLayer.second])
                if (setup_->psModule(outer->ttStubRef()))
                  if (outer->r() > inner->r())
                    seeds.emplace_back(inner, outer);
        vector<vector<StubSF*>> seedTracks;
        seedTracks.reserve(seeds.size());
        for (const pair<StubMHT*, StubMHT*>& seed : seeds) {
          const double r1 = seed.first->r() + dT;
          const double r2 = seed.second->r() + dT;
          const double z1 = seed.first->z();
          const double z2 = seed.second->z();
          const double r = (r2 + r1) / 2.;
          const double z = (z2 + z1) / 2.;
          const double dr = (r2 - r1);
          const double dz = (z2 - z1);
          const double cot = cot_.digi(dz / dr);
          const double zT = zT_.digi(z - cot * r);
          const double z0 = zT - cot * setup_->chosenRofZ();
          if (abs(z0) > setup_->beamWindowZ() + cot_.base() * setup_->chosenRofZ())
          //if (abs(z0) > setup_->beamWindowZ())
            continue;
          if (abs(zT) > dZT / 2. + zT_.base() / 2.)
          //if (abs(zT) > dZT / 2.)
            continue;
          vector<StubSF*> stubsSF;
          stubsSF.reserve(track.size());
          for (StubMHT* stub : track) {
            const double sr = stub->r() + dT;
            const double sz = stub->z();
            const double chi = sz - (zT + sr * cot);
            const double dZ = zT_.base() + cot_.base() * abs(sr) + setup_->dZ(stub->ttStubRef());
            if (abs(chi) < dZ / 2.) {
              const double v = setup_->v1(stub->ttStubRef(), cot + setup_->sectorCot(stub->sectorEta()));
              const double cov = pow(cot_.base() * stub->r(), 2) + pow(zT_.base(), 2);
              const double chi2 = pow(chi, 2) / (v + cov);
              stubsSF_.emplace_back(*stub, zT_.integer(zT), cot_.integer(cot), chi2);
              stubsSF.push_back(&stubsSF_.back());
            }
          }
          set<int> ids;
          for (StubSF* stub : stubsSF)
            ids.insert(stub->layer());
          if ((int)ids.size() >= setup_->sfMinLayers())
            seedTracks.push_back(stubsSF);
        }
        if (seedTracks.empty())
          continue;
        stable_sort(seedTracks.begin(), seedTracks.end(), smallerChi2);
        //stable_sort(seedTracks.begin(), seedTracks.end(), moreLayer);
        const vector<StubSF*>& stubs = seedTracks.front();
        copy(stubs.begin(), stubs.end(), back_inserter(output));
        if (debug) {
          cout << cot_.base() << " " << zT_.base() << endl;
          StubSF* stub = stubs.front();
          const double qOverPt = dataFormats_->format(Variable::qOverPt, Process::mht).floating(stub->qOverPt());
          const double phiT = dataFormats_->format(Variable::phiT, Process::mht).floating(stub->phiT());
          const double off = (stub->sectorPhi() + setup_->numSectorsPhi() * region_ - .5) * dataFormats_->format(Variable::phiT, Process::ht).range();
          cout << "m0SF = " << " " << -qOverPt << endl;
          cout << "c0SF = " << " " << deltaPhi(phiT + qOverPt * setup_->chosenRofPhi() + off) << endl;
          cout << "m1SF = " << " " << cot_.floating(stub->cot()) + setup_->sectorCot(stub->sectorEta()) << endl;
          cout << "c1SF = " << " " << zT_.floating(stub->zT()) - cot_.floating(stub->cot()) * setup_->chosenRofZ() << endl;
          cout << endl;
          for (const vector<StubSF*>& track : seedTracks) {
            TTBV hitPattern(0, setup_->numLayers());
            for (StubSF* stub : track)
              hitPattern.set(stub->layer());
            cout << hitPattern << " " << chi2(track) << " ";
            double chi2(0.);
            for (StubSF* stub : track) {
              chi2 += stub->chi2();
              cout << stub->chi2() << " ";
            }
            cout << chi2 << endl;
          }
          cout << endl;
        }
      }
      if (enableTruncation_ && ((int)output.size() > setup_->numFrames())) {
        const auto limit = next(output.begin(), setup_->numFrames());
        lostStream.reserve(distance(limit, output.end()));
        transform(limit, output.end(), back_inserter(lostStream), [](StubSF* stub){ return stub->frame(); });
        output.erase(limit, output.end());
      }
      acceptedStream.reserve(output.size());
      transform(output.begin(), output.end(), back_inserter(acceptedStream), [](StubSF* stub){ return stub->frame(); });
    }
  }

  //
  int SeedFilter::layerId(StubMHT* stub) const {
    int layer = stub->layer();
    bool barrel = setup_->barrel(stub->ttStubRef());
    int lay = -1;
    if (layer == 0)
      lay = 1;
    else if (layer == 1)
      lay = 2;
    else if (layer == 2)
      lay = barrel ? 6 : 11;
    else if (layer == 3)
      lay = barrel ? 5 : 12;
    else if (layer == 4)
      lay = barrel ? 4 : 13;
    else if (layer == 5)
      lay = 14;
    else if (layer == 6)
      lay = barrel ? 3 : 15;
    return lay;
  }

  // remove and return first element of deque, returns nullptr if empty
  template<class T>
  T* SeedFilter::pop_front(deque<T*>& ts) const {
    T* t = nullptr;
    if (!ts.empty()) {
      t = ts.front();
      ts.pop_front();
    }
    return t;
  }

} // namespace trackerTFP