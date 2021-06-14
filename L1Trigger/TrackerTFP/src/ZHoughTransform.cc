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

  static constexpr int numStages = 5;
  static constexpr int numBinsZT = 2;
  static constexpr int numBinsCot = 2;
  static constexpr int numCells = numBinsZT + numBinsCot;

  ZHoughTransform::ZHoughTransform(const ParameterSet& iConfig, const Setup* setup, const DataFormats* dataFormats, int region) :
    enableTruncation_(iConfig.getParameter<bool>("EnableTruncation")),
    setup_(setup),
    dataFormats_(dataFormats),
    region_(region),
    input_(dataFormats->numChannel(Process::mht)),
    stage_(0) {}

  // read in and organize input product (fill vector input_)
  void ZHoughTransform::consume(const StreamsStub& streams) {
    auto valid = [](int& sum, const FrameStub& frame){ return sum += (frame.first.isNonnull() ? 1 : 0); };
    const int offset = region_ * dataFormats_->numChannel(Process::mht);
    int nStubsMHT(0);
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      const StreamStub& stream = streams[offset + channel];
      nStubsMHT += accumulate(stream.begin(), stream.end(), 0, valid);
    }
    stubsZHT_.reserve(nStubsMHT * ( numCells * numStages + 1 ));
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      //if (channel != 2)
        //continue;
      const StreamStub& stream = streams[offset + channel];
      vector<StubSF*>& stubs = input_[channel];
      stubs.reserve(stream.size());
      // Store input stubs in vector, so rest of ZHT algo can work with pointers to them (saves CPU)
      for (const FrameStub& frame : stream) {
        StubSF* stub = nullptr;
        if (frame.first.isNonnull()) {
          StubMHT stubMHT(frame, dataFormats_);
          stubsZHT_.emplace_back(stubMHT, stubMHT.layer(), 0, 0, 0.);
          stub = &stubsZHT_.back();
        }
        stubs.push_back(stub);
      }
    }
    //sort(stubsZHT_.begin(), stubsZHT_.end(), [](const StubSF& lhs, const StubSF& rhs){ return lhs.z() > rhs.z(); });
  }

  // fill output products
  void ZHoughTransform::produce(StreamsStub& accepted, StreamsStub& lost) {
    //if (region_ != 1)
      //return;
    vector<deque<StubSF*>> streams(dataFormats_->numChannel(Process::mht));
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++)
      streams[channel] = deque<StubSF*>(input_[channel].begin(), input_[channel].end());
    vector<deque<StubSF*>> stubsCells(dataFormats_->numChannel(Process::mht) * numCells);
    for (stage_ = 0; stage_ < numStages; stage_++) {
      // fill ZHT cells
      for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++)
        fill(channel, streams[channel], stubsCells);
      // perform static load balancing
      for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
        vector<deque<StubSF*>> tmp(numCells);
        // gather streams to mux together: same ZHT cell of 4 adjacent ZHT input streams
        for (int k = 0; k < numCells; k++)
          //swap(tmp[k], stubsCells[(channel / numCells) * dataFormats_->numChannel(Process::mht) + channel % numCells + k * numCells]);
          swap(tmp[k], stubsCells[channel * numCells + k]);
        slb(tmp, streams[channel], lost[channel]);
      }
    }
    // fill output product
    for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      deque<StubSF*>& stubs = streams[channel];
      StreamStub& stream = accepted[region_ * dataFormats_->numChannel(Process::mht) + channel];
      merge(stubs, stream);
    }
    // fill output product
    /*for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
      const deque<StubSF*>& stubs = streams[channel];
      StreamStub& stream = accepted[region_ * dataFormats_->numChannel(Process::mht) + channel];
      stream.reserve(stubs.size());
      for (StubSF* stub : stubs)
        stream.emplace_back(stub ? stub->frame() : FrameStub());
      //for (StubSF* stub : stubs)
        //if (stub)
          //cout << stub->cot() << " " << stub->zT() << " " << stub->r() << " " << stub->z() << endl;
    }*/
  }

  // perform finer pattern recognition per track
  void ZHoughTransform::fill(int channel, const deque<StubSF*>& stubs, vector<deque<StubSF*>>& streams) {
    if (stubs.empty())
      return;
    const double baseZT = dataFormats_->format(Variable::zT, Process::sf).base() * pow(2, numStages - stage_);
    const double baseCot = dataFormats_->format(Variable::cot, Process::sf).base() * pow(2, numStages - stage_);
    int id;
    auto different = [&id](StubSF* stub){ return !stub || id != stub->trackId(); };
    for (auto it = stubs.begin(); it != stubs.end();) {
      if (!*it) {
        const auto begin = find_if(it, stubs.end(), [](StubSF* stub){ return stub; });
        const int nGaps = distance(it, begin);
        for (deque<StubSF*>& stream : streams)
          stream.insert(stream.end(), nGaps, nullptr);
        it = begin;
        continue;
      }
      const auto start = it;
      const double cotGlobal = (*start)->cotf() + setup_->sectorCot((*start)->sectorEta());
      id = (*it)->trackId();
      it = find_if(it, stubs.end(), different);
      const int size = distance(start, it);
      // create finer track candidates stub container
      vector<vector<StubSF*>> mhtCells(numCells);
      for (vector<StubSF*>& mhtCell : mhtCells)
        mhtCell.reserve(size);
      // fill finer track candidates stub container
      for (auto stub = start; stub != it; stub++) {
        const double r = (*stub)->r() + setup_->chosenRofPhi() - setup_->chosenRofZ();
        const double chi = (*stub)->z();
        const double dChi = setup_->dZ((*stub)->ttStubRef(), cotGlobal);
        //const double dChi = 0.;
        //cout << region_ << " " << stage_ << " " << (*stub)->cot() << " " << (*stub)->zT() << " " << r << " " << chi << " " << dChi << " ";
        // identify finer track candidates for this stub
        // 0 and 1 belong to the ZHT cells with smaller cot; 0 and 2 belong to those with smaller zT
        vector<int> cells;
        cells.reserve(numCells);
        const bool compA = 2. * abs(chi) < baseZT + dChi;
        const bool compB = 2. * abs(chi) < abs(r) * baseCot + dChi;
        const bool compC = 2. * abs(chi) < dChi;
        if (chi >= 0. && r >= 0.) { cells.push_back(1); if (compA) cells.push_back(3); if(compB) cells.push_back(0); if(compC) cells.push_back(2); }
        if (chi >= 0. && r <  0.) { cells.push_back(3); if (compA) cells.push_back(1); if(compB) cells.push_back(2); if(compC) cells.push_back(0); }
        if (chi <  0. && r >= 0.) { cells.push_back(2); if (compA) cells.push_back(0); if(compB) cells.push_back(3); if(compC) cells.push_back(1); }
        if (chi <  0. && r <  0.) { cells.push_back(0); if (compA) cells.push_back(2); if(compB) cells.push_back(1); if(compC) cells.push_back(3); }
        for (int cell : cells) {
          const double cot = (cell / numBinsZT - .5) * baseCot / 2.;
          const double zT = (cell % numBinsZT - .5) * baseZT / 2.;
          stubsZHT_.emplace_back(**stub, zT, cot, cell);
          mhtCells[cell].push_back(&stubsZHT_.back());
          //cout << cell << " ";
        }
        //cout << endl;
      }
      //cout << endl;
      // perform pattern recognition
      for (int sel = 0; sel < numCells; sel++) {
        deque<StubSF*>& stream = streams[channel * numCells + sel];
        vector<StubSF*>& mhtCell = mhtCells[sel];
        set<int> layers;
        auto toLayer = [](StubSF* stub){ return stub->layer(); };
        transform(mhtCell.begin(), mhtCell.end(), inserter(layers, layers.begin()), toLayer);
        if ((int)layers.size() < setup_->mhtMinLayers())
          mhtCell.clear();
        for (StubSF* stub : mhtCell)
          stream.push_back(stub);
        stream.insert(stream.end(), size - (int)mhtCell.size(), nullptr);
      }
    }
    for (int sel = 0; sel < numCells; sel++) {
      deque<StubSF*>& stream = streams[channel * numCells + sel];
      // remove all gaps between end and last stub
      for(auto it = stream.end(); it != stream.begin();)
        it = (*--it) ? stream.begin() : stream.erase(it);
      // read out fine track cannot start before rough track has read in completely, add gaps to take this into account
      int pos(0);
      for (auto it = stream.begin(); it != stream.end();) {
        if (!(*it)) {
          it = stream.erase(it);
          continue;
        }
        id = (*it)->trackId();
        const int s = distance(it, find_if(it, stream.end(), different));
        const int d = distance(stream.begin(), it);
        pos += s;
        if (d < pos) {
          const int diff = pos - d;
          it = stream.insert(it, diff, nullptr);
          it = next(it, diff);
        }
        else
          it = stream.erase(remove(next(stream.begin(), pos), it, nullptr), it);
        it = next(it, s);
      }
      // adjust stream start so that first output stub is in first place in case of quickest track
      if (!stream.empty())
        stream.erase(stream.begin(), next(stream.begin(), setup_->mhtMinLayers()));
    }
  }

  // Static load balancing of inputs: mux 4 streams to 1 stream
  void ZHoughTransform::slb(vector<deque<StubSF*>>& inputs, deque<StubSF*>& accepted, StreamStub& lost) const {
    accepted.clear();
    if (all_of(inputs.begin(), inputs.end(), [](const deque<StubSF*>& stubs){ return stubs.empty(); }))
      return;
    // input fifos
    vector<deque<StubSF*>> stacks(numCells);
    // helper for handshake
    TTBV empty(-1, numCells, true);
    TTBV enable(0, numCells);
    // clock accurate firmware emulation, each while trip describes one clock tick, one stub in and one stub out per tick
    while(!all_of(inputs.begin(), inputs.end(), [](const deque<StubSF*>& d){ return d.empty(); }) or
          !all_of(stacks.begin(), stacks.end(), [](const deque<StubSF*>& d){ return d.empty(); })) {
      // store stub in fifo
      for(int channel = 0; channel < numCells; channel++){
        StubSF* stub = pop_front(inputs[channel]);
        if (stub)
          stacks[channel].push_back(stub);
      }
      // identify empty fifos
      for (int channel = 0; channel < numCells; channel++)
        empty[channel] = stacks[channel].empty();
      // chose new fifo to read from if current fifo got empty
      const int iEnableOld = enable.plEncode();
      if (enable.none() || empty[iEnableOld]) {
        enable.reset();
        const int iNotEmpty = empty.plEncode(false);
        if (iNotEmpty < numCells)
          enable.set(iNotEmpty);
      }
      // read from chosen fifo
      const int iEnable = enable.plEncode();
      if (enable.any())
        accepted.push_back(pop_front(stacks[iEnable]));
      else
        // gap if no fifo has been chosen
        accepted.push_back(nullptr);

    }
    // perform truncation if desired
    if (enableTruncation_ && (int)accepted.size() > setup_->numFrames()) {
      const auto limit = next(accepted.begin(), setup_->numFrames());
      auto valid = [](int& sum, StubSF* stub){ return sum += stub ? 1 : 0; };
      const int nLost = accumulate(limit, accepted.end(), 0, valid);
      lost.reserve(nLost);
      for (auto it = limit; it != accepted.end(); it++)
        if (*it)
          lost.emplace_back((*it)->frame());
      accepted.erase(limit, accepted.end());
    }
    // cosmetics -- remove gaps at the end of stream
    for(auto it = accepted.end(); it != accepted.begin();)
      it = (*--it) == nullptr ? accepted.erase(it) : accepted.begin();
  }

  //
  void ZHoughTransform::merge(deque<StubSF*>& stubs, StreamStub& stream) const {
    stubs.erase(remove(stubs.begin(), stubs.end(), nullptr), stubs.end());
    set<int> candidates;
    const int weight = numCells * pow(2, numStages);
    for (const StubSF* stub : stubs)
      candidates.insert(stub->trackId() / weight);
    vector<deque<StubSF*>> tracks(candidates.size());
    for (auto it = stubs.begin(); it != stubs.end();) {
      const auto start = it;
      const int id = (*it)->trackId();
      const int candId = id / weight;
      const int pos = distance(candidates.begin(), find(candidates.begin(), candidates.end(), candId));
      deque<StubSF*>& track = tracks[pos];
      auto different = [id](const StubSF* stub){ return id != stub->trackId(); };
      it = find_if(it, stubs.end(), different);
      for (auto s = start; s != it; s++)
        if (find_if(track.begin(), track.end(), [s](const StubSF* stub){ return (*s)->ttStubRef() == stub->ttStubRef(); }) == track.end())
          track.push_back(*s);
    }
    const int size = accumulate(tracks.begin(), tracks.end(), 0, [](int& sum, const deque<StubSF*>& stubs){ return sum += (int)stubs.size(); });
    stream.reserve(size);
    for (deque<StubSF*>& track : tracks) {
      const double cot = (*track.begin())->cotf();
      const double zT = (*track.begin())->ztf();
      for (const StubSF* stub : track) {
        const StubSF stubZHT(*stub, zT, cot);
        stream.push_back(stubZHT.frame());
      }
    }
  }

  // remove and return first element of deque, returns nullptr if empty
  template<class T>
  T* ZHoughTransform::pop_front(deque<T*>& ts) const {
    T* t = nullptr;
    if (!ts.empty()) {
      t = ts.front();
      ts.pop_front();
    }
    return t;
  }

} // namespace trackerTFP