#include "L1Trigger/TrackerTFP/interface/KalmanFilter.h"

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

  KalmanFilter::KalmanFilter(const ParameterSet& iConfig, const Setup* setup, const DataFormats* dataFormats, KalmanFilterFormats* kalmanFilterFormats, int region, int numChannel) :
    enableTruncation_(iConfig.getParameter<bool>("EnableTruncation")),
    setup_(setup),
    dataFormats_(dataFormats),
    kalmanFilterFormats_(kalmanFilterFormats),
    region_(region),
    inputStubs_(numChannel * setup_->numLayers()),
    inputTracks_(numChannel),
    channelStubs_(numChannel),
    layer_(0),
    x0_(&kalmanFilterFormats_->format(VariableKF::x0)),
    x1_(&kalmanFilterFormats_->format(VariableKF::x1)),
    x2_(&kalmanFilterFormats_->format(VariableKF::x2)),
    x3_(&kalmanFilterFormats_->format(VariableKF::x3)),
    H00_(&kalmanFilterFormats_->format(VariableKF::H00)),
    H12_(&kalmanFilterFormats_->format(VariableKF::H12)),
    m0_(&kalmanFilterFormats_->format(VariableKF::m0)),
    m1_(&kalmanFilterFormats_->format(VariableKF::m1)),
    v0_(&kalmanFilterFormats_->format(VariableKF::v0)),
    v1_(&kalmanFilterFormats_->format(VariableKF::v1)),
    r0_(&kalmanFilterFormats_->format(VariableKF::r0)),
    r1_(&kalmanFilterFormats_->format(VariableKF::r1)),
    S00_(&kalmanFilterFormats_->format(VariableKF::S00)),
    S01_(&kalmanFilterFormats_->format(VariableKF::S01)),
    S12_(&kalmanFilterFormats_->format(VariableKF::S12)),
    S13_(&kalmanFilterFormats_->format(VariableKF::S13)),
    K00_(&kalmanFilterFormats_->format(VariableKF::K00)),
    K10_(&kalmanFilterFormats_->format(VariableKF::K10)),
    K21_(&kalmanFilterFormats_->format(VariableKF::K21)),
    K31_(&kalmanFilterFormats_->format(VariableKF::K31)),
    R00_(&kalmanFilterFormats_->format(VariableKF::R00)),
    R11_(&kalmanFilterFormats_->format(VariableKF::R11)),
    R00Rough_(&kalmanFilterFormats_->format(VariableKF::R00Rough)),
    R11Rough_(&kalmanFilterFormats_->format(VariableKF::R11Rough)),
    invR00Approx_(&kalmanFilterFormats_->format(VariableKF::invR00Approx)),
    invR11Approx_(&kalmanFilterFormats_->format(VariableKF::invR11Approx)),
    invR00Cor_(&kalmanFilterFormats_->format(VariableKF::invR00Cor)),
    invR11Cor_(&kalmanFilterFormats_->format(VariableKF::invR11Cor)),
    invR00_(&kalmanFilterFormats_->format(VariableKF::invR00)),
    invR11_(&kalmanFilterFormats_->format(VariableKF::invR11)),
    C00_(&kalmanFilterFormats_->format(VariableKF::C00)),
    C01_(&kalmanFilterFormats_->format(VariableKF::C01)),
    C11_(&kalmanFilterFormats_->format(VariableKF::C11)),
    C22_(&kalmanFilterFormats_->format(VariableKF::C22)),
    C23_(&kalmanFilterFormats_->format(VariableKF::C23)),
    C33_(&kalmanFilterFormats_->format(VariableKF::C33)) {}

  // read in and organize input stubs
  void KalmanFilter::consume(const TTDTC::Streams& stubs, const TTDTC::Streams& lost) {
    auto valid = [](int& sum, const TTDTC::Frame& frame){ return sum += (frame.first.isNonnull() ? 1 : 0); };
    int nStubs(0);
    for (int channel = 0; channel < (int)inputStubs_.size(); channel++) {
      const TTDTC::Stream& stream = stubs[region_ * inputStubs_.size() + channel];
      const TTDTC::Stream& streamLost = lost[region_ * inputStubs_.size() + channel];
      nStubs += accumulate(stream.begin(), stream.end(), 0, valid);
      nStubs += accumulate(streamLost.begin(), streamLost.end(), 0, valid);
    }
    stubs_.reserve(nStubs);
    for (int channel = 0; channel < (int)inputStubs_.size(); channel++) {
      const int layerId = channel % setup_->numLayers();
      const TTDTC::Stream& stream = stubs[region_ * inputStubs_.size() + channel];
      const TTDTC::Stream& streamLost = lost[region_ * inputStubs_.size() + channel];
      vector<StubKFin*>& input = inputStubs_[channel];
      input.reserve(stream.size());
      for (const TTDTC::Frame& frame : stream) {
        StubKFin* stub = nullptr;
        if (frame.first.isNonnull()) {
          stubs_.emplace_back(frame, dataFormats_, layerId);
          stub = &stubs_.back();
        }
        input.push_back(stub);
      }
      for (const TTDTC::Frame& frame : streamLost) {
        stubs_.emplace_back(frame, dataFormats_, layerId);
        input.push_back(&stubs_.back());
      }
    }
    for (int channel = 0; channel < (int)inputTracks_.size(); channel++) {
      const int offset = channel * setup_->numLayers();
      vector<StubKFin*>& stubs = channelStubs_[channel];
      int nStubs(0);
      for (int layer = 0; layer < setup_->numLayers(); layer++) {
        const vector<StubKFin*>& stream = inputStubs_[offset + layer];
        nStubs += accumulate(stream.begin(), stream.end(), 0, [](int& sum, StubKFin* stub){ return sum += stub ? 1 : 0; });
      }
      stubs.reserve(nStubs);
      for (int layer = 0; layer < setup_->numLayers(); layer++)
        for (StubKFin* stub : inputStubs_[offset + layer])
          if (stub)
            stubs.push_back(stub);
    }
  }

  // read in and organize input tracks
  void KalmanFilter::consume(const StreamsTrack& tracks) {
    auto valid = [](int& sum, const FrameTrack& frame){ return sum += (frame.first.isNonnull() ? 1 : 0); };
    int nTracks(0);
    for (int channel = 0; channel < (int)inputTracks_.size(); channel++) {
      const StreamTrack& stream = tracks[region_ * inputTracks_.size() + channel];
      nTracks += accumulate(stream.begin(), stream.end(), 0, valid);
    }
    tracks_.reserve(nTracks);
    for (int channel = 0; channel < (int)inputTracks_.size(); channel++) {
      const vector<StubKFin*>& stubs = channelStubs_[channel];
      const StreamTrack& stream = tracks[region_ * inputTracks_.size() + channel];
      vector<TrackKFin*>& tracks = inputTracks_[channel];
      tracks.reserve(stream.size());
      for (const FrameTrack& frame : stream) {
        TrackKFin* track = nullptr;
        if (frame.first.isNonnull()) {
          tracks_.emplace_back(frame, dataFormats_, stubs);
          track = &tracks_.back();
        }
        tracks.push_back(track);
      }
    }
  }

  // fill output products
  void KalmanFilter::produce(TTDTC::Streams& acceptedStubs, StreamsTrack& acceptedTracks, TTDTC::Streams& lostStubs, StreamsTrack& lostTracks) {
    deque<State*> statesLost;
    deque<State*> statesAccepted;
    auto put = [this](const deque<State*>& states, TTDTC::Streams& streamsStubs, StreamsTrack& streamsTracks) {
      const int offset = region_ * setup_->numLayers();
      StreamTrack& tracks = streamsTracks[region_];
      tracks.reserve(states.size());
      for (int layer = 0; layer < setup_->numLayers(); layer++)
        streamsStubs[offset + layer].reserve(states.size());
      for (State* state : states) {
        tracks.emplace_back(state->frame());
        for (const StubKF& stub : state->stubs())
          streamsStubs[offset + stub.layer()].emplace_back(stub.frame());
        // adding a gap to all layer without a stub
        for (int layer : state->hitPattern().ids(false))
          streamsStubs[offset + layer].emplace_back(TTDTC::Frame());
      }
    };
    auto print = [this](const deque<State*>& stream) {
      vector<stringstream> sss(15);
      for (State* state : stream) {
        static constexpr int w = 8;
        if (!state) {
          for (stringstream& ss : sss)
            ss << setw(w) << " " << " ";
          continue;
        }
        sss[0] << setw(w) << x0_->integer(state->x0()) << " ";
        sss[1] << setw(w) << x1_->integer(state->x1()) << " ";
        sss[2] << setw(w) << x2_->integer(state->x2()) << " ";
        sss[3] << setw(w) << x3_->integer(state->x3()) << " ";
        sss[4] << setw(w) << H00_->integer(state->H00()) << " ";
        sss[5] << setw(w) << m0_->integer(state->m0()) << " ";
        sss[6] << setw(w) << m1_->integer(state->m1()) << " ";
        sss[7] << setw(w) << dataFormats_->format(Variable::dPhi, Process::kf).integer(state->dPhi()) << " ";
        sss[8] << setw(w) << dataFormats_->format(Variable::dZ, Process::kf).integer(state->dZ()) << " ";
        sss[9] << setw(w) << C00_->integer(state->C00()) << " ";
        sss[10] << setw(w) << C01_->integer(state->C01()) << " ";
        sss[11] << setw(w) << C11_->integer(state->C11()) << " ";
        sss[12] << setw(w) << C22_->integer(state->C22()) << " ";
        sss[13] << setw(w) << C23_->integer(state->C23()) << " ";
        sss[14] << setw(w) << C33_->integer(state->C33()) << " ";
      }
      for (stringstream& ss : sss)
        cout << ss.str() << endl;
    };
    vector<deque<State*>> streams(inputTracks_.size());
    //for (int channel = 0; channel < (int)inputTracks_.size(); channel++) {
    {int channel = 1;
      deque<State*>& stream = streams[channel];
      // proto state creation
      for (TrackKFin* track : inputTracks_[channel]) {
        State* state = nullptr;
        if (track) {
          states_.emplace_back(dataFormats_, track);
          state = &states_.back();
        }
        stream.push_back(state);
      }
      // state propagation
      for (layer_ = 0; layer_ < setup_->numLayers(); layer_++) {
        print(stream);
        layer(stream);
        cout << endl;
        print(stream);
        throw cms::Exception("...");
      }
      // untruncated best state selection
      deque<State*> untruncatedStream = stream;
      accumulator(untruncatedStream);
      // apply truncation
      if (enableTruncation_ && (int)stream.size() > setup_->numFrames())
        stream.resize(setup_->numFrames());
      // best state per candidate selection
      accumulator(stream);
      // storing of best states missed due to truncation
      sort(untruncatedStream.begin(), untruncatedStream.end());
      sort(stream.begin(), stream.end());
      set_difference(untruncatedStream.begin(), untruncatedStream.end(), stream.begin(), stream.end(), back_inserter(statesLost));
    }
    // stream merging
    for (const deque<State*>& stream : streams)
      copy(stream.begin(), stream.end(), back_inserter(statesAccepted));
    // apply truncation originated from merging
    if (enableTruncation_ && (int)statesAccepted.size() > setup_->numFrames()) {
      const auto it = next(statesAccepted.begin(), setup_->numFrames());
      copy(it, statesAccepted.end(), back_inserter(statesLost));
      statesAccepted.erase(it, statesAccepted.end());
    }
    put(statesAccepted, acceptedStubs, acceptedTracks);
    put(statesLost, lostStubs, lostTracks);
  }

  // adds a layer to states
  void KalmanFilter::layer(deque<State*>& stream) {
    static const int latency = 5;
    // dynamic state container for clock accurate emulation
    deque<State*> streamOutput;
    deque<State*> stack;
    // static delay container
    vector<State*> delay(latency, nullptr);
    delay.reserve(latency + 1);
    // each trip corresponds to a f/w clock tick
    // done if no states to process left, taking as much time as needed
    while (!stream.empty() || !stack.empty() || !all_of(delay.begin(), delay.end(), [](const State* state){ return state == nullptr; })) {
      State* state = pop_front(stream);
      if (!state)
        state = pop_front(stack);
      streamOutput.push_back(state);
      if (!state || !state->stub() || state->layer() != layer_)
        state = nullptr;
      if (state != nullptr)
        comb(state);
      delay.push_back(state);
      state = pop_front(delay);
      if (state != nullptr)
        stack.push_back(state);
    }
    stream = streamOutput;
    for (State*& state : stream) {
      if (!state || !state->stub() || state->layer() != layer_ )
        continue;
      update(state);
    }
  }

  // repicks combinatoric stubs for state
  void KalmanFilter::comb(State*& state) {
    // picks next stub on layer
    StubKFin* stub = state->stub();
    const int layer = stub->layer();
    TrackKFin* track = state->track();
    const vector<StubKFin*>& stubs = track->layerStubs(layer);
    const TTBV& hitPattern = state->hitPattern();
    const int pos = distance(stubs.begin(), find(stubs.begin(), stubs.end(), stub)) + 1;
    StubKFin* stubNext = nullptr;
    if (pos != (int)stubs.size())
      stubNext = stubs[pos];
    // picks next stub on different layer, nullifies state if skipping layer is not valid
    else {
      bool valid(true);
      // having already maximum number of added layers
      if (hitPattern.count() == setup_->kfMaxLayers())
        valid = false;
      // not enough layers remain after skipping
      if (hitPattern.count() + track->hitPattern().count(layer + 1, setup_->numLayers()) < setup_->kfMaxLayers())
        valid = false;
      if (valid) {
        // pick next stub on next populated layer
        for (int nextLayer = layer_ + 1; nextLayer < setup_->numLayers(); nextLayer++) {
          if (track->hitPattern(nextLayer)) {
            stubNext = track->layerStub(nextLayer);
            break;
          }
        }
      }
    }
    if (stubNext) {
      states_.emplace_back(state, stubNext);
      state = &states_.back();
    } else
      state = nullptr;
  }

  // best state selection
  void KalmanFilter::accumulator(deque<State*>& stream) {
    // accumulator delivers contigious stream of best state per track
    // remove gaps and not final states
    stream.erase(remove_if(stream.begin(), stream.end(), [this](State* state){ return !state || state->hitPattern().count() < setup_->kfMaxLayers(); }), stream.end());
    // update chi2
    for (State* state : stream)
      state->finish();
    // sort in number of skipped layers
    auto lessSkippedLayers = [](State* lhs, State* rhs) { return lhs->numSkippedLayers() < rhs->numSkippedLayers(); };
    stable_sort(stream.begin(), stream.end(), lessSkippedLayers);
    // sort in number of consistent stubs
    auto moreConsistentLayers = [](State* lhs, State* rhs) { return lhs->numConsistentLayers() > rhs->numConsistentLayers(); };
    stable_sort(stream.begin(), stream.end(), moreConsistentLayers);
    // sort in track id
    stable_sort(stream.begin(), stream.end(), [](State* lhs, State* rhs){ return lhs->trackId() < rhs->trackId(); });
    // keep first state (best due to previous sorts) per track id
    stream.erase(unique(stream.begin(), stream.end(), [](State* lhs, State* rhs){ return lhs->track() == rhs->track(); }), stream.end());
  }

  // updates state
  void KalmanFilter::update(State*& state) {
    // R-Z plane
    static const double dH = H00_->digi(setup_->chosenRofPhi() - setup_->chosenRofZ());
    const double H00 = H00_->digi(state->H00());
    //const double H12 = H12_->digi(H00 + dH);
    const double H12 = H12_->digi(state->H12());
    const double m1 = m1_->digi(state->m1());
    const double v1 = v1_->digi(state->v1());
    double x2 = x2_->digi(state->x2());
    double x3 = x3_->digi(state->x3());
    double C22 = C22_->digi(state->C22());
    double C23 = C23_->digi(state->C23());
    double C33 = C33_->digi(state->C33());
    v1_->updateRangeActual(v1);
    C22_->updateRangeActual(C22);
    C23_->updateRangeActual(C23);
    C33_->updateRangeActual(C33);
    H12_->updateRangeActual(H12);
    m1_->updateRangeActual(m1);
    const double r1C = x3_->digi(m1  - x3);
    const double r1 = r1_->digi(r1C - x2 * H12);
    const double S12 = S12_->digi(C23 + H12 * C22);
    const double S13 = S13_->digi(C33 + H12 * C23);
    const double R11C = S13_->digi(v1 + S13);
    const double R11 = R11_->digi(R11C + H12 * S12);
    // imrpoved dynamic cancelling
    const int msbZ = (int)ceil(log2(R11 / R11_->base()));
    const double R11Shifted = R11_->digi(R11 * pow(2., 16 - msbZ));
    const double R11Rough = R11Rough_->digi(R11Shifted);
    const double invR11Approx = invR11Approx_->digi(1. / R11Rough);
    const double invR11Cor = invR11Cor_->digi(2. - invR11Approx * R11Shifted);
    const double invR11Shifted = invR11Approx * invR11Cor;
    const double invR11 = invR11_->digi(invR11Shifted * pow(2., 16 - msbZ));
    //const double invR11 = invR11_->digi(1. / R11);
    const double K21 = K21_->digi(S12 * invR11);
    const double K31 = K31_->digi(S13 * invR11);
    x2 = x2_->digi(x2 + r1 * K21);
    x3 = x3_->digi(x3 + r1 * K31);
    C22 = C22_->digi(C22 - S12 * K21);
    C23 = C23_->digi(C23 - S13 * K21);
    C33 = C33_->digi(C33 - S13 * K31);
    // R-Phi plane
    const double m0 = m0_->digi(state->m0());
    const double v0 = v0_->digi(state->v0());
    double x0 = x0_->digi(state->x0());
    double x1 = x1_->digi(state->x1());
    double C00 = C00_->digi(state->C00());
    double C01 = C01_->digi(state->C01());
    double C11 = C11_->digi(state->C11());
    v0_->updateRangeActual(v0);
    C00_->updateRangeActual(C00);
    C01_->updateRangeActual(C01);
    C11_->updateRangeActual(C11);
    H00_->updateRangeActual(H00);
    m0_->updateRangeActual(m0);
    const double r0C = x1_->digi(m0 - x1);
    const double r0 = r0_->digi(r0C - x0 * H00);
    const double S00 = S00_->digi(C01  + H00 * C00);
    const double S01 = S01_->digi(C11  + H00 * C01);
    const double R00C = S01_->digi(v0 + S01);
    const double R00 = R00_->digi(R00C + H00 * S00);
    // improved dynamic cancelling
    const int msbPhi = (int)ceil(log2(R00 / R00_->base()));
    const double R00Shifted = R00_->digi(R00 * pow(2., 16 - msbPhi));
    const double R00Rough = R00Rough_->digi(R00Shifted);
    const double invR00Approx = invR00Approx_->digi(1. / R00Rough);
    const double invR00Cor = invR00Cor_->digi(2. - invR00Approx * R00Shifted);
    const double invR00Shifted = invR00Approx * invR00Cor;
    const double invR00 = invR00_->digi(invR00Shifted * pow(2., 16 - msbPhi));
    //const double invR00 = invR00_->digi(1. / R00);
    const double K00 = K00_->digi(S00 * invR00);
    const double K10 = K10_->digi(S01 * invR00);
    x0 = x0_->digi(x0 + r0 * K00);
    x1 = x1_->digi(x1 + r0 * K10);
    C00 = C00_->digi(C00 - S00 * K00);
    C01 = C01_->digi(C01 - S01 * K00);
    C11 = C11_->digi(C11 - S01 * K10);
    // residual cut
    bool valid(true);
    if (!dataFormats_->format(Variable::phi, Process::sf).inRange(r0) ||
        !dataFormats_->format(Variable::z, Process::sf).inRange(r1))
      valid = false;
    /*if (state->hitPattern().count() == 2) {
      if (abs(r1) > (setup_->psModule(state->stub()->ttStubRef()) ? 1. : 10.))
        valid = false;
      if (abs(r0) > 5.e-3)
        valid = false;
    }*/
    // parameter cut
    if (!dataFormats_->format(Variable::inv2R, Process::sf).inRange(x0) ||
        !dataFormats_->format(Variable::phiT, Process::sf).inRange(x1) ||
        !dataFormats_->format(Variable::cot, Process::sf).inRange(x2) ||
        !dataFormats_->format(Variable::zT, Process::sf).inRange(x3))
      valid = false;
    if (valid) {
      S12_->updateRangeActual(S12);
      S13_->updateRangeActual(S13);
      K21_->updateRangeActual(K21);
      K31_->updateRangeActual(K31);
      R11_->updateRangeActual(R11);
      R11Rough_->updateRangeActual(R11Rough);
      invR11Approx_->updateRangeActual(invR11Approx);
      invR11Cor_->updateRangeActual(invR11Cor);
      invR11_->updateRangeActual(invR11Shifted);
      r1_->updateRangeActual(r1);
      C22_->updateRangeActual(C22);
      C23_->updateRangeActual(C23);
      C33_->updateRangeActual(C33);
      x2_->updateRangeActual(x2);
      x3_->updateRangeActual(x3);
      S00_->updateRangeActual(S00);
      S01_->updateRangeActual(S01);
      K00_->updateRangeActual(K00);
      K10_->updateRangeActual(K10);
      R00_->updateRangeActual(R00);
      R00Rough_->updateRangeActual(R00Rough);
      invR00Approx_->updateRangeActual(invR00Approx);
      invR00Cor_->updateRangeActual(invR00Cor);
      invR00_->updateRangeActual(invR00Shifted);
      r0_->updateRangeActual(r0);
      C00_->updateRangeActual(C00);
      C01_->updateRangeActual(C01);
      C11_->updateRangeActual(C11);
      x0_->updateRangeActual(x0);
      x1_->updateRangeActual(x1);
      states_.emplace_back(State(state, (initializer_list<double>){x0, x1, x2, x3, C00, C11, C22, C33, C01, C23}));
      state = &states_.back();
    } else
      state = nullptr;
    static constexpr int w = 8;
    cout << "v0            " << setw(w) << v0_->integer(v0) << endl;
    cout << "v1            " << setw(w) << v1_->integer(v1) << endl;
    cout << "S00           " << setw(w) << S00_->integer(S00) << endl;
    cout << "S01           " << setw(w) << S01_->integer(S01) << endl;
    cout << "S12           " << setw(w) << S12_->integer(S12) << endl;
    cout << "S12           " << setw(w) << S12 << endl;
    cout << "S13           " << setw(w) << S13_->integer(S13) << endl;
    cout << "r0            " << setw(w) << r0_->integer(r0) << endl;
    cout << "r1            " << setw(w) << r1_->integer(r1) << endl;
    cout << "R00           " << setw(w) << R00_->integer(R00) << endl;
    cout << "R11           " << setw(w) << R11_->integer(R11) << endl;
    cout << "Shift0        " << setw(2) << 16 - msbPhi << endl;
    cout << "Shift1        " << setw(2) << 16 - msbZ << endl;
    cout << "R00Shifted    " << setw(w) << R00_->integer(R00Shifted) << endl;
    cout << "R11Shifted    " << setw(w) << R11_->integer(R11Shifted) << endl;
    cout << "R00Rough      " << setw(w) << R00Rough_->integer(R00Rough) << endl;
    cout << "R11Rough      " << setw(w) << R11Rough_->integer(R11Rough) << endl;
    cout << "invR00Approx  " << setw(w) << invR00Approx_->integer(invR00Approx) << endl;
    cout << "invR11Approx  " << setw(w) << invR11Approx_->integer(invR11Approx) << endl;
    cout << "invR00Cor     " << setw(w) << invR00Cor_->integer(invR00Cor) << endl;
    cout << "invR11Cor     " << setw(w) << invR11Cor_->integer(invR11Cor) << endl;
    cout << "invR00Shifted " << setw(w) << floor(invR00Shifted / invR00Approx_->base() / invR00Cor_->base()) << endl;
    cout << "invR11Shifted " << setw(w) << floor(invR11Shifted / invR11Approx_->base() / invR11Cor_->base()) << endl;
    cout << "invR00        " << setw(w) << invR00_->integer(invR00) << endl;
    cout << "invR11        " << setw(w) << invR11_->integer(invR11) << endl;
    cout << "K00           " << setw(w) << K00_->integer(K00) << endl;
    cout << "K10           " << setw(w) << K10_->integer(K10) << endl;
    cout << "K21           " << setw(w) << K21_->integer(K21) << endl;
    cout << "K31           " << setw(w) << K31_->integer(K21) << endl;
  }

  // remove and return first element of deque, returns nullptr if empty
  template<class T>
  T* KalmanFilter::pop_front(deque<T*>& ts) const {
    T* t = nullptr;
    if (!ts.empty()) {
      t = ts.front();
      ts.pop_front();
    }
    return t;
  }

  // remove and return first element of vector, returns nullptr if empty
  template<class T>
  T* KalmanFilter::pop_front(vector<T*>& ts) const {
    T* t = nullptr;
    if (!ts.empty()) {
      t = ts.front();
      ts.erase(ts.begin());
    }
    return t;
  }

} // namespace trackerTFP