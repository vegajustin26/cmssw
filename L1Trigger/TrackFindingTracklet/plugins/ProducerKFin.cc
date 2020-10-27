#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"

#include "L1Trigger/TrackerDTC/interface/Setup.h"
#include "L1Trigger/TrackerTFP/interface/DataFormats.h"
#include "L1Trigger/TrackerTFP/interface/LayerEncoding.h"

#include <string>
#include <vector>
#include <deque>
#include <iterator>
#include <cmath>
#include <numeric>

using namespace std;
using namespace edm;
using namespace trackerDTC;
using namespace trackerTFP;

namespace trackFindingTracklet {

  /*! \class  trackFindingTracklet::ProducerKFin
   *  \brief  transforms hybrid TTTracks into KF input
   *  \author Thomas Schuh
   *  \date   2020, Oct
   */
  class ProducerKFin : public stream::EDProducer<> {
  public:
    explicit ProducerKFin(const ParameterSet&);
    ~ProducerKFin() override {}

  private:
    virtual void beginRun(const Run&, const EventSetup&) override;
    virtual void produce(Event&, const EventSetup&) override;
    virtual void endJob() {}

    // ED input token of TTTracks
    EDGetTokenT<TTTracks> edGetTokenTTTracks_;
    // ED output token for stubs
    EDPutTokenT<TTDTC::Streams> edPutTokenAcceptedStubs_;
    EDPutTokenT<TTDTC::Streams> edPutTokenLostStubs_;
    // ED output token for tracks
    EDPutTokenT<StreamsTrack> edPutTokenAcceptedTracks_;
    EDPutTokenT<StreamsTrack> edPutTokenLostTracks_;
    // Setup token
    ESGetToken<Setup, SetupRcd> esGetTokenSetup_;
    // DataFormats token
    ESGetToken<DataFormats, DataFormatsRcd> esGetTokenDataFormats_;
    // LayerEncoding token
    ESGetToken<LayerEncoding, LayerEncodingRcd> esGetTokenLayerEncoding_;
    // configuration
    ParameterSet iConfig_;
    // helper class to store configurations
    const Setup* setup_;
    // helper class to extract structured data from TTDTC::Frames
    const DataFormats* dataFormats_;
    // helper class to encode layer
    const LayerEncoding* layerEncoding_;
    //
    bool enableTruncation_;
  };

  ProducerKFin::ProducerKFin(const ParameterSet& iConfig) :
    iConfig_(iConfig)
  {
    const InputTag& inputTag = iConfig.getParameter<InputTag>("InputTag");
    const string& branchAccepted = iConfig.getParameter<string>("BranchAccepted");
    const string& branchLost = iConfig.getParameter<string>("BranchLost");
    // book in- and output ED products
    edGetTokenTTTracks_ = consumes<TTTracks>(inputTag);
    edPutTokenAcceptedStubs_ = produces<TTDTC::Streams>(branchAccepted);
    edPutTokenAcceptedTracks_ = produces<StreamsTrack>(branchAccepted);
    edPutTokenLostStubs_ = produces<TTDTC::Streams>(branchLost);
    edPutTokenLostTracks_ = produces<StreamsTrack>(branchLost);
    // book ES products
    esGetTokenSetup_ = esConsumes<Setup, SetupRcd, Transition::BeginRun>();
    esGetTokenDataFormats_ = esConsumes<DataFormats, DataFormatsRcd, Transition::BeginRun>();
    esGetTokenLayerEncoding_ = esConsumes<LayerEncoding, LayerEncodingRcd, Transition::BeginRun>();
    // initial ES products
    setup_ = nullptr;
    dataFormats_ = nullptr;
    layerEncoding_ = nullptr;
    //
    enableTruncation_ = iConfig.getParameter<bool>("EnableTruncation");
  }

  void ProducerKFin::beginRun(const Run& iRun, const EventSetup& iSetup) {
    // helper class to store configurations
    setup_ = &iSetup.getData(esGetTokenSetup_);
    if (!setup_->configurationSupported())
      return;
    // check process history if desired
    if (iConfig_.getParameter<bool>("CheckHistory"))
      setup_->checkHistory(iRun.processHistory());
    // helper class to extract structured data from TTDTC::Frames
    dataFormats_ = &iSetup.getData(esGetTokenDataFormats_);
    // helper class to encode layer
    layerEncoding_ = &iSetup.getData(esGetTokenLayerEncoding_);
  }

  void ProducerKFin::produce(Event& iEvent, const EventSetup& iSetup) {
    auto toFrameStub = [](StubKFin* stub) { return stub->frame(); };
    auto toFrameTrack = [](const TrackKFin& track){ return track.frame(); };
    const DataFormat& dfCot = dataFormats_->format(Variable::cot, Process::sf);
    const DataFormat& dfZT = dataFormats_->format(Variable::zT, Process::sf);
    const DataFormat& dfQoverPt = dataFormats_->format(Variable::qOverPt, Process::sf);
    const DataFormat& dfPhiT = dataFormats_->format(Variable::phiT, Process::sf);
    const DataFormat& dfPhi = dataFormats_->format(Variable::phi, Process::sf);
    const DataFormat& dfZ = dataFormats_->format(Variable::z, Process::sf);
    // empty KFin products
    TTDTC::Streams streamAcceptedStubs(dataFormats_->numStreams(Process::kf) * setup_->numLayers());
    StreamsTrack streamAcceptedTracks(dataFormats_->numStreams(Process::kf));
    TTDTC::Streams streamLostStubs(dataFormats_->numStreams(Process::kf) * setup_->numLayers());
    StreamsTrack streamLostTracks(dataFormats_->numStreams(Process::kf));
    // read in hybrid track finding product and produce KFin product
    if (setup_->configurationSupported()) {
      Handle<TTTracks> handleTTTracks;
      iEvent.getByToken<TTTracks>(edGetTokenTTTracks_, handleTTTracks);
      const TTTracks& ttTracks = *handleTTTracks.product();
      vector<TTTrackRefs> ttTrackRefsRegions(setup_->numRegions());
      vector<int> nTTTracksRegions(setup_->numRegions(), 0);
      for (const TTTrack<Ref_Phase2TrackerDigi_>& ttTrack : ttTracks)
        nTTTracksRegions[ttTrack.phiSector()]++;
      for (int region = 0; region < setup_->numRegions(); region++)
        ttTrackRefsRegions[region].reserve(nTTTracksRegions[region]);
      int i(0);
      for (const TTTrack<Ref_Phase2TrackerDigi_>& ttTrack : ttTracks)
        ttTrackRefsRegions[ttTrack.phiSector()].emplace_back(TTTrackRef(handleTTTracks, i++));
      for (int region = 0; region < setup_->numRegions(); region++) {
        const TTTrackRefs& ttTrackRefs = ttTrackRefsRegions[region];
        vector<TrackKFin> tracks;
        tracks.reserve(ttTrackRefs.size());
        vector<StubKFin> stubs;
        const int nStubs = accumulate(ttTrackRefs.begin(), ttTrackRefs.end(), 0, [](int& sum, const TTTrackRef& ttTrackRef){ return sum += ttTrackRef->getStubRefs().size(); });
        stubs.reserve(nStubs);
        int trackId(0);
        vector<int> numLayerStubs(setup_->numLayers(), 0);
        for (const TTTrackRef& ttTrackRef : ttTrackRefs) {
          // cut on more then 256 tracks
          if (trackId >= dataFormats_->format(Variable::trackId, Process::kfin).range())
            continue;
          // get rz parameter
          double cot = ttTrackRef->tanL();
          double zT = ttTrackRef->z0() + setup_->chosenRofZ() * cot;
          int binEta(-1);
          for (; binEta < setup_->numSectorsEta(); binEta++) {
            if (zT < sinh(setup_->boundarieEta(binEta + 1)) * setup_->chosenRofZ())
              break;
          }
          // cut on outer eta sector boundaries
          if (binEta == -1 || binEta == setup_->numSectorsEta())
            continue;
          cot -= setup_->sectorCot(binEta);
          zT -= setup_->sectorCot(binEta) * setup_->chosenRofZ();
          // cut on eta and |z0| < 15 cm
          if (!dfZT.inRange(zT) || !dfCot.inRange(cot))
            continue;
          const int binZT = dfZT.toUnsigned(dfZT.integer(zT));
          const int binCot = dfCot.toUnsigned(dfCot.integer(cot));
          const vector<int>& layerEncoding = layerEncoding_->layerEncoding(binEta, binZT, binCot);
          // get rphi parameter
          double qOverPt = ttTrackRef->rInv() / 2.; // times / over +-2?
          double phiT = deltaPhi(ttTrackRef->phi() - setup_->hybridChosenRofPhi() * qOverPt - region * setup_->baseRegion());
          const int sectorPhi = phiT < 0. ? 0 : 1; // dirty hack
          phiT -= (sectorPhi - .5) * setup_->baseSector();
          // cut on nonant size and pt > 2 GeV
          if (!dfPhiT.inRange(phiT) || !dfQoverPt.inRange(qOverPt))
            continue;
          // loop over stubs
          TTBV hitPattern(0, setup_->numLayers());
          vector<int> layerMap(setup_->numLayers(), 0);
          for (const TTStubRef& ttStubRef : ttTrackRef->getStubRefs()) {
            const GlobalPoint& gp = setup_->stubPos(ttStubRef);
            const double r = gp.perp() - setup_->hybridChosenRofPhi();
            const double phi = deltaPhi(gp.phi() - (ttTrackRef->phi() - qOverPt * gp.perp()));
            const double z = gp.z() - (ttTrackRef->z0() + ttTrackRef->tanL() * gp.perp());
            const double layer = distance(layerEncoding.begin(), find(layerEncoding.begin(), layerEncoding.end(), setup_->layerId(ttStubRef)));
            // cut on phi and z residuals
            if (!dfPhi.inRange(phi) || !dfZ.inRange(z))
              continue;
            hitPattern.set(layer);
            int& nLayerStubs = layerMap[layer];
            // cut on max 4 stubs per layer
            if (nLayerStubs == setup_->kfMaxStubsPerLayer())
              continue;
            nLayerStubs++;
            numLayerStubs[layer]++;
            stubs.emplace_back(ttStubRef, dataFormats_, r, phi, z, trackId, layer);
          }
          tracks.emplace_back(ttTrackRef, dataFormats_, hitPattern, setup_->layerMap(hitPattern, layerMap), phiT, qOverPt, zT, cot, sectorPhi, binEta, trackId++);
        }
        // truncate and transform stubs
        vector<vector<StubKFin*>> layerStubs(setup_->numLayers());
        for (int layer = 0; layer < setup_->numLayers(); layer++)
          layerStubs[layer].reserve(numLayerStubs[layer]);
        for (StubKFin& stub : stubs)
          layerStubs[stub.layer()].push_back(&stub);
        const int offset = region * setup_->numLayers();
        for (int layer = 0; layer < setup_->numLayers(); layer++) {
          TTDTC::Stream& acceptedStubs = streamAcceptedStubs[offset + layer];
          TTDTC::Stream& lostStubs = streamLostStubs[offset + layer];
          vector<StubKFin*>& stubsLayer = layerStubs[layer];
          if (enableTruncation_ && (int)stubsLayer.size() > setup_->numFrames()) {
            const auto it = next(stubsLayer.begin(), setup_->numFrames());
            lostStubs.reserve(setup_->numFrames() - (int)stubsLayer.size());
            transform(it, stubsLayer.end(), back_inserter(lostStubs), toFrameStub);
            stubsLayer.erase(it, stubsLayer.end());
          }
          transform(stubsLayer.begin(), stubsLayer.end(), back_inserter(acceptedStubs), toFrameStub);
        }
        // truncate and transform tracks
        StreamTrack& acceptedTracks = streamAcceptedTracks[region];
        StreamTrack& lostTracks = streamLostTracks[region];
        if (enableTruncation_ && (int)tracks.size() > setup_->numFrames()) {
          const auto it = next(tracks.begin(), setup_->numFrames());
          lostTracks.reserve(setup_->numFrames() - (int)tracks.size());
          transform(it, tracks.end(), back_inserter(lostTracks), toFrameTrack);
          tracks.erase(it, tracks.end());
        }
        transform(tracks.begin(), tracks.end(), back_inserter(acceptedTracks), toFrameTrack);
      }
    }
    // store products
    iEvent.emplace(edPutTokenAcceptedStubs_, move(streamAcceptedStubs));
    iEvent.emplace(edPutTokenAcceptedTracks_, move(streamAcceptedTracks));
    iEvent.emplace(edPutTokenLostStubs_, move(streamLostStubs));
    iEvent.emplace(edPutTokenLostTracks_, move(streamLostTracks));
  }

} // namespace trackFindingTracklet

DEFINE_FWK_MODULE(trackFindingTracklet::ProducerKFin);