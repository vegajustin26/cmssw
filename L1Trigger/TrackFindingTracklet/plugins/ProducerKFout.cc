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
#include <numeric>

using namespace std;
using namespace edm;
using namespace trackerDTC;
using namespace trackerTFP;

namespace trackFindingTracklet {

  /*! \class  trackFindingTracklet::ProducerKFout
   *  \brief  Converts KF output into TTTracks
   *  \author Thomas Schuh
   *  \date   2020, Oct
   */
  class ProducerKFout : public stream::EDProducer<> {
  public:
    explicit ProducerKFout(const ParameterSet&);
    ~ProducerKFout() override {}

  private:
    void beginRun(const Run&, const EventSetup&) override;
    void produce(Event&, const EventSetup&) override;
    void endJob() {}

    // ED input token of kf tracks
    EDGetTokenT<StreamsTrack> edGetToken_;
    // ED output token for TTTracks
    EDPutTokenT<TTTracks> edPutToken_;
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
    //
    const LayerEncoding* layerEncoding_;
    // used data formats
    const DataFormat* zT_;
    const DataFormat* cot_;
    const DataFormat* phi_;
    const DataFormat* z_;
  };

  ProducerKFout::ProducerKFout(const ParameterSet& iConfig) :
    iConfig_(iConfig)
  {
    const string& label = iConfig.getParameter<string>("LabelKF");
    const string& branch = iConfig.getParameter<string>("BranchAccepted");
    // book in- and output ED products
    edGetToken_ = consumes<StreamsTrack>(InputTag(label, branch));
    edPutToken_ = produces<TTTracks>(branch);
    // book ES products
    esGetTokenSetup_ = esConsumes<Setup, SetupRcd, Transition::BeginRun>();
    esGetTokenDataFormats_ = esConsumes<DataFormats, DataFormatsRcd, Transition::BeginRun>();
    esGetTokenLayerEncoding_ = esConsumes<LayerEncoding, LayerEncodingRcd, Transition::BeginRun>();
    // initial ES products
    setup_ = nullptr;
    dataFormats_ = nullptr;
    layerEncoding_ = nullptr;
    // used data formats
    zT_ = nullptr;
    cot_ = nullptr;
    phi_ = nullptr;
    z_ = nullptr;
  }

  void ProducerKFout::beginRun(const Run& iRun, const EventSetup& iSetup) {
    // helper class to store configurations
    setup_ = &iSetup.getData(esGetTokenSetup_);
    if (!setup_->configurationSupported())
      return;
    // check process history if desired
    if (iConfig_.getParameter<bool>("CheckHistory"))
      setup_->checkHistory(iRun.processHistory());
    // helper class to extract structured data from TTDTC::Frames
    dataFormats_ = &iSetup.getData(esGetTokenDataFormats_);
    //
    layerEncoding_ = &iSetup.getData(esGetTokenLayerEncoding_);
    // used data formats
    zT_ = &dataFormats_->format(Variable::zT, Process::kfin);
    cot_ = &dataFormats_->format(Variable::cot, Process::kfin);
    phi_ = &dataFormats_->format(Variable::phi, Process::sf);
    z_ = &dataFormats_->format(Variable::z, Process::sf);
  }

  void ProducerKFout::produce(Event& iEvent, const EventSetup& iSetup) {
    // empty KFout product
    TTTracks ttTracks;
    // read in KF Product and produce KFout product
    if (setup_->configurationSupported()) {
      Handle<StreamsTrack> handle;
      iEvent.getByToken<StreamsTrack>(edGetToken_, handle);
      const StreamsTrack& streams = *handle.product();
      int nTracks(0);
      for (const StreamTrack& stream : streams)
        nTracks += accumulate(stream.begin(), stream.end(), 0, [](int& sum, const FrameTrack& frame){ return sum += frame.first.isNonnull() ? 1 : 0; });
      ttTracks.reserve(nTracks);
      for (const StreamTrack& stream : streams) {
        for (const FrameTrack& frame : stream) {
          if (frame.first.isNull())
            continue;
          TrackKF track(frame, dataFormats_);
          const TTTrackRef& ttTrackRef = track.ttTrackRef();
          // get rz parameter
          double cot = ttTrackRef->tanL();
          double zT = ttTrackRef->z0() + setup_->chosenRofZ() * cot;
          const int binEta = track.sectorEta();
          cot -= setup_->sectorCot(binEta);
          zT -= setup_->sectorCot(binEta) * setup_->chosenRofZ();
          const int binZT = zT_->toUnsigned(zT_->integer(zT));
          const int binCot = cot_->toUnsigned(cot_->integer(cot));
          const vector<int>& layerEncoding = layerEncoding_->layerEncoding(binEta, binZT, binCot);
          vector<vector<TTStubRef>> layerStubs(setup_->numLayers());
          for (vector<TTStubRef>& stubs : layerStubs)
            stubs.reserve(setup_->kfMaxStubsPerLayer());
          vector<int> layerMap(setup_->numLayers(), 0);
          for (const TTStubRef& ttStubRef : track.frame().first->getStubRefs()) {
            const GlobalPoint& gp = setup_->stubPos(ttStubRef);
            const double phi = deltaPhi(gp.phi() - (ttTrackRef->phi() -  ttTrackRef->rInv() / 2. * gp.perp()));
            const double z = gp.z() - (ttTrackRef->z0() + ttTrackRef->tanL() * gp.perp());
            // cut on phi and z residuals
            if (!phi_->inRange(phi) || !z_->inRange(z))
              continue;
            const int layerId = setup_->layerId(ttStubRef);
            const int layerKF = distance(layerEncoding.begin(), find(layerEncoding.begin(), layerEncoding.end(), layerId));
            // cut on max 4 stubs per layer
            int& nLayerStubs = layerMap[layerKF];
            if (nLayerStubs == setup_->kfMaxStubsPerLayer())
              continue;
            layerStubs[layerKF].push_back(ttStubRef);
            nLayerStubs++;
          }
          vector<TTStubRef> ttstubRefs;
          ttstubRefs.reserve(track.hitPattern().count());
          for (int layer = 0; layer < setup_->numLayers(); layer++)
            if (track.hitPattern(layer))
              ttstubRefs.push_back(layerStubs[layer][track.layerMap(layer)]);
          track.ttStubRefs(ttstubRefs);
          ttTracks.emplace_back(track.ttTrack());
        }
      }
    }
    // store products
    iEvent.emplace(edPutToken_, move(ttTracks));
  }

} // namespace trackFindingTracklet

DEFINE_FWK_MODULE(trackFindingTracklet::ProducerKFout);