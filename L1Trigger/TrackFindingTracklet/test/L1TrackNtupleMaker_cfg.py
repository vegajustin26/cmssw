############################################################
# define basic process
############################################################

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os
process = cms.Process("L1TrackNtuple")

############################################################
# edit options here
############################################################

# GEOMETRY = "D49"
GEOMETRY = "D76"
# Set L1 tracking algorithm:
# 'HYBRID' (baseline, 4par fit) or 'HYBRID_DISPLACED' (extended, 5par fit).
# (Or legacy algos 'TMTT' or 'TRACKLET').
L1TRKALGO = 'HYBRID_DISPLACED'

WRITE_DATA = False

############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.L1track = dict(limit = -1)
process.MessageLogger.Tracklet = dict(limit = -1)

if GEOMETRY == "D49":
    print("using geometry " + GEOMETRY + " (tilted)")
    process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
    process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
elif GEOMETRY == "D76":
    print("using geometry " + GEOMETRY + " (tilted)")
    process.load('Configuration.Geometry.GeometryExtended2026D76Reco_cff')
    process.load('Configuration.Geometry.GeometryExtended2026D76_cff')
else:
    print("this is not a valid geometry!!!")

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')


############################################################
# input and output
############################################################

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

#--- To use MCsamples scripts, defining functions get*data*(),
#--- follow instructions https://cernbox.cern.ch/index.php/s/enCnnfUZ4cpK7mT

#from MCsamples.Scripts.getCMSdata_cfi import *
#from MCsamples.Scripts.getCMSlocaldata_cfi import *

if GEOMETRY == "D49":
  # Read data from card files (defines getCMSdataFromCards()):
  #from MCsamples.RelVal_1120.PU200_TTbar_14TeV_cfi import *
  #inputMC = getCMSdataFromCards()

  # Or read .root files from directory on local computer:
  #dirName = "$myDir/whatever/"
  #inputMC=getCMSlocaldata(dirName)

  # Or read specified dataset (accesses CMS DB, so use this method only occasionally):
  #dataName="/RelValTTbar_14TeV/CMSSW_11_2_0_pre5-PU25ns_110X_mcRun4_realistic_v3_2026D49PU200-v1/GEN-SIM-DIGI-RAW"
  #inputMC=getCMSdata(dataName)

  # Or read specified .root file:
  #inputMC = ["/store/relval/CMSSW_11_3_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_113X_mcRun4_realistic_v3_2026D49PU200_rsb-v1/00000/00260a30-734a-4a3a-a4b0-f836ce5502c6.root"]

    # ttbar, PU=200, (100-500 events)
    # inputMC = ["/store/relval/CMSSW_11_3_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_113X_mcRun4_realistic_v3_2026D49PU200_rsb-v1/00000/00260a30-734a-4a3a-a4b0-f836ce5502c6.root",
    # "/store/relval/CMSSW_11_3_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_113X_mcRun4_realistic_v3_2026D49PU200_rsb-v1/00000/05943dfd-6493-450b-9ebd-b3344eb599e4.root",
    # "/store/relval/CMSSW_11_3_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_113X_mcRun4_realistic_v3_2026D49PU200_rsb-v1/00000/08d638b7-d2dc-4a8c-b0f7-30aecd83232c.root",
    # "/store/relval/CMSSW_11_3_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_113X_mcRun4_realistic_v3_2026D49PU200_rsb-v1/00000/0abc4f74-c1e7-4c89-b453-71cc9fea5f47.root",
    # "/store/relval/CMSSW_11_3_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_113X_mcRun4_realistic_v3_2026D49PU200_rsb-v1/00000/103cf801-a176-4d51-b9f8-e81420f59abf.root"]

    # single electron, 2-100, 100k events
    # inputMC = ["/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/0299f361-5e4c-4c25-a570-965bbdce8224.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/04aba858-3887-44ab-84c7-f7dd5b2b8fef.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/08c80f17-82bb-42ed-9c29-03be51b92cae.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/0a9db0be-9e9d-4576-890d-75bcef49271c.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/169e5654-bdf8-4a7f-a033-6e65e46b14a1.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/191ef11f-cc12-42d4-b538-62772258b72d.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/1aa5f87c-55ac-4e5f-bf06-136ebba8ae21.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/20abb02b-c8cf-454d-825f-a2388f94c3f3.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/27644dba-5104-43d2-94ab-f2042e9766f9.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/2af0c9a7-65d7-466d-88f0-325bb93a03cb.root",
    # ]

    # single muon, 2-100, 100k events
    # inputMC = ["/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/0210f667-5334-4592-84c8-097840e0991c.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/04d4d021-449d-4b25-9a38-9d7448ba8f00.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/073f02c2-cfde-487f-a5bf-29435025025c.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/0ad9414a-64ca-4aaf-897b-1691f5068d67.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/0b75d04e-8b00-42d0-a9b3-a6b410079756.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/0fff19a9-40c3-4c61-9413-126764dbba0f.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/175a4823-96ae-4ebb-9663-96ada5f13fbb.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/1868a4af-6b24-48ec-a5a6-0088d8f75bf4.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/18ec7825-7c15-4b33-a754-312d7862f170.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/10000/1ba20f6b-b806-4540-a950-39cf9c2ab9fc.root"]

    # displaced muon, pt-range 2-100
    inputMC = ["/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/00000/0203c9f6-76c3-4ab1-bb66-fd0a2a164906.root",
    "/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/00000/032ae701-c3fe-40c5-b827-ff4443ec4e8c.root",
    "/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/00000/078d458c-8b89-4054-9af4-40cca40e73b5.root",
    "/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/00000/0872d9f5-4cae-4e43-be04-12b64929b1fc.root",
    "/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/00000/0976f5f1-d5c8-47bf-91a3-137721a9c048.root",
    "/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/00000/0f5c3da8-1bce-46f9-b286-488e11d27914.root",
    "/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/00000/114f2976-26f5-40e5-8772-a23438dd8c04.root",
    "/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/00000/11f0f85d-4073-4211-874a-4083b1feac1a.root",
    "/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/00000/133d0761-24c2-433a-8c79-4a0418e1a3ea.root",
    "/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D49noPU-v1/00000/1408932b-9b3e-4cfc-b4c6-e73f75f1a22b.root"]



elif GEOMETRY == "D76":
    # inputMC = ["/store/relval/CMSSW_11_3_0_pre6/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_113X_mcRun4_realistic_v6_2026D76PU200-v1/00000/00026541-6200-4eed-b6f8-d3a1fd720e9c.root"]

    # ttbar, PU=200, (100-500 events)
    # inputMC = ["/store/relval/CMSSW_11_3_0_pre6/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_113X_mcRun4_realistic_v6_2026D76PU200-v1/00000/00026541-6200-4eed-b6f8-d3a1fd720e9c.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_113X_mcRun4_realistic_v6_2026D76PU200-v1/00000/013d0125-8f6e-496b-8335-614398c9210d.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_113X_mcRun4_realistic_v6_2026D76PU200-v1/00000/058bd134-86de-47e1-bcde-379ed9b79e1b.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_113X_mcRun4_realistic_v6_2026D76PU200-v1/00000/0915d66c-cbd4-4ef6-9971-7dd59e198b56.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_113X_mcRun4_realistic_v6_2026D76PU200-v1/00000/09823c8d-e443-4066-8347-8c704929cb2b.root"]

    # single electron, 2-100, 100k events
    # inputMC = ["/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v2/10000/00315da1-a37c-42ac-b927-5bf795a55be3.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v2/10000/0062d3f8-2357-4431-a4b5-c2b130dd4cbd.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v2/10000/01e0064c-4741-4096-8b4c-0213f1b5b769.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v2/10000/0301021d-ed9c-4828-a282-d7e5fd694daf.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v2/10000/04c29edf-6cd3-431c-bce5-7beae2c07e25.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v2/10000/071f5296-6f7b-4e4f-9740-e464767d64c7.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v2/10000/07be45ea-8af2-4bd8-be2a-746d6a13d3e5.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v2/10000/08b19213-1d38-4ba7-8002-2e5a8c415c6b.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v2/10000/09fa5be1-de1f-42d3-a651-608875bb008a.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v2/10000/0c43b034-39c9-4c75-856a-7e6187ecb7a5.root",
    # ]
    #
    # #
    # # single muon, 2-100, 100k events
    # inputMC = ["/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/10000/05f802b7-b0b3-4cca-8b70-754682c3bb4c.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/10000/0b69ed0a-66e9-403a-88f0-fb3115615461.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/10000/0f4dea68-7574-43bb-97c3-5382d68a2704.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/10000/156b3ca6-c74a-4f46-ae5e-03d9b01acd4c.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/10000/16727f1d-2922-4e0a-8239-82e1ffecd43b.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/10000/1af620bf-1f6d-4d5a-8170-4135ac798581.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/10000/1dc513d9-75fc-44c0-b8e0-e2925323416b.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/10000/2010f402-2133-4c3a-851b-1ae68fe23eb3.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/10000/228dfbba-3d5c-42b9-b827-cfa8f11a2f38.root",
    # "/store/relval/CMSSW_11_3_0_pre6/RelValSingleMuFlatPt2To100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/10000/27d006b1-d023-4775-8430-382e6962149c.root",
    # ]

    # displaced muon, pt-range 2-100, 100,000 events
    inputMC = ["/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/011da61a-9524-4a96-b91f-03e8690af3bd.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/01f4d0af-9c8f-4454-aa4d-01fbcdeb0927.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/02ddce4d-7e8a-4bfc-8c92-87d07e734830.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/063029da-53c4-49a0-bfde-a003f1b80d67.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/08bedc23-a5b0-4378-8fc3-75146bf6aeb0.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/08dfc14e-839f-4c1d-af31-ded2d611ee50.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/0a81a697-6c4d-4dc5-ba45-81244d727c0a.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/0b32b93f-72bb-4b0c-8aff-2f3657f05556.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/0d04a4cd-5c72-4fd2-bde6-0e4452a6724e.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/0d38a7c2-d06c-4cd3-8127-94c19fb386d7.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/1423771d-373b-4509-a94e-82c83e3cb163.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/146c3bd7-e5e6-4cc4-8c3d-7133cf271d13.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/1eb4c4d9-449b-4d89-a82c-2a69db5a78bc.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/1f4afb4c-9dce-4d6d-8c96-efbfe78a3b67.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/21c4e5a5-18e6-46bb-8c88-98052f468c54.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/2569e74e-fa6d-4601-9428-9cc38a0f5ee7.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/2c96526e-796d-49f8-8649-ca1b784dc77e.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/2dc285d6-8bf1-4373-a433-8398f499c8fa.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/3225ed3f-599e-4d4a-a790-bc4f2b75c2e8.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/32290383-26c0-4827-9d84-ddf9b0788297.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/33a93827-488a-4760-aab9-80ea7ae5e6e7.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/34d4d931-01d7-423c-abca-7ff60dbd0877.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/3538233d-fc3a-46e4-89b6-3edcdfff1613.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/35c1422d-4f08-4e20-8eff-8fdb80110e3b.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/3669b229-8cc1-4e6d-8f6c-06909f7d8ade.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/367b140b-5b5f-47ae-8049-5703f1b70cd4.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/368e04f5-7c73-43e0-92f7-63ec4bd1de04.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/374cbac3-f31b-462e-8b6d-10f5560da7f6.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/3835320f-b4cf-4a4f-9bd9-428b8aa8e8d5.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/38ae643a-aa42-47a6-8c80-f76660e3b685.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/397c81a1-57b5-414a-8eb5-0a5b6b303f04.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/3a4325db-90ed-41e3-be0b-b2171714bff8.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/422e3331-0e59-4e93-8048-eeedb7da7e2d.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/45c69eb4-de5e-4b4c-ad91-5dc6ce03e3fd.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/48691fde-12f5-4dc6-9be8-080152fa18a2.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/486bb4a5-4496-4a67-bc3c-4bfec515b5c2.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/49acd988-c608-42e6-877a-6ceb2a89ae55.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/4aa13315-f749-484b-9b15-fa6823d46aee.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/4ba7c863-2d75-48cf-8763-3d5d373f9102.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/527869fa-0757-4375-abab-6d20c89807a0.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/5459a821-4afa-4418-80fe-e2b501c0898e.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/555b0ea4-fcfa-4e1c-88d5-68030ba8c579.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/56bd9bf1-51ce-4e33-b5bc-377ae348ac73.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/5741e2cb-b93b-475a-bce3-5d692bb2fad6.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/5acc8d62-95c9-498d-beea-6f1d3c024e69.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/5ca138fb-6d01-4bef-ab89-88cfd04b7a61.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/5dbe8104-868e-438a-b4aa-e7dfc7bc2e21.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/61b3bcd6-f05f-4e82-8c4e-a627d4d91107.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/639a564d-f7b9-4116-9f6d-bdc18197267b.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/6a012e3b-0507-4f42-82d7-09080786a076.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/6a505f52-67e1-4c87-b8c2-342692a5ed38.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/6d225f51-3626-4a91-8e6b-211848fe29b3.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/7ab8bb10-734f-44c9-b2d1-e6abfb3c607b.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/81c8b9ca-27ce-4fd7-a8ef-f404ed9d96ec.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/827d8361-c965-48e5-9d73-46bf31f348e6.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/82bd5997-e420-4986-a22e-21cb8adaee31.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/881d9c6b-4749-4255-90a4-977dbf1d0c2e.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/8cb09f0e-da92-4443-9c5e-3b2418312769.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/8fa43c21-fdaa-46d1-b964-50af426b6c5a.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/92f4efca-506c-4acb-a677-5927b07e087d.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/96f5d762-2f9c-42ab-9734-08183b3c1fbf.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/9843c74f-903c-4d1b-bf28-e22abffd688e.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/99471f00-8e80-4548-bd22-f1c5cd5ddcce.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/9a71ab73-afd1-43b3-9649-05a30ea0d5ea.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/9a81008e-1f40-408f-b693-389f967485d6.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/9ab93b66-566a-4c52-9a8a-65a0009d3536.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/9d6518f7-8ea8-4c93-8c8d-01b66fd9598f.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/a18a95b1-33c3-4bef-8fc0-265923984c07.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/a1bdfa18-8cd9-400c-ad25-9c81116d1a66.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/a1fbb6b4-6171-49b9-a88c-830cb2dd97f3.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/a27982b0-cfe5-40fe-9580-bff3535e647c.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/a6f8eb00-ac13-4dd3-af44-16d59acfbe0d.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/a7304d38-cf4f-4db4-b58e-4d69c76181ad.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/a8ceccba-6f59-4b92-a1b2-b55cf4ec8fc5.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/a958122a-d57d-458b-97c1-e87dd9bc591e.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/aa7dda72-8d59-469f-b904-fe6effc5419d.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/ad83edf5-ac1b-43ed-985d-849e608d8cad.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/b1973881-1585-48e3-ba36-5b749912a460.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/b948c794-7d18-4206-9ff7-d6a9d1a98338.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/b9b6f4ff-090a-4a56-9102-06e3222c0075.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/bcaefb95-d4e1-48cc-a303-068008ab6d10.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/bd8402f3-f787-49a2-9371-e2b787f65c67.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/bf1306c0-aa9b-420c-9554-8e8e364746f8.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/c50b1521-7073-4d0c-9323-a52baea4a12f.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/c541720d-b278-40b4-adeb-f5cb22e4c7b5.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/c8ebf167-566c-48a7-950c-f9bb66849fae.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/d12be334-3d17-4706-a288-1834b9350c21.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/d163f0b7-57f5-4a47-b42b-856a30baa2a9.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/d26f4954-bb32-4630-8e8b-860c8a32d961.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/e22f14e7-35db-4dce-8370-e03b243770b8.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/e27fd868-1e64-4ade-81ef-9dbc3e8837dc.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/e461a74b-fc95-4cf1-a4b4-b086b8f29bb3.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/e67a2f5d-172b-45c5-a816-53da920868f2.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/eba4d307-884c-47e3-94fb-709f94d276c6.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/eebe6b3d-53b8-47c9-8ee7-1f27f17f2168.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/f2058ea7-a744-450f-9e05-8ef67e60003f.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/f3984d44-f158-4186-a743-dd65eeec1a27.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/f4fcf7e1-e9f7-45bd-a97f-43a4997a7c9f.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/fb546594-9ab4-4659-95f0-ec7644a2cd85.root",
"/store/relval/CMSSW_11_3_0_pre6/RelValDisplacedMuPt2To100Dxy100/GEN-SIM-DIGI-RAW/113X_mcRun4_realistic_v6_2026D76noPU-v1/00000/fd41cb1d-6164-41d0-9254-6fdff3b5e952.root"]

else:
  print("this is not a valid geometry!!!")

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(*inputMC))

process.TFileService = cms.Service("TFileService", fileName = cms.string('disp_muon_PU0_'+GEOMETRY+'.root'), closeFileFast = cms.untracked.bool(True))
# process.TFileService = cms.Service("TFileService", fileName = cms.string('muon_PU0_'+GEOMETRY+'.root'), closeFileFast = cms.untracked.bool(True))
process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True))


############################################################
# L1 tracking: stubs / DTC emulation
############################################################

process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')

# remake stubs?
#from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *
#process.load("SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff")

#from SimTracker.TrackTriggerAssociation.TTClusterAssociation_cfi import *
#TTClusterAssociatorFromPixelDigis.digiSimLinks = cms.InputTag("simSiPixelDigis","Tracker")

#process.TTClusterStub = cms.Path(process.TrackTriggerClustersStubs)
#process.TTClusterStubTruth = cms.Path(process.TrackTriggerAssociatorClustersStubs)


# DTC emulation
process.load('L1Trigger.TrackerDTC.ProducerES_cff')
process.load('L1Trigger.TrackerDTC.ProducerED_cff')

# load code that analyzes DTCStubs
#process.load('L1Trigger.TrackerDTC.Analyzer_cff')

# modify default cuts
#process.TrackTriggerSetup.FrontEnd.BendCut = 5.0
#process.TrackTriggerSetup.Hybrid.MinPt = 1.0

process.dtc = cms.Path(process.TrackerDTCProducer)#*process.TrackerDTCAnalyzer)


############################################################
# L1 tracking
############################################################

process.load("L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff")

# HYBRID: prompt tracking
if (L1TRKALGO == 'HYBRID'):
    process.TTTracksEmulation = cms.Path(process.L1HybridTracks)
    process.TTTracksEmulationWithTruth = cms.Path(process.L1HybridTracksWithAssociators)
    NHELIXPAR = 4
    L1TRK_NAME  = "TTTracksFromTrackletEmulation"
    L1TRK_LABEL = "Level1TTTracks"
    L1TRUTH_NAME = "TTTrackAssociatorFromPixelDigis"

# HYBRID: extended tracking
elif (L1TRKALGO == 'HYBRID_DISPLACED'):
    process.TTTracksEmulation = cms.Path(process.L1ExtendedHybridTracks)
    process.TTTracksEmulationWithTruth = cms.Path(process.L1ExtendedHybridTracksWithAssociators)
    NHELIXPAR = 5
    L1TRK_NAME  = "TTTracksFromExtendedTrackletEmulation"
    L1TRK_LABEL = "Level1TTTracks"
    L1TRUTH_NAME = "TTTrackAssociatorFromPixelDigisExtended"

# LEGACY ALGORITHM (EXPERTS ONLY): TRACKLET
elif (L1TRKALGO == 'TRACKLET'):
    print("\n WARNING: This is not the baseline algorithm! Prefer HYBRID or HYBRID_DISPLACED!")
    print("\n To run the Tracklet-only algorithm, ensure you have commented out 'CXXFLAGS=-DUSEHYBRID' in BuildFile.xml & recompiled! \n")
    process.TTTracksEmulation = cms.Path(process.L1HybridTracks)
    process.TTTracksEmulationWithTruth = cms.Path(process.L1HybridTracksWithAssociators)
    NHELIXPAR = 4
    L1TRK_NAME  = "TTTracksFromTrackletEmulation"
    L1TRK_LABEL = "Level1TTTracks"
    L1TRUTH_NAME = "TTTrackAssociatorFromPixelDigis"

# LEGACY ALGORITHM (EXPERTS ONLY): TMTT
elif (L1TRKALGO == 'TMTT'):
    print("\n WARNING: This is not the baseline algorithm! Prefer HYBRID or HYBRID_DISPLACED! \n")
    process.load("L1Trigger.TrackFindingTMTT.TMTrackProducer_Ultimate_cff")
    L1TRK_PROC  =  process.TMTrackProducer
    L1TRK_NAME  = "TMTrackProducer"
    L1TRK_LABEL = "TML1TracksKF4ParamsComb"
    L1TRUTH_NAME = "TTTrackAssociatorFromPixelDigis"
    NHELIXPAR = 4
    L1TRK_PROC.EnableMCtruth = cms.bool(False) # Reduce CPU use by disabling internal histos.
    L1TRK_PROC.EnableHistos  = cms.bool(False)
    process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
    process.load("SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff")
    process.TTTrackAssociatorFromPixelDigis.TTTracks = cms.VInputTag( cms.InputTag(L1TRK_NAME, L1TRK_LABEL) )
    process.TTTracksEmulation = cms.Path(process.offlineBeamSpot*L1TRK_PROC)
    process.TTTracksEmulationWithTruth = cms.Path(process.offlineBeamSpot*L1TRK_PROC*process.TrackTriggerAssociatorTracks)

else:
    print("ERROR: Unknown L1TRKALGO option")
    exit(1)


############################################################
# Define the track ntuple process, MyProcess is the (unsigned) PDGID corresponding to the process which is run
# e.g. single electron/positron = 11
#      single pion+/pion- = 211
#      single muon+/muon- = 13
#      pions in jets = 6
#      taus = 15
#      all TPs = 1
############################################################

process.L1TrackNtuple = cms.EDAnalyzer('L1TrackNtupleMaker',
                                       MyProcess = cms.int32(1),
                                       DebugMode = cms.bool(False),      # printout lots of debug statements
                                       SaveAllTracks = cms.bool(True),   # save *all* L1 tracks, not just truth matched to primary particle
                                       SaveStubs = cms.bool(False),      # save some info for *all* stubs
                                       L1Tk_nPar = cms.int32(NHELIXPAR), # use 4 or 5-parameter L1 tracking?
                                       L1Tk_minNStub = cms.int32(4),     # L1 tracks with >= 4 stubs
                                       TP_minNStub = cms.int32(4),       # require TP to have >= X number of stubs associated with it
                                       TP_minNStubLayer = cms.int32(4),  # require TP to have stubs in >= X layers/disks
                                       TP_minPt = cms.double(2.0),       # only save TPs with pt > X GeV
                                       TP_maxEta = cms.double(2.5),      # only save TPs with |eta| < X
                                       TP_maxZ0 = cms.double(30.0),      # only save TPs with |z0| < X cm
                                       L1TrackInputTag = cms.InputTag(L1TRK_NAME, L1TRK_LABEL),         # TTTrack input
                                       MCTruthTrackInputTag = cms.InputTag(L1TRUTH_NAME, L1TRK_LABEL),  # MCTruth input
                                       # other input collections
                                       L1StubInputTag = cms.InputTag("TTStubsFromPhase2TrackerDigis","StubAccepted"),
                                       MCTruthClusterInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"),
                                       MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
                                       TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                       TrackingVertexInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                       # tracking in jets (--> requires AK4 genjet collection present!)
                                       TrackingInJets = cms.bool(False),
                                       GenJetInputTag = cms.InputTag("ak4GenJets", "")
                                       )

process.ana = cms.Path(process.L1TrackNtuple)


############################################################
# final schedule of what is to be run
############################################################

# use this if you want to re-run the stub making
# process.schedule = cms.Schedule(process.TTClusterStub,process.TTClusterStubTruth,process.dtc,process.TTTracksEmulationWithTruth,process.ana)

# use this if cluster/stub associators not available
# process.schedule = cms.Schedule(process.TTClusterStubTruth,process.dtc,process.TTTracksEmulationWithTruth,process.ana)

# use this to only run tracking + track associator
process.schedule = cms.Schedule(process.dtc,process.TTTracksEmulationWithTruth,process.ana)


############################################################
# write output dataset?
############################################################

if (WRITE_DATA):
  process.writeDataset = cms.OutputModule("PoolOutputModule",
      splitLevel = cms.untracked.int32(0),
      eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
      outputCommands = process.RAWSIMEventContent.outputCommands,
      fileName = cms.untracked.string('output_dataset.root'), ## ADAPT IT ##
      dataset = cms.untracked.PSet(
          filterName = cms.untracked.string(''),
          dataTier = cms.untracked.string('GEN-SIM')
      )
  )
  process.writeDataset.outputCommands.append('keep  *TTTrack*_*_*_*')
  process.writeDataset.outputCommands.append('keep  *TTStub*_*_*_*')

  process.pd = cms.EndPath(process.writeDataset)
  process.schedule.append(process.pd)
