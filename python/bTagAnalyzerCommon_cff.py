import FWCore.ParameterSet.Config as cms
from SimTracker.TrackHistory.TrackClassifier_cff import *
bTagAnalyzerCommon = cms.PSet(
    trackClassifier,
    runFatJets               = cms.bool(False),
    runSubJets               = cms.bool(False),
    allowJetSkipping         = cms.bool(True),
    useSelectedTracks        = cms.bool(True),
    useTrackHistory          = cms.bool(False),
    produceJetTrackTree      = cms.bool(False),
    produceJetTrackTruthTree = cms.bool(False),
    produceAllTrackTree      = cms.bool(False),
    producePtRelTemplate     = cms.bool(False),
    storeEventInfo           = cms.bool(True),
    storePatMuons            = cms.bool(False),
    storeTagVariables        = cms.bool(False),
    storeTagVariablesSubJets = cms.bool(False),
    storeCSVTagVariables     = cms.bool(False),
    storeCSVTagVariablesSubJets = cms.bool(False),
    storeCTagVariables       = cms.bool(False),
    doCTag                   = cms.bool(False),
    fillsvTagInfo            = cms.bool(False),
    fillPU                   = cms.bool(True),
    fillGenPruned            = cms.bool(True),
    fillQuarks               = cms.bool(False),
    selTagger                = cms.int32(2),
    MaxEta                   = cms.double(2.5),
    MinPt                    = cms.double(20.0),
    src                      = cms.InputTag('generator'),
    primaryVertexColl        = cms.InputTag('offlinePrimaryVertices'),
    tracksColl               = cms.InputTag('generalTracks'),
    BranchNamePrefix         = cms.string(''),
    Jets                     = cms.InputTag('selectedPatJets'),
    SubJets                  = cms.VInputTag(),
    SubJetLabels             = cms.vstring(),
    muonCollectionName       = cms.InputTag('muons'),
    patMuonCollectionName    = cms.InputTag('selectedPatMuons'),
    triggerTable             = cms.InputTag('TriggerResults'),
    genParticles             = cms.InputTag('genParticles'),
    prunedGenParticles       = cms.InputTag('prunedGenParticlesBoost'),
    candidates               = cms.InputTag("particleFlow"),
    maxDeltaR                = cms.double(0.4),
    explicitJTA              = cms.bool(False),
    subJetMaxDeltaR          = cms.double(0.4),
    subJetExplicitJTA        = cms.bool(True),
    use_ttbar_filter         = cms.bool(False),
    ttbarproducer            = cms.InputTag("ttbarselectionproducer"),
    clusterTPMap			 = cms.InputTag("tpClusterProducer"),
    rho                      = cms.InputTag("fixedGridRhoFastjetAll"),
    beta                     = cms.double(1.0),
    R0                       = cms.double(0.8),
    maxSVDeltaRToJet         = cms.double(0.7),
    weightFile               = cms.FileInPath('RecoBTag/PerformanceMeasurements/data/BoostedDoubleSV_AK8_BDT_v3.weights.xml.gz'),
    trackPairV0Filter        = cms.PSet(k0sMassWindow = cms.double(0.03)),
    distJetAxis              = cms.double(0.07),
    decayLength              = cms.double(5.0),
    deltaR                   = cms.double(0.3),
    distJetAxisSubJets       = cms.double(0.07),
    decayLengthSubJets       = cms.double(5.0),
    deltaRSubJets            = cms.double(0.3),
    TriggerPathNames = cms.vstring(
        # based on https://cmswbm.web.cern.ch/cmswbm/cmsdb/servlet/HLTSummary?RUN=284044&NAME=/cdaq/physics/Run2016/25ns15e33/v4.2.3/HLT/V2
        # PF Jets
        "HLT_PFJet40_v*",
        "HLT_PFJet60_v*",
        "HLT_PFJet80_v*",
        "HLT_PFJet140_v*",
        "HLT_PFJet200_v*",
        "HLT_PFJet260_v*",
        "HLT_PFJet320_v*",
        "HLT_PFJet400_v*",
        "HLT_PFJet450_v*",
        "HLT_PFJet500_v*",
        # Di PF Jets
        "HLT_DiPFJetAve40_v*",
        "HLT_DiPFJetAve60_v*",
        "HLT_DiPFJetAve80_v*",
        "HLT_DiPFJetAve140_v*",
        "HLT_DiPFJetAve200_v**",
        "HLT_DiPFJetAve260_v*",
        "HLT_DiPFJetAve320_v*",
        "HLT_DiPFJetAve400_v**",
        "HLT_DiPFJetAve500_v*",
        # AK 8 PF Jets
        #"HLT_AK8PFJet40_v*",
        #"HLT_AK8PFJet60_v*",
        #"HLT_AK8PFJet80_v*",
        #"HLT_AK8PFJet140_v*",
        #"HLT_AK8PFJet200_v*",
        #"HLT_AK8PFJet260_v*",
        #"HLT_AK8PFJet320_v*",
        #"HLT_AK8PFJet400_v*",
        #"HLT_AK8PFJet450_v*",
        #"HLT_AK8PFJet500_v*",
        # HT
        "HLT_HT200_v*",
        "HLT_HT275_v*",
        "HLT_HT325_v*",
        "HLT_HT425_v*",
        "HLT_HT575_v*",
        "HLT_HT650_v*",
        # HT range
        "HLT_HT430to450_v*",
        "HLT_HT450to470_v*",
        "HLT_HT470to500_v*",
        "HLT_HT500to550_v*",
        "HLT_HT550to650_v*",
        # PF HT
        #"HLT_PFHT125_v*",
        #"HLT_PFHT200_v*",
        #"HLT_PFHT250_v*",
        #"HLT_PFHT300_v*",
        #"HLT_PFHT350_v*",
        #"HLT_PFHT400_v*",
        #"HLT_PFHT475_v*",
        #"HLT_PFHT600_v*",
        #"HLT_PFHT650_v*",
        #"HLT_PFHT900_v*",
        # BTagMu
        "HLT_BTagMu_DiJet20_Mu5*",
        "HLT_BTagMu_DiJet70_Mu5*",
        "HLT_BTagMu_DiJet110_Mu5*",
        "HLT_BTagMu_DiJet110_L1FastJet_Mu5*",
        "HLT_BTagMu_DiJet170_Mu5_v*",
        "HLT_BTagMu_Jet300_L1FastJet_Mu5*",
        "HLT_BTagMu_AK8Jet300_Mu5_v*"
    ),
    TTbarTriggerPathNames = cms.vstring(
        # trigger for ttbar: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel#Triggers
        "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",  # 2-electron case
        "HLT_Mu17_Mu8_v*",                                                                          # 2-muon case1
        "HLT_Mu17_TkMu8_v*",                                                                        # 2-muon case2
        "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",                                      # muon + electron case1
        "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"                                       # muon + electron case2
        ),
    PFJet80TriggerPathNames = cms.vstring("HLT_PFJet80_v*")
)
