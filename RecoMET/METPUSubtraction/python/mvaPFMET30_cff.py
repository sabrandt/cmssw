import FWCore.ParameterSet.Config as cms

from RecoMET.METProducers.PFMET_cfi import pfMet
from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from RecoMET.METPUSubtraction.mvaPFMET_leptons_cfi import *
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
## CV: importing mvaPFMET_leptons_cfi breaks produceAndDiscriminateHPSPFTaus sequence
##    (hpsPFTauDiscriminationByDecayModeFinding module gets overwritten by None,
##     in case RecoTauTag/Configuration/python/RecoPFTauTag_cff.py is loaded by
##     by top-level cfg.py file before RecoMET/METPUSubtraction/python/mvaPFMET_cff.py gets loaded)
##from RecoJets.JetProducers.PileupJetIDCutParams_cfi import full_53x_wp

calibratedAK4PFJetsForPFMVAMEt = cms.EDProducer('PFJetCorrectionProducer',
    src = cms.InputTag('ak4PFJets'),
    correctors = cms.vstring("ak4PFL1FastL2L3") # NOTE: use "ak5PFL1FastL2L3" for MC / "ak5PFL1FastL2L3Residual" for Data
)

from RecoJets.JetProducers.pileupjetidproducer_cfi import pileupJetIdEvaluator
from RecoJets.JetProducers.PileupJetIDParams_cfi import JetIdParams
puJetIdForPFMVAMEt = pileupJetIdEvaluator.clone(
    algos = cms.VPSet(
        cms.PSet(
        tmvaVariables = cms.vstring(
            "nvtx",
            "jetPt",
            "jetEta",
            "jetPhi",
            "dZ",
            "beta",
            "betaStar",
            "nCharged",
            "nNeutrals",
            "dR2Mean",
            "ptD",
            "frac01",
            "frac02",
            "frac03",
            "frac04",
            "frac05"
            ),
        tmvaWeights = cms.string("RecoJets/JetProducers/data/TMVAClassificationCategory_JetID_MET_53X_Dec2012.weights.xml.gz"),
        tmvaMethod = cms.string("JetID"),
        tmvaSpectators = cms.vstring(),
        JetIdParams = JetIdParams,                  
        impactParTkThreshold = cms.double(0.),
        version = cms.int32(-1),
        cutBased = cms.bool(False), 
        label = cms.string("full")
        )
        ),
    produceJetIds = cms.bool(True),
    runMvas = cms.bool(True),
    jets = cms.InputTag("calibratedAK4PFJetsForPFMVAMEt"),#calibratedAK4PFJetsForPFMVAMEt
    applyJec = cms.bool(True),
    inputIsCorrected = cms.bool(True),
    jec     = cms.string("AK4PF"),
)

pfMVAMEt = cms.EDProducer("PFMETProducerMVA",
    srcCorrJets = cms.InputTag('calibratedAK4PFJetsForPFMVAMEt'),
    srcUncorrJets = cms.InputTag('ak4PFJets30'),
    srcMVAPileupJetId = cms.InputTag('puJetIdForPFMVAMEt','fullDiscriminant'),
    srcPFCandidates = cms.InputTag('packedPFCandidates30'),
    srcVertices = cms.InputTag('offlinePrimaryVertices'),
    srcLeptons = cms.VInputTag("isomuons","isoelectrons"),#,"isotaus"), # NOTE: you need to set this to collections of electrons, muons and tau-jets
                                 #                                             passing the lepton reconstruction & identification criteria applied in your analysis
    minNumLeptons = cms.int32(0),                     
    srcRho = cms.InputTag('fixedGridRhoFastjetAll'),
    globalThreshold = cms.double(-1.),#pfMet.globalThreshold,
    minCorrJetPt = cms.double(-1.),
    inputFileNames = cms.PSet(
        U     = cms.FileInPath('RecoMET/METPUSubtraction/data/my_training/gbru_7_4_X_miniAOD_50NS_July2015.root'),
        DPhi  = cms.FileInPath('RecoMET/METPUSubtraction/data/my_training/gbrphi_7_4_X_miniAOD_50NS_July2015.root'),
        CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/my_training/gbru1cov_7_4_X_miniAOD_50NS_July2015.root'),
        CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/my_training/gbru2cov_7_4_X_miniAOD_50NS_July2015.root')
        #U     = cms.FileInPath('RecoMET/METPUSubtraction/data/Max250_Weights/RecoilCor_13TeV_Unity.root'),
        #DPhi  = cms.FileInPath('RecoMET/METPUSubtraction/data/Max250_Weights/PhiCor_13TeV_Unity.root'),
        #CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/Max250_Weights/CovU1_13TeV_Unity.root'),
        #CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/Max250_Weights/CovU2_13TeV_Unity.root')
        #U     = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmet_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        #DPhi  = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        #CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        #CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root')
    ),
    inputRecords = cms.PSet(
        U     = cms.string("RecoilCor"),
        DPhi  = cms.string("PhiCor"), 
        CovU1 = cms.string("CovU1"),
        CovU2 = cms.string("CovU2")
    ),
    loadMVAfromDB = cms.bool(False),                             

    corrector = cms.string("ak4PFL1Fastjet"),
    useType1  = cms.bool(True), 
    dZcut     = cms.double(0.1),
       
    verbosity = cms.int32(0)
)

packedPFCandidates30 = cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("particleFlow"),
    cut = cms.string("abs(eta) < 3.0"))

#process.AK4PFJetsPuppi30             = process.ak4PFJets.clone(src = cms.InputTag("puppi30"))
pfMet30                      = pfMet.clone();
pfMet30.src                  = cms.InputTag('packedPFCandidates30')
    #process.pfJetMETcorrPuppi30          = process.pfJetMETcorrPuppi.clone(src = 'AK4PFJetsPuppi30')
    #process.pfType1PuppiCorrectedMet30   = process.pfType1PuppiCorrectedMet.clone(src='pfMetPuppi30')
    #process.pfType1PuppiCorrectedMet30.srcType1Corrections = cms.VInputTag(cms.InputTag('pfJetMETcorrPuppi30', 'type1'))
ak4PFJets30                  = ak4PFJets.clone(src = cms.InputTag("packedPFCandidates30"))
#pfJetMETcorr30               = pfJetMETcorr.clone(src = 'ak4PFJets30')
#pfType1CorrectedMet30        = process.pfType1CorrectedMet.clone(src='pfMet30')
    #pfType1CorrectedMet30.srcType1Corrections = cms.VInputTag(cms.InputTag('pfJetMETcorr30', 'type1'))
    #allMET30           = cms.Sequence(process.pfMet30*process.pfMetPuppi30*process.pfJetMETcorr30*process.pfJetMETcorrPuppi30*process.pfType1CorrectedMet30*process.pfType1PuppiCorrectedMet30 )

    #MVA Met 3.0
calibratedAK4PFJetsForPFMVAMEt30 = calibratedAK4PFJetsForPFMVAMEt.clone( src = cms.InputTag('ak4PFJets30'))
puJetIdForPFMVAMEt30            = puJetIdForPFMVAMEt.clone(jets = cms.InputTag("calibratedAK4PFJetsForPFMVAMEt30"))
pfMVAMEt30                      = pfMVAMEt.clone()
pfMVAMEt30.srcMVAPileupJetId    = cms.InputTag('puJetIdForPFMVAMEt30','fullDiscriminant')
pfMVAMEt30.srcPFCandidates      = cms.InputTag("packedPFCandidates30")
pfMVAMEt30.srcUncorrJets        = "ak4PFJets30"
pfMVAMEt30.srcCorrJets          = "calibratedAK4PFJetsForPFMVAMEt30"

pfMVAMEt30Sequence  = cms.Sequence(
    (isomuonseq+isoelectronseq)*
    packedPFCandidates30*
    ak4PFJets30*
    calibratedAK4PFJetsForPFMVAMEt30*
    puJetIdForPFMVAMEt30*
    pfMVAMEt30
)
