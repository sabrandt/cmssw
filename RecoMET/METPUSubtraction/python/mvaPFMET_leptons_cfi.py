import FWCore.ParameterSet.Config as cms

# Single muon for Wjets
isomuons = cms.EDFilter(
    "MuonSelector",
    src = cms.InputTag('muons'),
    cut = cms.string(
        "isGlobalMuon && isTrackerMuon && abs(eta) < 2.5 && pt > 9.5 " +
        "&& globalTrack.normalizedChi2 < 10 "+
        "&& globalTrack.hitPattern.numberofValidMuonHits > 0 "+
        "&& numberOfMatchedStations > 1 "+
        "&& abs(muonBestTrack.dxy(vertex.position)) < 0.2 "+
        "&& abs(muonBestTrack.dz(vertex.position)) < 0.5 "+
        "&& innerTrack.hitPattern.numberOfValidPixelHits > 0 "+
        "&& innerTrack.trackerLayersWithMeasurement > 5 "+
        "&& (pfIsolationR04.sumChargedHadronPt+pfIsolationR04.sumNeutralHadronEt+pfIsolationR04.sumPhotonEt)/pt < 0.3"
        ),
    filter = cms.bool(False)     
    )

# compared to w/z selection:
# missing conversion veto,
# isolation requirement not the same
# missing hits requirement

isoelectrons = cms.EDFilter(
    "GsfElectronSelector",
    src = cms.InputTag('gedGsfElectrons'),
    cut = cms.string(
        "abs(eta) < 2.5 && pt > 9.5 && ecalDrivenSeed" +
        "&& (isolationVariables03.tkSumPt)/et < 0.3 " +
        #barrel
        "&& ((abs(eta) < 1.4442 " +
        "&& sigmaIetaIeta                                  < 0.009996 " +
        "&& abs(deltaPhiSuperClusterTrackAtVtx)            < 0.035973 " +
        "&& abs(deltaEtaSuperClusterTrackAtVtx)            < 0.008925 " +
        "&& hcalOverEcal                                   < 0.050537 " +
        "&& abs(1./superCluster.energy - 1./p)             < 0.091942 " +
        "&& abs(dxy(vertex.position))                      < 0.012235 " + 
        "&& abs(dz(vertex.position))                       < 0.042020) " + 
        #endcap
        "|| (abs(eta)  > 1.566 "+
        "&& sigmaIetaIeta                                  < 0.030135 " +
        "&& abs(deltaPhiSuperClusterTrackAtVtx)            < 0.067879 " +
        "&& abs(deltaEtaSuperClusterTrackAtVtx)            < 0.007429 " +
        "&& hcalOverEcal                                   < 0.067778 " +
        "&& abs(1./superCluster.energy - 1./p)             < 0.098919 " +
        "&& abs(dxy(vertex.position))                      < 0.027984 " + 
        "&& abs(dz(vertex.position))                       < 0.133431))"
        ),
    filter = cms.bool(False)
    )

from RecoJets.Configuration.RecoPFJets_cff import kt6PFJets as dummy
kt6PFJetsForRhoComputationVoronoiMet = dummy.clone(
        doRhoFastjet = True,
        voronoiRfact = 0.9
        )

from RecoTauTag.RecoTau.PFRecoTauDiscriminationByHPSSelection_cfi import hpsSelectionDiscriminator
hpsPFTauDiscriminationByDecayModeFinding = hpsSelectionDiscriminator.clone(
        PFTauProducer = cms.InputTag('hpsPFTauProducer')
            )

from RecoTauTag.RecoTau.TauDiscriminatorTools import requireLeadTrack
# Define decay mode prediscriminant
requireDecayMode = cms.PSet(
        BooleanOperator = cms.string("and"),
            decayMode = cms.PSet(
            Producer = cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding'),
                    cut = cms.double(0.5)
                )
        )

from RecoTauTag.Configuration.HPSPFTaus_cff import hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits

#hpsPFTauDiscriminationAgainstMuon2 = cms.EDProducer("PFRecoTauDiscriminationAgainstMuon2",
                                                         #PFTauProducer = cms.InputTag('hpsPFTauProducer'),
                                                         #Prediscriminants = requireDecayMode.clone(),
                                                         #discriminatorOption = cms.string('loose'), # available options are: 'loose', 'medium', 'tight'
                                                         #HoPMin = cms.double(0.2)
                                                     #)


#hpsPFTauDiscriminationByMVAIsolation = cms.EDProducer(
    #"PFRecoTauDiscriminationByMVAIsolation",
            #PFTauProducer = cms.InputTag('hpsPFTauProducer'),
            #rhoProducer = cms.InputTag('kt6PFJetsForRhoComputationVoronoiMet','rho'),
            #Prediscriminants = requireDecayMode.clone(),
            #gbrfFilePath = cms.FileInPath('RecoTauTag/RecoTau/data/gbrfTauIso_v2.root'),
            #returnMVA = cms.bool(False),
            #mvaMin = cms.double(0.8),
            #)

isotaus = cms.EDFilter(
    "PFTauSelector",
    src = cms.InputTag('hpsPFTauProducer'),
    BooleanOperator = cms.string("and"),
    discriminators = cms.VPSet(
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),       selectionCut=cms.double(0.5)),
    #cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByMVAIsolation"),           selectionCut=cms.double(0.5)),
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits"),           selectionCut=cms.double(0.5)),
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseElectronRejection"), selectionCut=cms.double(0.5)),
    #cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationAgainstMuon2"),             selectionCut=cms.double(0.5)) 
    ),
    cut = cms.string("abs(eta) < 2.3 && pt > 19.0 "),
    filter = cms.bool(False)
    )

isomuonseq     = cms.Sequence(isomuons)
isoelectronseq = cms.Sequence(isoelectrons)
isotauseq      = cms.Sequence(
    hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits*
     #kt6PFJetsForRhoComputationVoronoiMet*
     #hpsPFTauDiscriminationByMVAIsolation*
     #hpsPFTauDiscriminationAgainstMuon2*
     isotaus
    )

leptonSelection = cms.PSet(
    SelectEvents = cms.PSet(
    SelectEvents = cms.vstring(
    'isomuonseq',
    'isoelectronseq',
    'isotauseq')
    )
    )
