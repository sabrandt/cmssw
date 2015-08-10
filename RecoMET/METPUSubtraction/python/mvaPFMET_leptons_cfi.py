import FWCore.ParameterSet.Config as cms

# Single muon for Wjets
isomuons = cms.EDFilter(
    "MuonSelector",
    src = cms.InputTag('muons'),
    cut = cms.string(    "(isGlobalMuon) && abs(eta) < 2.4 && pt > 25"+#17. "+ cuts1&9
                         "&& isPFMuon"+#cuts5
                         "&& globalTrack.isNonnull"+#cuts5
                         "&& innerTrack.hitPattern.numberOfValidPixelHits > 0"+#cuts6
                         #"&& innerTrack.normalizedChi2 < 10"+#cuts7
                         "&& globalTrack.normalizedChi2 < 10"+#cuts8
                         "&& numberOfMatchedStations > 1"+
                         "&& innerTrack.hitPattern.numberOfValidTrackerHits>5"+#cuts6
                         "&& globalTrack.hitPattern.numberOfValidMuonHits>0"+#cuts6
                         #"&& (pfIsolationR04.sumChargedHadronPt+pfIsolationR04.sumNeutralHadronEt+pfIsolationR04.sumPhotonEt)/pt < 0.3"+#cuts4
                         "&&( (pfIsolationR04.sumChargedHadronPt+pfIsolationR04.sumNeutralHadronEt+pfIsolationR04.sumPhotonEt - 0.5*(pfIsolationR04.sumPUPt))/pt < 0.12"+
                         "&& (pfIsolationR04.sumChargedHadronPt)/pt < 0.12)"+ #cuts10
                         "&& abs(innerTrack().dxy)<0.2"#cuts2
                         #"&& abs(muonBestTrack().dz) < 0.5" #cuts3 #?????? removed
                         ),
    filter = cms.bool(False)
    )

isoelectrons = cms.EDFilter(
    "GsfElectronSelector",
    src = cms.InputTag('gedGsfElectrons'),
    cut = cms.string(
        "abs(eta) < 2.5 && pt > 25" +
        "&& ecalDrivenSeed" +
        #"&& gsfTrack.hitPattern.numberOfHits == 0"   + #comment by steph
        "&& (isolationVariables03.tkSumPt)/et              < 0.3"  +
        "&& (abs(eta) < 1.4442 " +
        "&& (pfIsolationVariables.sumChargedHadronPt+pfIsolationVariables.sumNeutralHadronEt+pfIsolationVariables.sumPhotonEt)/pt < 0.107587"  +
        "&&(pfIsolationVariables.sumChargedHadronPt)/pt < 0.107587"+
        "&& sigmaIetaIeta                                  < 0.009996 " +
        "&& abs(deltaPhiSuperClusterTrackAtVtx)            < 0.035973 " +
        "&& abs(deltaEtaSuperClusterTrackAtVtx)            < 0.008925 " +
        "&& hcalOverEcal                                   < 0.050537 " +
        "&& abs(1./superCluster.energy - 1./p)             < 0.091942 " +
        "&& abs(gsfTrack.dxy)                              < 0.012235) " +
        #endcap
        "|| (abs(eta)  > 1.566 "+
        "&& (pfIsolationVariables.sumChargedHadronPt+pfIsolationVariables.sumNeutralHadronEt+pfIsolationVariables.sumPhotonEt)/pt < 0.113254"  +
        "&&(pfIsolationVariables.sumChargedHadronPt)/pt < 0.113254"+
        "&& sigmaIetaIeta                                  < 0.030135 " +
        "&& abs(deltaPhiSuperClusterTrackAtVtx)            < 0.067879 " +
        "&& abs(deltaEtaSuperClusterTrackAtVtx)            < 0.007429 " +
        "&& hcalOverEcal                                   < 0.067778 " +
        "&& abs(1./superCluster.energy - 1./p)             < 0.098919 " +
        "&& abs(gsfTrack.dxy)                              < 0.027984) " 
        #"&& abs(gsfTrack.dz)                               < 0.133431)"
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
