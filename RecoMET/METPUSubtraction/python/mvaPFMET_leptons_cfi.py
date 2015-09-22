import FWCore.ParameterSet.Config as cms

isomuons = cms.EDFilter(
    "MuonSelector",
    src = cms.InputTag('muons'),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    charged_hadron_iso = cms.InputTag(''),
    neutral_hadron_iso = cms.InputTag(''),
    photon_iso = cms.InputTag(''),
    typeID = cms.string("Tight"),
    ptCut = cms.double(20),
    etaCut = cms.double(2.4),
    typeIso = cms.string("dBeta"),
    relativeIsolationCut = cms.double(0.12),
    #find the MiniAOD cuts for muons
    cut = cms.string(    "(isGlobalMuon) && abs(eta) < 2.4 && pt > 20"+
                         "&& isPFMuon"+
                         "&& globalTrack.isNonnull"+#1
                         "&& innerTrack.hitPattern.numberOfValidPixelHits > 0"+#1
                         "&& globalTrack.normalizedChi2 < 10"+#1
                         "&& numberOfMatchedStations > 1"+#2
                         "&& innerTrack.hitPattern.numberOfValidTrackerHits>5"+#2
                         "&& globalTrack.hitPattern.numberOfValidMuonHits>0"+#2
                         "&&( (pfIsolationR04.sumChargedHadronPt+pfIsolationR04.sumNeutralHadronEt+pfIsolationR04.sumPhotonEt - 0.5*(pfIsolationR04.sumPUPt))/pt < 0.12"+#3
                         "&& (pfIsolationR04.sumChargedHadronPt)/pt < 0.12)"+ #3
                         "&& abs(innerTrack().dxy)<0.2"#3
                         )
    )
    
isoelectrons = cms.EDFilter(
    "GsfElectronSelector",
    src = cms.InputTag('gedGsfElectrons'),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    charged_hadron_iso = cms.InputTag(''),
    neutral_hadron_iso = cms.InputTag(''),
    photon_iso = cms.InputTag(''),
    typeID = cms.string("Tight"),
    #ptCut = cms.double(10),
    #etaCut = cms.double(2.5),
    typeIso = cms.string("dBeta"),
    relativeIsolationCut = cms.double(0.12),
    cut = cms.string("abs(superCluster.eta)<2.5 && pt > 20"+#1
                    "&& (abs(superCluster.eta) < 1.4442 " +#3
                    "&&(pfIsolationVariables.sumChargedHadronPt)/pt < 0.107587"+#12
                    #"&& sigmaIetaIeta                                  < 0.000 " + 
                    "&& abs(deltaPhiSuperClusterTrackAtVtx)            < 0.035973 " +
                    "&& abs(deltaEtaSuperClusterTrackAtVtx)            < 0.008925 " + 
                    "&& hcalOverEcal                                   < 0.050537 " + 
                    "&& abs(1./superCluster.energy - 1./p)             < 0.091942) " + 
                    #endcap
                    "|| (abs(superCluster.eta)  > 1.566 "+#3
                    "&&(pfIsolationVariables.sumChargedHadronPt)/pt < 0.113254"+ 
                    #"&& sigmaIetaIeta                                  < 0.00 " + 
                    "&& abs(deltaPhiSuperClusterTrackAtVtx)            < 0.067879 " +
                    "&& abs(deltaEtaSuperClusterTrackAtVtx)            < 0.007429 " +
                    "&& hcalOverEcal                                   < 0.086782  " + 
                    "&& abs(1./superCluster.energy - 1./p)             < 0.100683) "  
    
                    )
    #filter = cms.bool(False)
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
