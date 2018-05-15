from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.Heppy.physicsobjects.Electron import Electron
from PhysicsTools.Heppy.physicsobjects.Muon import Muon
from PhysicsTools.Heppy.physicsobjects.Photon import Photon
#from CMGTools.TTHAnalysis.tools.EfficiencyCorrector import EfficiencyCorrector

from PhysicsTools.HeppyCore.utils.deltar import bestMatch
#from CMGTools.TTHAnalysis.electronCalibrator import ElectronCalibrator
import PhysicsTools.HeppyCore.framework.config as cfg
from PhysicsTools.HeppyCore.utils.deltar import * 
from PhysicsTools.Heppy.physicsutils.genutils import *


import math, os
from ROOT import heppy, TLorentzVector
import math
import copy


class LeptonCleaner( Analyzer ):

    def __init__(self, cfg_ana, cfg_comp, looperName ):
        super(LeptonCleaner,self).__init__(cfg_ana,cfg_comp,looperName)

        self.lepPtMin = getattr(self.cfg_ana, 'minLepPt', -1)
        self.lepSelCut = getattr(self.cfg_ana, 'lepSelCut', lambda lep : True)
        self.eleGammaDR =  getattr(self.cfg_ana, 'eleGammaDR', 1.0)
        self.muGammaDR =  getattr(self.cfg_ana, 'muGammaDR', 0.5)

        if not hasattr(self.cfg_ana ,"collectionPostFix"):self.cfg_ana.collectionPostFix=""


    #----------------------------------------
    # DECLARATION OF HANDLES OF LEPTONS STUFF   
    #----------------------------------------
        

    # def declareHandles(self):
    #     super(LeptonCleaner, self).declareHandles()

    #     #leptons
    #     self.handles['muons'] = AutoHandle(self.cfg_ana.muons,"std::vector<pat::Muon>")            
    #     self.handles['electrons'] = AutoHandle(self.cfg_ana.electrons,"std::vector<pat::Electron>")            

    def beginLoop(self, setup):
        super(LeptonCleaner,self).beginLoop(setup)
        self.counters.addCounter('events')
        count = self.counters.counter('events')
        count.register('all events')
    

    # def makeAllMuons(self, event):
    #     """
    #            make a list of all muons, and apply basic corrections to them
    #     """
    #     # Start from all muons
    #     allmuons = map( Muon, self.handles['muons'].product() )

    #     # mt2 gg specific: clean photons of muons in dr of 0.5
    #     photonCleanedMu = []
    #     photons = []
    #     if hasattr(event, 'selectedPhotons'):
    #         photons = [ g for g in event.selectedPhotons ] 

    #     if len(photons)>0:
    #         for mu in enumerate(allmuons):
    #             isNotPhoton = True
    #             for gamma in photons:
    #                 dr = deltaR(gamma.eta(),gamma.phi(), mu.eta(),mu.phi())
    #                 if (dr < 0.5): 
    #                     isNotPhoton = False
    #             if isNotPhoton:
    #                 photonCleanedMu.append(mu)
    #         allmuons = photonCleanedMu

    #         return allmuons

    # def makeAllElectrons(self, event):
    #     """
    #            make a list of all electrons, and apply basic corrections to them
    #     """
    #     allelectrons = map( Electron, self.handles['electrons'].product() )

    #     # mt2 gg specific: clean photons of electrons in dr of 0.5
    #     photonCleanedEle = []
    #     photons = []
    #     if hasattr(event, 'selectedPhotons'):
    #         photons = [ g for g in event.selectedPhotons ] 

    #     print "Doing the electron removal if its near a photon"
    #     print len(photons)

    #     if len(photons)>0:
    #         print "Doing the electron removal if its near a photon"
    #         print len(photons)
    #         for ele in enumerate(allelectrons):
    #             isNotPhoton = True
    #             for gamma in photons:
    #                 dr = deltaR(gamma.eta(),gamma.phi(), ele.eta(),ele.phi())
    #                 if ( dr < 1. ): 
    #                     isNotPhoton = False
    #             if isNotPhoton:
    #                 photonCleanedEle.append(ele)
    #         allelectrons = photonCleanedEle

    #         return allelectron

    # #------------------
    # # MAKE LEPTON LISTS
    # #------------------

    
    # def makeLeptons(self, event):
    #     ### inclusive leptons = all leptons that could be considered somewhere in the analysis, with minimal requirements (used e.g. to match to MC)
    #     event.inclusiveLeptons = []
    #     ### selected leptons = subset of inclusive leptons passing some basic id definition and pt requirement
    #     ### other    leptons = subset of inclusive leptons failing some basic id definition and pt requirement
    #     event.selectedLeptons = []
    #     event.selectedMuons = []
    #     event.selectedElectrons = []
    #     event.otherLeptons = []

    #     #muons
    #     allmuons = self.makeAllMuons(event)

    #     #electrons        
    #     allelectrons = self.makeAllElectrons(event)


    def process(self, event):
        self.readCollections( event.input )
        self.counters.counter('events').inc('all events')

        photonCleanedLeptons = []
        photons = []
        if hasattr(event, 'selectedPhotons'):
            photons = [ g for g in event.selectedPhotons ] 

        leptons = []
        if hasattr(event, 'selectedLeptons'):
            leptons = [ l for l in event.selectedLeptons if l.pt() > self.lepPtMin and self.lepSelCut(l) ]

        if len(photons)>0:
            for lep in leptons:
                isNotPhoton = True
                for gamma in photons:
                    dr = deltaR(gamma.eta(),gamma.phi(), lep.eta(),lep.phi())
                    if ( dr < 1.0 and abs(lep.pdgId())==11 ): #electron
                        isNotPhoton = False
                    if ( dr < 0.5 and abs(lep.pdgId())==13):  #muon
                        isNotPhoton = False
                if isNotPhoton:
                    photonCleanedLeptons.append(lep)
    
        #self.selectedLeptons  = [ l for l in event.selectedLeptons if l not in self.discardedLeptons ]
        self.selectedLeptons  = [ l for l in photonCleanedLeptons ]
#        self.selectedLeptons = photonCleanedLeptons
#        self.selectedLeptons = photonCleanedLeptons

# #        #call the leptons functions
# #        self.makeLeptons(event)

        setattr(event,"selectedLeptons"          +self.cfg_ana.collectionPostFix, self.selectedLeptons          ) 

            
        return True



#A default config
setattr(LeptonCleaner,"defaultConfig",cfg.Analyzer(
#    verbose=False,
    class_object=LeptonCleaner,
    # input collections

    minLepPt = 10,
    lepSelCut = lambda lep : True,

    eleGammaDR = 1.0,
    muGammaDR = 0.5,

    collectionPostFix = ""

#     rhoMuon= 'fixedGridRhoFastjetAll',
#     rhoElectron = 'fixedGridRhoFastjetAll',
# ##    photons='slimmedPhotons',
#     # energy scale corrections and ghost muon suppression (off by default)
#     # inclusive very loose muon selection
#     inclusive_muon_id  = "POG_ID_Loose",
#     inclusive_muon_pt  = 3,
#     inclusive_muon_eta = 2.4,
#     inclusive_muon_dxy = 0.5,
#     inclusive_muon_dz  = 1.0,
#     muon_dxydz_track   = "muonBestTrack",
#     # loose muon selection
#     loose_muon_id     = "POG_ID_Loose",
#     loose_muon_pt     = 5,
#     loose_muon_eta    = 2.4,
#     loose_muon_dxy    = 0.05,
#     loose_muon_dz     = 0.2,
#     loose_muon_relIso = 0.4,
#     # loose_muon_isoCut = lambda muon :muon.miniRelIso < 0.2 
#     # inclusive very loose electron selection
#     inclusive_electron_id  = "",
#     inclusive_electron_pt  = 5,
#     inclusive_electron_eta = 2.5,
#     inclusive_electron_dxy = 0.5,
#     inclusive_electron_dz  = 1.0,
#     inclusive_electron_lostHits = 1.0,
#     # loose electron selection
#     loose_electron_id     = "", #POG_MVA_ID_NonTrig_full5x5",
#     loose_electron_pt     = 7,
#     loose_electron_eta    = 2.4,
#     loose_electron_dxy    = 0.05,
#     loose_electron_dz     = 0.2,
#     loose_electron_relIso = 0.4,
#     # loose_electron_isoCut = lambda electron : electron.miniRelIso < 0.1
#     loose_electron_lostHits = 1.0,
#     # muon isolation correction method (can be "rhoArea" or "deltaBeta")
    )
)
