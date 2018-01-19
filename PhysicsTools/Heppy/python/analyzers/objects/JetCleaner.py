import math, os
from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.Heppy.physicsobjects.PhysicsObjects import Jet
from PhysicsTools.HeppyCore.utils.deltar import deltaR2, deltaPhi, matchObjectCollection, matchObjectCollection2, bestMatch,matchObjectCollection3
from PhysicsTools.Heppy.physicsutils.JetReCalibrator import JetReCalibrator
import PhysicsTools.HeppyCore.framework.config as cfg

from PhysicsTools.Heppy.physicsutils.QGLikelihoodCalculator import QGLikelihoodCalculator

import copy
def cleanNearestJetOnly(jets,leptons,deltaR):
    dr2 = deltaR**2
    good = [ True for j in jets ]
    for l in leptons:
        ibest, d2m = -1, dr2
        for i,j in enumerate(jets):
            d2i = deltaR2(l.eta(),l.phi(), j.eta(),j.phi())
            if d2i < d2m:
                ibest, d2m = i, d2i
        if ibest != -1: good[ibest] = False
    return [ j for (i,j) in enumerate(jets) if good[i] == True ] 

def cleanJetsAndLeptons(jets,leptons,deltaR,arbitration):
    dr2 = deltaR**2
    goodjet = [ True for j in jets ]
    goodlep = [ True for l in leptons ]
    for il, l in enumerate(leptons):
        ibest, d2m = -1, dr2
        for i,j in enumerate(jets):
            d2i = deltaR2(l.eta(),l.phi(), j.eta(),j.phi())
            if d2i < dr2:
                choice = arbitration(j,l)
                if choice == j:
                   # if the two match, and we prefer the jet, then drop the lepton and be done
                   goodlep[il] = False
                   break 
                elif choice == (j,l) or choice == (l,j):
                   # asked to keep both, so we don't consider this match
                   continue
            if d2i < d2m:
                ibest, d2m = i, d2i
        # this lepton has been killed by a jet, then we clean the jet that best matches it
        if not goodlep[il]: continue 
        if ibest != -1: goodjet[ibest] = False
    return ( [ j for (i ,j) in enumerate(jets)    if goodjet[i ] == True ], 
             [ l for (il,l) in enumerate(leptons) if goodlep[il] == True ] )



class JetCleaner( Analyzer ):
    """Taken from RootTools.JetAnalyzer, simplified, modified, only cleaning    """
    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(JetCleaner,self).__init__(cfg_ana, cfg_comp, looperName)
        self.doPuId = getattr(self.cfg_ana, 'doPuId', True)
        self.jetLepDR = getattr(self.cfg_ana, 'jetLepDR', 0.4)
        self.jetLepArbitration = getattr(self.cfg_ana, 'jetLepArbitration', lambda jet,lepton: lepton) 
        self.lepPtMin = getattr(self.cfg_ana, 'minLepPt', -1)
        self.lepSelCut = getattr(self.cfg_ana, 'lepSelCut', lambda lep : True)
        self.jetGammaDR =  getattr(self.cfg_ana, 'jetGammaDR', 0.4)
        self.jetGammaLepDR =  getattr(self.cfg_ana, 'jetGammaLepDR', 0.4)
        self.cleanFromLepAndGammaSimultaneously = getattr(self.cfg_ana, 'cleanFromLepAndGammaSimultaneously', False)
        if self.cleanFromLepAndGammaSimultaneously:
            if hasattr(self.cfg_ana, 'jetGammaLepDR'):
                self.jetGammaLepDR =  self.jetGammaLepDR 
            elif (self.jetGammaDR == self.jetLepDR):
                self.jetGammaLepDR = self.jetGammaDR
            else:
                raise RuntimeError, "DR for simultaneous cleaning of jets from leptons and photons is not defined, and dR(gamma, jet)!=dR(lep, jet)"
        if not hasattr(self.cfg_ana ,"collectionPostFix"):self.cfg_ana.collectionPostFix=""
    
    def beginLoop(self, setup):
        super(JetCleaner,self).beginLoop(setup)

    def process(self, event):
        leptons = []
        if hasattr(event, 'selectedLeptons'):
            leptons = [ l for l in event.selectedLeptons if l.pt() > self.lepPtMin and self.lepSelCut(l) ]
        if self.cfg_ana.cleanJetsFromTaus and hasattr(event, 'selectedTaus'):
            leptons = leptons[:] + event.selectedTaus
        if self.cfg_ana.cleanJetsFromIsoTracks and hasattr(event, 'selectedIsoCleanTrack'):
            leptons = leptons[:] + event.selectedIsoCleanTrack

        jetsEtaCut = [j for j in event.jets if abs(j.eta()) <  self.cfg_ana.jetEta ]
        self.cleanJetsAll, cleanLeptons = cleanJetsAndLeptons(jetsEtaCut, leptons, self.jetLepDR, self.jetLepArbitration)

        self.cleanJets    = [j for j in self.cleanJetsAll if abs(j.eta()) <  self.cfg_ana.jetEtaCentral ]
        self.cleanJetsFwd = [j for j in self.cleanJetsAll if abs(j.eta()) >= self.cfg_ana.jetEtaCentral ]
        self.discardedJets = [j for j in event.jets if j not in self.cleanJetsAll]
        if hasattr(event, 'selectedLeptons') and self.cfg_ana.cleanSelectedLeptons:
            self.discardedLeptons = [ l for l in leptons if l not in cleanLeptons ]
            self.selectedLeptons  = [ l for l in event.selectedLeptons if l not in self.discardedLeptons ]
        for lep in leptons:
            if hasattr(lep, "jetOverlap"):
                if lep.jetOverlap in self.cleanJetsAll:
                    #print "overlap reco", lep.p4().pt(), lep.p4().eta(), lep.p4().phi(), lep.jetOverlap.p4().pt(), lep.jetOverlap.p4().eta(), lep.jetOverlap.p4().phi()
                    lep.jetOverlapIdx = self.cleanJetsAll.index(lep.jetOverlap)
                elif lep.jetOverlap in self.discardedJets:
                    #print "overlap discarded", lep.p4().pt(), lep.p4().eta(), lep.p4().phi(), lep.jetOverlap.p4().pt(), lep.jetOverlap.p4().eta(), lep.jetOverlap.p4().phi()
                    lep.jetOverlapIdx = 1000 + self.discardedJets.index(lep.jetOverlap)

        ## First cleaning, then Jet Id
        self.noIdCleanJetsAll, cleanLeptons = cleanJetsAndLeptons(event.jetsAllNoID, leptons, self.jetLepDR, self.jetLepArbitration)
        self.noIdCleanJets = [j for j in self.noIdCleanJetsAll if abs(j.eta()) <  self.cfg_ana.jetEtaCentral ]
        self.noIdCleanJetsFwd = [j for j in self.noIdCleanJetsAll if abs(j.eta()) >=  self.cfg_ana.jetEtaCentral ]
        self.noIdDiscardedJets = [j for j in event.jetsAllNoID if j not in self.noIdCleanJetsAll]

        ## Clean Jets from photons (first cleaning, then Jet Id)
        photons = []
        if hasattr(event, 'selectedPhotons'):
            if self.cfg_ana.cleanJetsFromFirstPhoton:
                photons = event.selectedPhotons[:1]
            else:
                photons = [ g for g in event.selectedPhotons ] 

        self.gamma_cleanJetaAll = []
        self.gamma_noIdCleanJetsAll = []

        if self.cleanFromLepAndGammaSimultaneously:
            self.gamma_cleanJetsAll = cleanNearestJetOnly(jetsEtaCut, photons+leptons, self.jetGammaLepDR)
            self.gamma_noIdCleanJetsAll = cleanNearestJetOnly(event.jetsAllNoID, photons+leptons, self.jetGammaLepDR)
        else:
            self.gamma_cleanJetsAll = cleanNearestJetOnly(self.cleanJetsAll, photons, self.jetGammaDR)
            self.gamma_noIdCleanJetsAll = cleanNearestJetOnly(self.noIdCleanJetsAll, photons, self.jetGammaDR)

        self.gamma_cleanJets    = [j for j in self.gamma_cleanJetsAll if abs(j.eta()) <  self.cfg_ana.jetEtaCentral ]
        self.gamma_cleanJetsFwd = [j for j in self.gamma_cleanJetsAll if abs(j.eta()) >= self.cfg_ana.jetEtaCentral ]

        self.gamma_noIdCleanJets    = [j for j in self.gamma_noIdCleanJetsAll if abs(j.eta()) <  self.cfg_ana.jetEtaCentral ]
        self.gamma_noIdCleanJetsFwd = [j for j in self.gamma_noIdCleanJetsAll if abs(j.eta()) >= self.cfg_ana.jetEtaCentral ]
        ###

        ##And now for the cleaning from all photons
        ###
        gg_photons = [ g for g in event.selectedPhotons ] 
        self.gg_cleanJetsAll = []
        self.gg_noIdCleanJetsAll = []
        if self.cleanFromLepAndGammaSimultaneously:
            self.gg_cleanJetsAll = cleanNearestJetOnly(jetsEtaCut, gg_photons+leptons, self.jetGammaLepDR)
            self.gg_noIdCleanJetsAll = cleanNearestJetOnly(event.jetsAllNoID, gg_photons+leptons, self.jetGammaLepDR)
        else:
            self.gg_cleanJetsAll = cleanNearestJetOnly(self.cleanJetsAll, gg_photons, self.jetGammaDR)
            self.gg_noIdCleanJetsAll = cleanNearestJetOnly(self.noIdCleanJetsAll, gg_photons, self.jetGammaDR)

        self.gg_cleanJets    = [j for j in self.gg_cleanJetsAll if abs(j.eta()) <  self.cfg_ana.jetEtaCentral ]
        self.gg_cleanJetsFwd = [j for j in self.gg_cleanJetsAll if abs(j.eta()) >= self.cfg_ana.jetEtaCentral ]

        self.gg_noIdCleanJets    = [j for j in self.gg_noIdCleanJetsAll if abs(j.eta()) <  self.cfg_ana.jetEtaCentral ]
        self.gg_noIdCleanJetsFwd = [j for j in self.gg_noIdCleanJetsAll if abs(j.eta()) >= self.cfg_ana.jetEtaCentral ]
        ###

        if self.cfg_ana.alwaysCleanPhotons:
            self.cleanJets = self.gamma_cleanJets
            self.cleanJetsAll = self.gamma_cleanJetsAll
            self.cleanJetsFwd = self.gamma_cleanJetsFwd
            #
            self.noIdCleanJets = self.gamma_noIdCleanJets
            self.noIdCleanJetsAll = self.gamma_noIdCleanJetsAll
            self.noIdCleanJetsFwd = self.gamma_noIdCleanJetsFwd

        ## Jet Id, after jet/lepton cleaning
        self.cleanJetsFailIdAll = []
        for jet in self.noIdCleanJetsAll:
            if not self.testJetID( jet ):
                self.cleanJetsFailIdAll.append(jet)
        
        self.cleanJetsFailId = [j for j in self.cleanJetsFailIdAll if abs(j.eta()) <  self.cfg_ana.jetEtaCentral ]
        
        ## Jet Id, after jet/photon cleaning
        self.gamma_cleanJetsFailIdAll = []
        for jet in self.gamma_noIdCleanJetsAll:
            if not self.testJetID( jet ):
                self.gamma_cleanJetsFailIdAll.append(jet)

        self.gamma_cleanJetsFailId = [j for j in self.gamma_cleanJetsFailIdAll if abs(j.eta()) <  self.cfg_ana.jetEtaCentral ] 

        ## Jet Id, after jet/photon cleaning
        self.gg_cleanJetsFailIdAll = []
        for jet in self.gg_noIdCleanJetsAll:
            if not self.testJetID( jet ):
                self.gg_cleanJetsFailIdAll.append(jet)

        self.gg_cleanJetsFailId = [j for j in self.gg_cleanJetsFailIdAll if abs(j.eta()) <  self.cfg_ana.jetEtaCentral ]
       
        if self.cfg_comp.isMC:
            self.deltaMetFromJetSmearing = [0, 0]
            for j in self.cleanJetsAll:
                if hasattr(j, 'deltaMetFromJetSmearing'):
                    self.deltaMetFromJetSmearing[0] += j.deltaMetFromJetSmearing[0]
                    self.deltaMetFromJetSmearing[1] += j.deltaMetFromJetSmearing[1]

            self.cleanGenJets = cleanNearestJetOnly(event.genJets, leptons, self.jetLepDR)
            
            if self.cfg_ana.cleanGenJetsFromPhoton:
                if self.cleanFromLepAndGammaSimultaneously:
                    self.cleanGenJets = cleanNearestJetOnly(self.cleanGenJets, photons, self.jetLepDR)
                if self.cleanFromLepAndGammaSimultaneously:
                    self.cleanGenJets = cleanNearestJetOnly(self.cleanGenJets, photons, self.jetLepDR)

            if self.cfg_ana.do_mc_match:
                self.jetFlavour(event)

#        if hasattr(event,"cleanJetsAll"+self.cfg_ana.collectionPostFix): raise RuntimeError("Event already contains a jet collection with the following postfix: "+self.cfg_ana.collectionPostFix)

        setattr(event,"cleanJetsAll"           +self.cfg_ana.collectionPostFix, self.cleanJetsAll           ) 
        setattr(event,"cleanJets"              +self.cfg_ana.collectionPostFix, self.cleanJets              ) 
        setattr(event,"cleanJetsFwd"           +self.cfg_ana.collectionPostFix, self.cleanJetsFwd           ) 
        setattr(event,"cleanJetsFailIdAll"           +self.cfg_ana.collectionPostFix, self.cleanJetsFailIdAll           ) 
        setattr(event,"cleanJetsFailId"              +self.cfg_ana.collectionPostFix, self.cleanJetsFailId              ) 
        setattr(event,"discardedJets"          +self.cfg_ana.collectionPostFix, self.discardedJets          ) 
        setattr(event,"discardedLeptons"          +self.cfg_ana.collectionPostFix, self.discardedLeptons          ) 
        setattr(event,"selectedLeptons"          +self.cfg_ana.collectionPostFix, self.selectedLeptons          ) 
        setattr(event,"gamma_cleanJetsAll"     +self.cfg_ana.collectionPostFix, self.gamma_cleanJetsAll     ) 
        setattr(event,"gamma_cleanJets"        +self.cfg_ana.collectionPostFix, self.gamma_cleanJets        ) 
        setattr(event,"gamma_cleanJetsFwd"     +self.cfg_ana.collectionPostFix, self.gamma_cleanJetsFwd     ) 
        setattr(event,"gamma_cleanJetsFailIdAll"     +self.cfg_ana.collectionPostFix, self.gamma_cleanJetsFailIdAll     ) 
        setattr(event,"gamma_cleanJetsFailId"        +self.cfg_ana.collectionPostFix, self.gamma_cleanJetsFailId        ) 
        setattr(event,"gg_cleanJetsAll"     +self.cfg_ana.collectionPostFix, self.gg_cleanJetsAll     ) 
        setattr(event,"gg_cleanJets"        +self.cfg_ana.collectionPostFix, self.gg_cleanJets        ) 
        setattr(event,"gg_cleanJetsFwd"     +self.cfg_ana.collectionPostFix, self.gg_cleanJetsFwd     ) 
        setattr(event,"gg_cleanJetsFailIdAll"     +self.cfg_ana.collectionPostFix, self.gg_cleanJetsFailIdAll     ) 
        setattr(event,"gg_cleanJetsFailId"        +self.cfg_ana.collectionPostFix, self.gg_cleanJetsFailId        ) 


        if self.cfg_comp.isMC:
            setattr(event,"deltaMetFromJetSmearing"+self.cfg_ana.collectionPostFix, self.deltaMetFromJetSmearing) 
            setattr(event,"cleanGenJets"           +self.cfg_ana.collectionPostFix, self.cleanGenJets           )
 
        return True

        

    def testJetID(self, jet):
        jet.puJetIdPassed = jet.puJetId() 
        jet.pfJetIdPassed = jet.jetID('POG_PFID_Loose') 
        if self.cfg_ana.relaxJetId:
            return True
        else:
            return jet.pfJetIdPassed and (jet.puJetIdPassed or not(self.doPuId)) 
        
    def testJetNoID( self, jet ):
        # 2 is loose pile-up jet id
        return jet.pt() > self.cfg_ana.jetPt and \
               abs( jet.eta() ) < self.cfg_ana.jetEta;

    def jetFlavour(self,event):
        def isFlavour(x,f):
            id = abs(x.pdgId())
            if id > 999: return (id/1000)%10 == f
            if id >  99: return  (id/100)%10 == f
            return id % 100 == f


        self.partons   = [ p for p in event.genParticles if ((p.status() == 23 or p.status() == 3) and abs(p.pdgId())>0 and (abs(p.pdgId()) in [1,2,3,4,5,21]) ) ]
        match = matchObjectCollection2(self.cleanJetsAll,
                                       self.partons,
                                       deltaRMax = 0.3)

        for jet in self.cleanJetsAll:
            parton = match[jet]
            jet.partonId = (parton.pdgId() if parton != None else 0)
            jet.partonMotherId = (parton.mother(0).pdgId() if parton != None and parton.numberOfMothers()>0 else 0)
        

setattr(JetCleaner,"defaultConfig", cfg.Analyzer(
    class_object = JetCleaner,
    jetPt = 25.,
    jetEta = 4.7,
    jetEtaCentral = 2.4,
    jetLepDR = 0.4,
    jetLepArbitration = (lambda jet,lepton : lepton), # you can decide which to keep in case of overlaps; e.g. if the jet is b-tagged you might want to keep the jet
    cleanSelectedLeptons = True, #Whether to clean 'selectedLeptons' after disambiguation. Treat with care (= 'False') if running Jetanalyzer more than once
    minLepPt = 10,
    lepSelCut = lambda lep : True,
    relaxJetId = False,  
    doPuId = False, # Not commissioned in 7.0.X
    cleanJetsFromFirstPhoton = False,
    cleanJetsFromTaus = False,
    cleanJetsFromIsoTracks = False,
    alwaysCleanPhotons = False,
    do_mc_match=True,
    cleanGenJetsFromPhoton = False,
    jetGammaDR=0.4,
    cleanFromLepAndGammaSimultaneously = False,
    jetGammaLepDR=0.4,
    collectionPostFix = ""
    )
)
