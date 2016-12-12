from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.physicsutils.genutils import isNotFromHadronicShower, realGenMothers, realGenDaughters, motherRef

import PhysicsTools.HeppyCore.framework.config as cfg

        
class TauDecayModeAnalyzer( Analyzer ):
    """Classify and filter events according to Tau lepton decays

       Reads:
         event.gentaus

       Creates in the event:
         event.gentau_decayMode  = 11 for tau-->enunu
                                   13 for tau-->mununu
                                   1 for 1-prong hadronic decays
                                   3 for 3-prong hadronic decays
         event.gentau_neutralDaughters = number of neutral (non-neutrino) daughters     
         event.gentau_leadTrackPt = pT of the leading charged decay product
    """
    def __init__(self, cfg_ana, cfg_comp, looperName ):
        super(TauDecayModeAnalyzer,self).__init__(cfg_ana,cfg_comp,looperName)
    #---------------------------------------------
        

    def declareHandles(self):
        super(TauDecayModeAnalyzer, self).declareHandles()

    def beginLoop(self, setup):
        super(TauDecayModeAnalyzer,self).beginLoop(setup)

    def process(self, event):
        self.readCollections( event.input )

        # if not MC, nothing to do
        if not self.cfg_comp.isMC: 
            return True

#        event.gentaus.decayMode        = []
#        event.gentaus.leadTrackPt      = []
#        event.gentaus.neutralDaughters = []
        
        event.ngenTau1Prong = 0
        event.ngenTau3Prong = 0

        for index,p in enumerate(event.gentaus):
             
            leadTrackPt = float(0.)
            decayMode   = 0
            neutralDaughters = 0
            
            if p.status() !=2 and p.status()!=23:
                continue;
            if abs(p.motherId) != 25 and abs(p.motherId) != 24 and abs(p.motherId) != 23 and abs(p.motherId) != 21:
                continue;

            if p.status() == 23:
                for p2 in xrange(p.numberOfDaughters()):
                    if p.daughter(p2).status()==2 and abs(p.daughter(p2).pdgId())==15:
                        for d in xrange(p.daughter(p2).numberOfDaughters()):
                            if p.daughter(p2).daughter(d).status()!=1 and p.daughter(p2).daughter(d).status()!=2:
                                continue
                            if abs(p.daughter(p2).daughter(d).pdgId())==15:
                                continue
                            if abs(p.daughter(p2).daughter(d).pdgId())==11 or abs(p.daughter(p2).daughter(d).pdgId())==13:
                                decayMode   = p.daughter(p2).daughter(d).pdgId()
                                leadTrackPt = p.daughter(p2).daughter(d).pt()
                                break

                            elif abs(p.daughter(p2).daughter(d).threeCharge())>0.:
                                decayMode = decayMode+1
                                if p.daughter(p2).daughter(d).pt() > leadTrackPt:
                                    leadTrackPt = p.daughter(p2).daughter(d).pt()
                            
                            elif p.daughter(p2).daughter(d).threeCharge()==0 and abs(p.daughter(p2).daughter(d).pdgId())!=12 and abs(p.daughter(p2).daughter(d).pdgId())!=14 and abs(p.daughter(p2).daughter(d).pdgId())!=16:
                                neutralDaughters = neutralDaughters+1
                        break
                    else:
                        continue

            elif p.status() == 2:
                for d in xrange(p.numberOfDaughters()):
                    if p.daughter(d).status()!=1 and p.daughter(d).status()!=2:
                        continue
                    if abs(p.daughter(d).pdgId())==15:
                        continue
                    if abs(p.daughter(d).pdgId())==11 or abs(p.daughter(d).pdgId())==13:
                        decayMode   = p.daughter(d).pdgId()
                        leadTrackPt = p.daughter(d).pt()
                        break
                    
                    elif abs(p.daughter(d).threeCharge())>0:
                        decayMode = decayMode+1
                        if p.daughter(d).pt() > leadTrackPt:
                            leadTrackPt = p.daughter(d).pt()
                            
                    elif p.daughter(d).threeCharge()==0 and abs(p.daughter(d).pdgId())!=12 and abs(p.daughter(d).pdgId())!=14 and abs(p.daughter(d).pdgId())!=16:
                        neutralDaughters = neutralDaughters+1

            else:
                continue
            
            
            p.decayMode = decayMode
            p.leadTrackPt = leadTrackPt
            p.neutralDaughters = neutralDaughters
            
            if p.decayMode == 1:
                event.ngenTau1Prong = event.ngenTau1Prong + 1
            elif p.decayMode == 3:
                event.ngenTau3Prong = event.ngenTau3Prong + 1
            
        return True

setattr(TauDecayModeAnalyzer,"defaultConfig",
    cfg.Analyzer(TauDecayModeAnalyzer,
    )
)
