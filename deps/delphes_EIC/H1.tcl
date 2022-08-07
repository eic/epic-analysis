set ExecutionPath {
  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency


  ChargedHadronMomentumSmearing
  ElectronMomentumSmearing
  TrackMerger
  

  ECal
  HCal

  Calorimeter
  EFlowMerger
  EFlowFilter

  PhotonEfficiency
  PhotonIsolation

  ElectronFilter
  ElectronEfficiency
  ElectronIsolation

  ChargedHadronFilter
  MissingET
  NeutrinoFilter
  GenJetFinder
  GenMissingET
  FastJetFinder

  JetEnergyScale
  ElectronEnergyScale
    
  JetFlavorAssociation
  GenJetFlavorAssociation

  UniqueObjectFinder

  ScalarHT

  

  TreeWriter
}


module ParticlePropagator ParticlePropagator {
    set InputArray Delphes/stableParticles
    set OutputArray stableParticles
    set ChargedHadronOutputArray chargedHadrons
    set ElectronOutputArray electrons

    # radius of the magnetic field coverage, in m                                                          
    set Radius 1.5
    # half-length of the magnetic field coverage, in m                                                     
    set HalfLength 1.20
    # magnetic field                                                                                       
    set Bz 1.16
}

## H1 central tracker covered                                                                              
set CommonTrackingEfficiency {
    (abs(eta) <= 2.0) * (pt > 0.100)                     * (1.0) +
    0.0
}

set CommonTrackingResolution {
    (abs(eta)<=1.0)  * (sqrt( (2.0e-2)^2 + (pt*cosh(eta)*5e-3)^2  ) ) +
    (abs(eta)>1.0 && abs(eta)<2.0)  * (sqrt( (10.0e-2)^2 + (pt*cosh(eta)*1e-2)^2  ) )
}


module Efficiency ChargedHadronTrackingEfficiency {
  set InputArray ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons
  set EfficiencyFormula $CommonTrackingEfficiency
}

##############################                                                                             
# Electron tracking efficiency                                                                             
##############################                                                                             

module Efficiency ElectronTrackingEfficiency {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons
  set EfficiencyFormula $CommonTrackingEfficiency

}
########################################                                                                   
# Momentum resolution for charged tracks                                                                   
########################################                                                                   

module MomentumSmearing ChargedHadronMomentumSmearing {
  set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons
  set ResolutionFormula  $CommonTrackingResolution
}
###################################                                                                        
# Momentum resolution for electrons                                                                        
###################################                                                                        
module MomentumSmearing ElectronMomentumSmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons
  set ResolutionFormula $CommonTrackingResolution
}


##############                                                                                             
# Track merger                                                                                             
##############                                                                                             

module Merger TrackMerger {
# add InputArray InputArray                                                                                
  add InputArray ChargedHadronMomentumSmearing/chargedHadrons
  add InputArray ElectronMomentumSmearing/electrons

  set OutputArray tracks
}

################################                                                                           
# Track impact parameter smearing                                                                          
################################                                                                           

#module TrackSmearing TrackSmearing {

# set InputArray TrackMerger/tracks
#  set BeamSpotInputArray BeamSpotFilter/beamSpotParticle                                                  
#  set OutputArray tracks
#  set ApplyToPileUp true                                                                                  
  # magnetic field                                                                                         
#  set Bz 1.16
#  set PResolutionFormula { 0.0 }
#  set CtgThetaResolutionFormula { 0.0 }
#  set PhiResolutionFormula { 0.0 }
#  set D0ResolutionFormula "0.02"
#  set DZResolutionFormula "0.02 "

#}



#############                                                                                              
#   ECAL                                                                                                   
#############                                                                                              

module SimpleCalorimeter ECal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray TrackMerger/tracks

  set TowerOutputArray ecalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowPhotons

  set IsEcal true
  set EnergyMin 0.10



  set EnergySignificanceMin 1.0

  set SmearTowerCenter true
  set pi [expr {acos(-1)}]
  set PhiBins {}
  for {set i -30} {$i <=30} {incr i} {
      add PhiBins [expr {$i * $pi/30.0}]
  }
  for {set i -10} {$i <=10} {incr i} {
      set eta [expr {$i * 0.1}]
      add EtaPhiBins $eta $PhiBins
  }
  set PhiBins {}
  for {set i -30} {$i <=30} {incr i} {
    add PhiBins [expr {$i * $pi/30.0}]
  }

  foreach eta {-3.3 -3.26996837 -3.14642305 -3.03653567 -2.93760447 -2.84766006 -2.76522251 -2.68915144 \
-2.61854952 -2.55269788 -2.49101173 -2.43300894 -2.3782873  -2.3265078  -2.27738197 -2.23066235 -2.1861350\
3 -2.14361383 -2.10293569 -2.063957 -2.02655061 -1.99060337 -1.95601417 -1.92269228 -1.89055593 -1.8595312\
  -1.82955102 -1.80055436 -1.77248548 -1.74529337 -1.71893119 -1.69335587 -1.66852765 -1.64440978 -1.62096\
821 -1.59817135 -1.57598979 -1.55439612 -1.53336478 -1.51287184 -1.4928949  -1.47341295 -1.45440623 -1.435\
85618 -1.41774529 -1.40005705 -1.38277588 -1.36588703 -1.34937654 -1.33323117 -1.31743839 -1.30198626 -1.2\
8686345 -1.27205918 -1.25756317 -1.24336562 -1.22945719 -1.21582897 -1.20247241 -1.18937936 -1.17654201 -1\
.16395288 -1.15160481 -1.13949092 -1.12760462 -1.11593955 -1.10448965 -1.09324904 -1.08221211     -1.07137\
342 -1.06072776 -1.0502701  -1.03999558} {
    add EtaPhiBins $eta $PhiBins
}
  foreach eta {1.0 1.0502701  1.06072776 1.07137342 1.08221211 1.09324904 1.10448965 1.11593955 1.127604\
62 1.13949092 1.15160481 1.16395288 1.17654201 1.18937936 1.20247241 1.21582897 1.22945719 1.24336562 1.25\
756317 1.27205918 1.28686345 1.30198626 1.31743839 1.33323117 1.34937654 1.36588703 1.38277588 1.40005705 \
1.41774529 1.43585618 1.45440623 1.47341295 1.4928949  1.51287184 1.53336478 1.55439612 1.57598979 1.59817\
135 1.62096821 1.64440978 1.66852765 1.69335587 1.71893119 1.74529337 1.77248548 1.80055436 1.82955102 1.8\
595312 1.89055593 1.92269228 1.95601417 1.99060337 2.02655061 2.063957 2.10293569 2.14361383 2.18613503 2.\
23066235 2.27738197 2.3265078 2.3782873  2.43300894 2.49101173 2.55269788 2.61854952 2.68915144 2.76522251\
 2.84766006 2.93760447 3.03653567 3.14642305 3.26996837 3.3} {
    add EtaPhiBins $eta $PhiBins
}
    
  add EnergyFraction {0} {0.0}
  # energy fractions for e, gamma and pi0                                                                  
  add EnergyFraction {11} {1.0}
  add EnergyFraction {22} {1.0}
  add EnergyFraction {111} {1.0}
  # energy fractions for muon, neutrinos and neutralinos                                                   
  add EnergyFraction {12} {0.0}
  add EnergyFraction {13} {0.0}
  add EnergyFraction {14} {0.0}
  add EnergyFraction {16} {0.0}
  add EnergyFraction {1000022} {0.0}
  add EnergyFraction {1000023} {0.0}
  add EnergyFraction {1000025} {0.0}
  add EnergyFraction {1000035} {0.0}
  add EnergyFraction {1000045} {0.0}
# energy fractions for K0short and Lambda                                                                
  add EnergyFraction {310} {0.3}
  add EnergyFraction {3122} {0.3}


   set ResolutionFormula {
       ( eta> -1.46 && eta < 3.35 )  * sqrt(energy^2*0.025^2 + energy*0.11^2 )
     + ( eta> -3.35 && eta < -1.46 ) * sqrt(energy^2*0.030^2 + energy*0.10^2 )
   }

}

#############                                                                                              
#   HCAL                                                                                                   
#############                                                                                              

module SimpleCalorimeter HCal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray ECal/eflowTracks

  set TowerOutputArray hcalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowNeutralHadrons

  set IsEcal false

  ##Assumes noise 100 MeV per tower.                                                                       
  set EnergyMin 0.5
  set EnergySignificanceMin 1.0

  set SmearTowerCenter true

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  set PhiBins {}
  for {set i -30} {$i <=30} {incr i} {
      add PhiBins [expr {$i * $pi/30.0}]
  }
  for {set i -10} {$i <=10} {incr i} {
      set eta [expr {$i * 0.1}]
      add EtaPhiBins $eta $PhiBins
  }

  for {set i -30} {$i <=30} {incr i} {
      add PhiBins [expr {$i * $pi/30.0}]
  }

  foreach eta {-3.3 -2.95880652 -2.68264484 -2.46773612 -2.29224349 -2.14432155 -2.01681569 -1.90506801 \
-1.80587261 -1.71692581 -1.63651428 -1.56332731 -1.49633825 -1.43472677 -1.37782606 -1.325086   -1.2760468\
4 -1.23031998 -1.18757364 -1.14752205 -1.10991713 -1.07454199 -1.04120583 -1.00} {
      add EtaPhiBins $eta $PhiBins
  }

  foreach eta {1.0 1.04 1.075 1.1099 1.14752205 1.18757364 1.23031998 1.27604684 1.325086 1.37782606 1.4\
3472677 1.49633825 1.56332731 1.63651428 1.71692581 1.80587261 1.90506801 2.01681569 2.14432155 2.29224349\
 2.46773612 2.68264484 2.95880652 3.3} {
      add EtaPhiBins $eta $PhiBins
  }
    
  add EnergyFraction {0} {1.0}
  # energy fractions for e, gamma and pi0                                                                  
  add EnergyFraction {11} {0.0}
  add EnergyFraction {22} {0.0}
  add EnergyFraction {111} {0.0}
  # energy fractions for muon, neutrinos and neutralinos                                                   
  add EnergyFraction {12} {0.0}
  add EnergyFraction {13} {0.0}
  add EnergyFraction {14} {0.0}
  add EnergyFraction {16} {0.0}
  add EnergyFraction {1000022} {0.0}
  add EnergyFraction {1000023} {0.0}
  add EnergyFraction {1000025} {0.0}
  add EnergyFraction {1000035} {0.0}
  add EnergyFraction {1000045} {0.0}

  add EnergyFraction {310} {0.7}
  add EnergyFraction {3122} {0.7}

  # set HCalResolutionFormula {resolution formula as a function of eta and energy}                         

  set ResolutionFormula {
       ( eta> -0.64 && eta <  3.20 )  * sqrt(energy^2*0.20^2 + energy*0.50^2)
     + ( eta>  3.20 && eta <  3.35 )  * sqrt(energy^2*0.40^2 + energy*0.90^2)
     + ( eta> -0.97 && eta < -0.64 )  * sqrt(energy^2*0.40^2 + energy*0.90^2)
  }
}

#################                                                                                          
# Electron filter                                                                                          
#################                                                                                          

module PdgCodeFilter ElectronFilter {
  set InputArray HCal/eflowTracks
  set OutputArray electrons
  set Invert true
  add PdgCode {11}
  add PdgCode {-11}
}
######################                                                                                     
# ChargedHadronFilter                                                                                      
######################                                                                                     

module PdgCodeFilter ChargedHadronFilter {
  set InputArray HCal/eflowTracks
  set OutputArray chargedHadrons

  add PdgCode {11}
  add PdgCode {-11}
  add PdgCode {13}
  add PdgCode {-13}
}


###################################################                                                        
# Tower Merger (in case not using e-flow algorithm)                                                        
###################################################                                                        

module Merger Calorimeter {
# add InputArray InputArray                                                                                
  add InputArray ECal/ecalTowers
  add InputArray HCal/hcalTowers
  set OutputArray towers
}
####################                                                                                       
# Energy flow merger                                                                                       
####################                                                                                       

module Merger EFlowMerger {
# add InputArray InputArray                                                                                
  add InputArray HCal/eflowTracks
  add InputArray ECal/eflowPhotons
  add InputArray HCal/eflowNeutralHadrons
  set OutputArray eflow
}

######################                                                                                     
# EFlowFilter                                                                                              
######################                                                                                     

module PdgCodeFilter EFlowFilter {
  set InputArray EFlowMerger/eflow
  set OutputArray eflow

  add PdgCode {11}
  add PdgCode {-11}
  add PdgCode {13}
  add PdgCode {-13}
}
###################                                                                                        
# Photon efficiency                                                                                        
###################                                                                                        

module Efficiency PhotonEfficiency {
  set InputArray ECal/eflowPhotons
  set OutputArray photons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}                                 

  # efficiency formula for photons                                                                         
    set EfficiencyFormula { 1}
}
module Isolation PhotonIsolation {
  set CandidateInputArray PhotonEfficiency/photons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray photons

  set DeltaRMax 0.5

  set PTMin 0.5

  set PTRatioMax 0.12
}
module Efficiency ElectronEfficiency {
  set InputArray ElectronFilter/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}                                 

  # efficiency formula for electrons                                                                       
    set EfficiencyFormula {1}
}

####################                                                                                       
# Electron isolation                                                                                       
####################                                                                                       

module Isolation ElectronIsolation {
  set CandidateInputArray ElectronEfficiency/electrons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray electrons

  set DeltaRMax 0.5

  set PTMin 0.5
  set PTRatioMax 0.12
}

###################                                                                                        
# Missing ET merger                                                                                        
###################                                                                                        

module Merger MissingET {
# add InputArray InputArray                                                                                
  add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
}
##################                                                                                         
# Scalar HT merger                                                                                         
##################                                                                                         

module Merger ScalarHT {
# add InputArray InputArray                                                                                
  add InputArray UniqueObjectFinder/jets
  add InputArray UniqueObjectFinder/electrons
  add InputArray UniqueObjectFinder/photons

  set EnergyOutputArray energy
}
#####################                                                                                      
# Neutrino Filter                                                                                          
#####################                                                                                      
module PdgCodeFilter NeutrinoFilter {
  set InputArray Delphes/stableParticles
  set OutputArray filteredParticles
  set PTMin 0.0
  add PdgCode {12}
  add PdgCode {14}
  add PdgCode {16}
  add PdgCode {-12}
  add PdgCode {-14}
  add PdgCode {-16}

}
#####################                                                                                      
# MC truth jet finder                                                                                      
#####################                                                                                      

module FastJetFinder GenJetFinder {
  set InputArray NeutrinoFilter/filteredParticles

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt                      
  set JetAlgorithm 6
  set ParameterR 1.0

  set JetPTMin 3.0
}
#########################                                                                                  
# Gen Missing ET merger                                                                                    
########################                                                                                   

module Merger GenMissingET {
# add InputArray InputArray                                                                                
  add InputArray NeutrinoFilter/filteredParticles
  set MomentumOutputArray momentum
}
############                                                                                               
# Jet finder                                                                                               
############                                                                                               

module FastJetFinder FastJetFinder {
#  set InputArray Calorimeter/towers                                                                       
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt                      
  set JetAlgorithm 6
  set ParameterR 1.0

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4
set ComputeTrimming 1
  set RTrim 0.4
  set PtFracTrim 0.20
  #set PtFracTrim 0.05                                                                                     

  set ComputePruning 1
  set ZcutPrun 0.1
  set RcutPrun 0.5
  set RPrun 0.8

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.8

 set JetPTMin 3.0}
##################                                                                                         
# Jet Energy Scale                                                                                         
##################                                                                                         

module EnergyScale JetEnergyScale {
  set InputArray FastJetFinder/jets
  set OutputArray jets

  # scale formula for jets (do not apply it)                                                               
  set ScaleFormula {1.0}
}

module EnergyScale ElectronEnergyScale {
  set InputArray ElectronIsolation/electrons
  set OutputArray electrons
  set ScaleFormula {0.995}
}


########################                                                                                   
# Jet Flavor Association                                                                                   
########################                                                                                   

module JetFlavorAssociation JetFlavorAssociation {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScale/jets

  set DeltaR 0.5
  set PartonPTMin 4.0
  set PartonEtaMax 4.0

}
module JetFlavorAssociation GenJetFlavorAssociation {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray GenJetFinder/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 4.0

}
#####################################################                                                      
# Find uniquely identified photons/electrons/tau/jets                                                      
#####################################################                                                      

module UniqueObjectFinder UniqueObjectFinder {
# earlier arrays take precedence over later ones                                                           
# add InputArray InputArray OutputArray                                                                    
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronEnergyScale/electrons electrons
  add InputArray JetEnergyScale/jets jets
}
module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass                                                            
  add Branch Delphes/allParticles Particle GenParticle

  add Branch TrackMerger/tracks Track Track
  add Branch Calorimeter/towers Tower Tower

  add Branch HCal/eflowTracks EFlowTrack Track
  add Branch ECal/eflowPhotons EFlowPhoton Tower
  add Branch HCal/eflowNeutralHadrons EFlowNeutralHadron Tower

  add Branch GenJetFinder/jets GenJet Jet
  add Branch GenMissingET/momentum GenMissingET MissingET

  add Branch UniqueObjectFinder/jets Jet Jet
  add Branch UniqueObjectFinder/electrons Electron Electron
  add Branch UniqueObjectFinder/photons Photon Photon

  add Branch MissingET/momentum MissingET MissingET
  add Branch ScalarHT/energy ScalarHT ScalarHT
}


