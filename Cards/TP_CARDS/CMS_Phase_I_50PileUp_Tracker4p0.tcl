#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {

  PileUpMerger

  ModifyBeamSpot
  GenBeamSpotFilter

  ParticlePropagator

  StatusPid

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  ChargedHadronMomentumSmearing
  ElectronEnergySmearing
  MuonMomentumSmearing

  TrackMerger

  Calorimeter

  TrackPileUpSubtractor

  EFlowMerger

  GlobalRhoKt4
  GlobalRhoGridFastJet

  RhoKt4
  RhoGridFastJet

  FastJetFinder
  TrackJetFinder
  NeutrinoFilter
  GenJetFinderNoNu

  JetPileUpSubtractor
  JetFlavourAssociation

  BTagging

  PileUpJetID

  PhotonEfficiency
  PhotonIsolation

  ElectronEfficiency
  ElectronIsolation

  MuonEfficiency
  MuonIsolation

  RunPUPPI
  PuppiRhoKt4
  PuppiRhoGridFastJet
  PuppiJetFinder

  PuppiJetPileUpSubtractor
  PuppiJetFlavourAssociation

  PuppiBTagging

  PuppiPileUpJetID

  GenMissingET
  MissingET
  PuppiMissingET

  TreeWriter

}

#  JetPileUpSubtractorGrid
#  JetPileUpSubtractor4VArea
#  PuppiJetPileUpSubtractorGrid
#  PuppiJetPileUpSubtractor4VArea

#### remove the module which do the filter of jet constituent                                                                                                                    
# ConstituentFilter                                                                                                                                                          
# PuppiConstituentFilter                                                                                                                                                               
## unique object finder actually removed from the sequence                                                                                                                           
# UniqueObjectFinderGJ                                                                                                                                                                  
# UniqueObjectFinderEJ                                                                                                                                                               
# UniqueObjectFinderMJ                                                                                                                                                              


#################                                                                                                                                                                       
# PileUp Merger #                                                                                                                                                                      
#################                                                                                                                                                                        
module PileUpMerger PileUpMerger {
 ## inputs are status 1 HEPMC particles --> real final state one                                                                                                                         
 set InputArray  Delphes/stableParticles
 ## output array is called stable particles                                                                                                                                             
 set OutputArray stableParticles
 ## store NPU                                                                                                                                                                              
 set NPUOutputArray NPU
 # Get rid of beam spot from http://red-gridftp11.unl.edu/Snowmass/MinBias100K_14TeV.pileup ...                                                                                           
 set InputBSX 2.44
 set InputBSY 3.39
 # replace it with beam spot from CMSSW files                                                                                                                                             
 set OutputBSX 0.24
 set OutputBSY 0.39
 # pre-generated minbias input file --> change this dummy name <random access with unifor number between 0 and NEntries>                                                                 
 set PileUpFile /afs/cern.ch/user/s/spigazzi/work/public/PU_14TeV/MinBias_14TeV_100k_TunePP15.pileup
 #average expected pile up <poissonian generation>                                                                                                                                       
 set MeanPileUp 50
 # spread in the beam direction in m (assumes gaussian) ;                                                                                                                                
 set ZVertexSpread 0.053
}

##################                                                                                                                                                           
# ModifyBeamSpot #                                                                                                                                                      
##################                                                                                                                                                                 
module ModifyBeamSpot ModifyBeamSpot {
  set ZVertexSpread 0.053
  set InputArray    PileUpMerger/stableParticles
  set OutputArray   stableParticles
  set PVOutputArray PV
}

#####################                                                                                                                                                          
# GenBeamSpotFilter #                                                                                                                                                             
#####################                                                                                                                                                              
module GenBeamSpotFilter GenBeamSpotFilter {
    set InputArray ModifyBeamSpot/stableParticles
    set OutputArray beamSpotParticles
}

#####################################################################################                                                                                             
# Propagate particles in cylinder and divide charged particles in different classes #                                                                                       
#####################################################################################                                                                                                
module ParticlePropagator ParticlePropagator {
  ## take particles after beam spot smearing                                                                                                                                         
  set InputArray ModifyBeamSpot/stableParticles
  ## produce independent output collection: all particles, only charged hadrons, electrons and muons                                                                                     
  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons
  ## radius of the magnetic field coverage, in m                                                                                                                                    
  set Radius 1.29
  ## half-length of the magnetic field coverage, in m                                                                                                                            
  set HalfLength 3.00
  ## magnetic field                                                                                                                                                             
  set Bz 3.8
}

###############################################################################################################                                                                     
# StatusPidFilter: this module removes all generated particles except electrons, muons, taus, and status == 3 #                                                                     
###############################################################################################################                                                                         
module StatusPidFilter StatusPid {
    ## take the particles from Pythia8 not adding pile-up                                                                                                                            
    set InputArray  Delphes/allParticles
    set OutputArray filteredParticles
    set PTMin 0.35
}


####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
  set InputArray ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons

  # add EfficiencyFormula {efficiency formula as a function of eta and pt} - Phase II
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.85) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0)                  * (0.97) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.85) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0)                  * (0.90) + \
                         (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.8) + \
                         (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0) * (0.82) + \
                         (abs(eta) > 4.0) * (0.00)
  }
}

##############################
# Electron tracking efficiency 
##############################

module Efficiency ElectronTrackingEfficiency {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons
  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
  # tracking efficiency formula for electrons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.85) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0)                  * (0.97) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.85) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0)                  * (0.90) + \
                         (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.8) + \
                         (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0 && pt <= 10.0) * (0.8+pt*0.01) + \
                         (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 10.0) * (0.85) + \
                         (abs(eta) > 4.0) * (0.00)
  }
}

##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
  # tracking efficiency formula for muons
    set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) + \
                                            (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.998) + \
                                            (abs(eta) <= 1.5) * (pt > 1.0)                 * (0.9998) + \
                           (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.98) + \
                           (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0)                  * (0.98) + \
                           (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.95) + \
                           (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0) * (0.95) + \
                           (abs(eta) > 4.0) * (0.00)
  }
}

########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons

  # set ResolutionFormula {resolution formula as a function of eta and pt}
  set ResolutionFormula {                  (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e1) * (0.013) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.02) + \
                                           (abs(eta) <= 1.5) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e1) * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.04) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 0.1   && pt <= 1.0)   * (0.017) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 1.0   && pt <= 10.0)  * (0.03) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 10.0  && pt <= 100.0) * (0.05) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 100.0)                * (0.30) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 0.1   && pt <= 1.0)   * (0.02) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 1.0   && pt <= 10.0)  * (0.04) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 10.0  && pt <= 100.0) * (0.07) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 100.0)                * (0.30) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 0.1   && pt <= 1.0)   * (0.025) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 1.0   && pt <= 10.0)  * (0.05) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 10.0  && pt <= 100.0) * (0.20) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 100.0)                * (0.80)
    }

}

#################################
# Energy resolution for electrons
#################################

module EnergySmearing ElectronEnergySmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons

  # set ResolutionFormula {resolution formula as a function of eta and energy}
  set ResolutionFormula {  (abs(eta) <= 2.5) * (energy > 0.1   && energy <= 2.5e1) * (energy*0.025) + \
                           (abs(eta) <= 2.5) * (energy > 2.5e1)                    * (energy*0.035) + \
                           (abs(eta) > 2.5 && abs(eta) <= 3.0)                     * (energy*0.035) + \
         		   (abs(eta) > 3.0 && abs(eta) <= 5.0)                     * (energy*0.07)
    }
}

###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons

  # set ResolutionFormula {resolution formula as a function of eta and pt}
  # resolution formula for muons
  set ResolutionFormula {  (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e1) * (0.012) + \
                           (abs(eta) <= 1.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.015) + \
                           (abs(eta) <= 1.5) * (pt > 2.0e2)                * (0.03) + \
                           (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                           (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e1) * (0.015) + \
                           (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.025) + \
			   (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 2.0e2)                * (0.03)  + \
                           (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 0.1   && pt <= 1.0)   * (0.017) + \
                           (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 1.0   && pt <= 10.0)  * (0.03) + \
                           (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 10.0  && pt <= 100.0) * (0.05) + \
                           (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 100.0)                * (0.30) + \
                           (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 0.1   && pt <= 1.0)   * (0.02) + \
                           (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 1.0   && pt <= 10.0)  * (0.04) + \
                           (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 10.0  && pt <= 100.0) * (0.07) + \
                           (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 100.0)                * (0.30) + \
                           (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 0.1   && pt <= 1.0)   * (0.025) + \
                           (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 1.0   && pt <= 10.0)  * (0.05) + \
                           (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 10.0  && pt <= 100.0) * (0.20) + \
                           (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 100.0)                * (0.80)
  }
  
}


#################################################################                                                                                                              
# Track merger : merge two collection of object in a output one #                                                                                                                
#################################################################                                                                                                                   
module Merger TrackMerger {
  ## take smeared charged hadron and electrons                                                                                                                                          
  add InputArray ChargedHadronMomentumSmearing/chargedHadrons
  add InputArray ElectronEnergySmearing/electrons
  set OutputArray tracks
}

#############
# Calorimeter
#############

module Calorimeter Calorimeter {
  ## particle from the propagation without any efficiency or smearing (for neutrals)                                                                                                  
  set ParticleInputArray ParticlePropagator/stableParticles
  ## track after smearing and efficiency: used for charged particles                                                                                                                     
  set TrackInputArray   TrackMerger/tracks
  ## output collections                                                                                                                                                                  
  set TowerOutputArray  towers
  set PhotonOutputArray photons
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowTowers

  set pi [expr {acos(-1)}]

  ## Granularity for |eta| < 1.65                                                                                                                                                       
  set PhiBins {}
  for {set i -36} {$i <= 36} {incr i} {
      add PhiBins [expr {$i * $pi/36.0}]
  }
  foreach eta {-1.566 -1.479 -1.392 -1.305 -1.218 -1.131 -1.044 -0.957 -0.87 -0.783 -0.696 -0.609 -0.522 -0.435 -0.348 -0.261 -0.174 -0.087 0 0.087 0.174 0.261 0.348 0.435 0.522 0.609 0.696 0.783 0.87 0.957 1.044 1.131 1.218 1.305 1.392 1.479 1.566 1.653} {
    add EtaPhiBins $eta $PhiBins
  }

  ## Granularity for 1.65 < |eta| < 4.5                                                                                                                                                  
  set PhiBins {}
  for {set i -18} {$i <= 18} {incr i} {
      add PhiBins [expr {$i * $pi/18.0}]
  }
  foreach eta {-4.35 -4.175 -4 -3.825 -3.65 -3.475 -3.3 -3.125 -2.95 -2.868 -2.65 -2.5 -2.322 -2.172 -2.043 -1.93 -1.83 -1.74 -1.653 1.74 1.83 1.93 2.043 2.172 2.322 2.5 2.65 2.868 2.95 3.125 3.3 3.475 3.65 3.825 4 4.175 4.35 4.525} {
    add EtaPhiBins $eta $PhiBins
  }

  ## Granularity for 4.5 < |eta| < 5                                                                                                                                                   
  set PhiBins {}
  for {set i -9} {$i <= 9} {incr i} {
      add PhiBins [expr {$i * $pi/9.0}]
  }
  foreach eta {-5 -4.7 -4.525 4.7 5} {
    add EtaPhiBins $eta $PhiBins
  }



  # default energy fractions {abs(PDG code)} {Fecal Fhcal}
  add EnergyFraction {0} {0.0 1.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {1.0 0.0}
  add EnergyFraction {22} {1.0 0.0}
  add EnergyFraction {111} {1.0 0.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0 0.0}
  add EnergyFraction {13} {0.0 0.0}
  add EnergyFraction {14} {0.0 0.0}
  add EnergyFraction {16} {0.0 0.0}
  add EnergyFraction {1000022} {0.0 0.0}
  add EnergyFraction {1000023} {0.0 0.0}
  add EnergyFraction {1000025} {0.0 0.0}
  add EnergyFraction {1000035} {0.0 0.0}
  add EnergyFraction {1000045} {0.0 0.0}
  # energy fractions for K0short and Lambda
  add EnergyFraction {310} {0.3 0.7}
  add EnergyFraction {3122} {0.3 0.7}

  # set ECalResolutionFormula {resolution formula as a function of eta and energy}                                                                                                      
  set ECalResolutionFormula { (abs(eta) <= 3.0) * sqrt(energy^2*0.005^2 + energy*0.027^2 + 0.15^2) + \
			      (abs(eta) > 3.0 && abs(eta) <= 5.0) * sqrt(energy^2*0.08^2 + energy*1.97^2)
  }


  # set HCalResolutionFormula {resolution formula as a function of eta and energy}                                                                                                       
  set HCalResolutionFormula {                  (abs(eta) <= 1.7) * sqrt(energy^2*0.0302^2 + energy*0.5205^2 + 1.59^2) + \
                             (abs(eta) > 1.7 && abs(eta) <= 3.2) * sqrt(energy^2*0.050^2 + energy*0.706^2) + \
			     (abs(eta) > 3.0 && abs(eta) <= 4.9) * sqrt(energy^2*0.05^2 + energy*1.00^2)
  }

}

####################################################################                                                                                                                    
## Track pile-up subtractor: apply CHS on top of track collection ##                                                                                                                    
####################################################################                                                                                                                     
module TrackPileUpSubtractor TrackPileUpSubtractor {
  ## take tracks from calorimeter module, smeared electrons and smeared muon. Take the PV from the ModifyBeamSpot --> pileup subtraction or CHS                                          
  ## Is not useful to run this module on top of NoPU collections                                                                                                                        
  add InputArray Calorimeter/eflowTracks eflowTracks
  add InputArray ElectronEnergySmearing/electrons electrons
  add InputArray MuonMomentumSmearing/muons muons
  set PVInputArray  ModifyBeamSpot/PV
  # assume perfect pile-up subtraction for tracks with |z| > fZVertexResolution in m                                                                                                    
  set ZVertexResolution 0.0001
}

########################                                                                                                                                                                 
## Energy flow merger ##                                                                                                                                                                 
########################                                                                                                                                                                 
module Merger EFlowMerger {
  ## charged particles after having applied CHS                                                                                                                                          
  add InputArray TrackPileUpSubtractor/eflowTracks
  ## calorimeter towers to get also photons and neutral hadrons for the whole detector coverage                                                                                         
  add InputArray Calorimeter/eflowTowers
  ## muons that are not included in the pileup subtractor                                                                                                                               
  add InputArray MuonMomentumSmearing/muons
  set OutputArray eflow
}

##############################################                                                                                                                                         
## Calculate Rho using KtClustering or grid ##                                                                                                                                   
##############################################                                                                                                                                           
module FastJetFinder GlobalRhoKt4 {
  ## take as input the particle flow particles                                                                                                                                            
  set InputArray EFlowMerger/eflow
  ## output name                                                                                                                                                                         
  set RhoOutputArray rho
  ## compute rho clustering the event                                                                                                                                                     
  set ComputeRho     true
  ## not compute rho using the grid                                                                                                                                                       
  set ComputeRhoGrid false
  ## area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area                                           
  set AreaAlgorithm 1
  ## jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt                                                                                               
  set JetAlgorithm 4
  ## Clustering and Ghost parameter                                                                                                                                                       
  set ParameterR  0.4
  set GhostEtaMax 5.0
  set RhoEtaMax   5.0
  ## eta bins for rho evaluation                                                                                                                                                          
  add RhoEtaRange 0.0 5.0
  set JetPTMin 0.0
}


module FastJetFinder GlobalRhoGridFastJet {
  ## take as input the particle flow particles                                                                                                                                             
  set InputArray EFlowMerger/eflow
  ## compute rho clustering the event                                                                                                                                                      
  set ComputeRho     false
  ## not compute rho using the grid                                                                                                                                                        
  set ComputeRhoGrid true
  set RhoOutputArray rho
  ## area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area                                              
  set AreaAlgorithm 1
  ## jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt                                                                                                 
  set JetAlgorithm 4
  ## Clustering and Ghost parameter                                                                                                                                                        
  set ParameterR  0.4
  set GhostEtaMax 5.0
  set RhoEtaMax   5.0
  ## eta bins for rho evaluation                                                                                                                                                           
  add RhoEtaRange 0.0 5.0
  set JetPTMin 0.0
}

#######################################                                                                                                                                                    
## Rho pile-up in different eta bins ##                                                                                                                                                    
#######################################                                                                                                                                                    
module FastJetFinder RhoKt4 {
  # input particles                                                                                                                                                                        
  set InputArray EFlowMerger/eflow
  # output name                                                                                                                                                                            
  set RhoOutputArray rho
  ## compute rho clustering the event                                                                                                                                                      
  set ComputeRho     true
  ## not compute rho using the grid                                                                                                                                                        
  set ComputeRhoGrid false
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area                                              
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt                                                                                                 
  set JetAlgorithm 4
  set ParameterR   0.4
  set GhostEtaMax  5.0
  set RhoEtaMax    5.0

  add RhoEtaRange 0.0 2.5
  add RhoEtaRange 2.5 4.0
  add RhoEtaRange 4.0 5.0

  set JetPTMin 0.0
}

module FastJetFinder RhoGridFastJet {
  ## take as input the particle flow particles                                                                                                                                             
  set InputArray EFlowMerger/eflow
  ## output name                                                                                                                                                                           
  set RhoOutputArray rho
  ## compute rho clustering the event                                                                                                                                                      
  set ComputeRho     false
  ## not compute rho using the grid                                                                                                                                                        
  set ComputeRhoGrid true
  ## area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area                                              
  set AreaAlgorithm 1
  ## Clustering and Ghost parameter                                                                                                                                                        
  set GhostEtaMax 5.0
  set RhoEtaMax   5.0
  ## eta bins for rho evaluation                                                                                                                                                           
  add RhoEtaRange 0.0 2.5
  add RhoEtaRange 2.5 4.0
  add RhoEtaRange 4.0 5.0

  set JetPTMin 0.0
}

################                                                                                                                                                                          
## Jet finder ##                                                                                                                                                                          
################                                                                                                                                                                           

module FastJetFinder FastJetFinder {
  set InputArray EFlowMerger/eflow
  set OutputArray jets
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area                                               
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt                                                                                                  
  set JetAlgorithm 6
  set ParameterR   0.4

  set JetPTMin 10.0
}

##################                                                                                                                                                                         
### Track jets ###                                                                                                                                                                         
##################                                                                                                                                                                         
module FastJetFinder TrackJetFinder {
  set InputArray  TrackPileUpSubtractor/eflowTracks
  set OutputArray jets
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area                                               
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt                                                                                                  
  set JetAlgorithm 6
  set ParameterR   0.4
  set JetPTMin      2.0
  set ParticlePTMin 0.35
}


#######################                                                                                                                                                                    
# MC truth jet finder #                                                                                                                                                                    
#######################                                                                                                                                                                    
module FastJetFinder GenJetFinder {
  ## generator level particle, no smearing and no efficiecny                                                                                                                               
  set InputArray Delphes/stableParticles
  set OutputArray jets
  set AreaAlgorithm 1
  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt                                                                                                      
  set JetAlgorithm 6
  set ParameterR   0.4
  set JetPTMin 5.0
}

############################################                                                                                                                                               
## Jet Pile-Up Subtraction: L1 correction ##                                                                                                                                               
############################################                                                                                                                                               

module JetPileUpSubtractor JetPileUpSubtractor { ## make the rho correction                                                                                                                
  ## input jets                                                                                                                                                                            
  set JetInputArray FastJetFinder/jets
  ## take the Rho from cluster the event with kt jets (decide to use this or the median grid or the safeAreaSbtraction)                                                                    
  set RhoInputArray RhoKt4/rho
  ## output jets                                                                                                                                                                           
  set OutputArray jets
  set doSafe4VAreaSubtraction false
  set JetPTMin 10.0
}

module JetPileUpSubtractor JetPileUpSubtractorGrid { ## make the rho correction                                                                                                            
  ## input jets                                                                                                                                                                            
  set JetInputArray FastJetFinder/jets
  ## take the Rho from cluster the event with kt jets (decide to use this or the median grid or the safeareasbtraction)                                                                    
  set RhoInputArray RhoGridFastJet/rho
  ## output jets                                                                                                                                                                           
  set OutputArray jets
  set doSafe4VAreaSubtraction false
  set JetPTMin 10.0
}

module JetPileUpSubtractor JetPileUpSubtractor4VArea { ## make the rho correction using safe 4V subtraction                                                                                 
  ## input jets                                                                                                                                                                            
  set JetInputArray FastJetFinder/jets
  ## not used when doSafe4VAreaSubtraction is true                                                                                                                                         
  set RhoInputArray RhoGridFastJet/rho
  ## output jets                                                                                                                                                                           
  set OutputArray jets
  set JetPTMin 10.0

  ## options for 4V safe subtracion                                                                                                                                                        
  set doSafe4VAreaSubtraction true
  ## use this info only if doSafe4VAreaSubtraction is set to true                                                                                                                          
  set InputArray EFlowMerger/eflow
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area                                               
  set AreaAlgorithmRho 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt                                                                                                  
  set JetAlgorithmRho 4
  set ParameterRRho   0.4
  set GhostEtaMaxRho  5.0
  set RhoEtaMaxRho    5.0
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area                                               
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt                                                                                                  
  set JetAlgorithm 6
  set ParameterR   0.4
  ## eta bins for rho evaluation                                                                                                                                                           
  add RhoEtaRange 0.0 2.5
  add RhoEtaRange 2.5 4.0
  add RhoEtaRange 4.0 5.0

}

module JetFlavourAssociation  JetFlavourAssociation {

  set PartonInputArray    Delphes/partons
  set ParticleInputArray  Delphes/allParticles
  set LHEPartonInputArray Delphes/LHEParticles
  set JetInputArray       JetPileUpSubtractor/jets

  set DeltaR 0.4
  set PartonPTMin 0.5
  set PartonEtaMax 4.0

}

################################################################################                                                                                                           
### Neutrino Filter on generated particles  of status 1 without any smearing ###                                                                                                           
################################################################################                                                                                                           
module NeutrinoFilter NeutrinoFilter {
  set InputArray  Delphes/stableParticles
  set OutputArray stableParticles
}

### make GenJets without neutrino                                                                                                                                                          
module FastJetFinder GenJetFinderNoNu {
  ## input particles                                                                                                                                                                       
  set InputArray NeutrinoFilter/stableParticles
  ## output name                                                                                                                                                                           
  set OutputArray jets
  set AreaAlgorithm 1
  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt                                                                                                      
  set JetAlgorithm 6
  set ParameterR 0.4
  set JetPTMin 10.0
}

### -sum of all particles after filtering neutrinos                                                                                                                                       
module Merger GenMissingET {
  add InputArray NeutrinoFilter/stableParticles
  set MomentumOutputArray momentum
}

###########################                                                                                                                                                                
### Run the puppi code  ###                                                                                                                                                                
###########################                                                                                                                                                                 
module RunPUPPI RunPUPPI {
  ## input information                                                                                                                                                                     
  set TrackInputArray   Calorimeter/eflowTracks
  set NeutralInputArray Calorimeter/eflowTowers
  set PVInputArray      ModifyBeamSpot/PV
  ## min puppi weight and use dZ vertex option                                                                                                                                             
  set MinPuppiWeight 0.01
  set UseExp         false
  ## define puppi algorithm parameters (more than one for the same eta region is possible)                                                                                                 
  add EtaMinBin           0.   2.5   3.0  2.5  3.0
  add EtaMaxBin           2.5  3.0   10.0 3.0  10.0
  add PtMinBin            0.1  0.5   1.0  0.5  1.0
  add ConeSizeBin         0.3  0.3   0.3  0.3  0.3
  add RMSPtMinBin         0.1  0.5   0.5  0.5  0.5
  add RMSScaleFactorBin   1.0  1.0   1.0  1.0  1.0
  add NeutralMinEBin      0.2  0.2   0.2  0.2  0.2
  add NeutralPtSlope      0.02 0.02  0.02 0.02 0.02
  add ApplyCHS            true true  true true true
  add UseCharged          true false false false false
  add ApplyLowPUCorr      true true  true  true  true
  add MetricId            5    0     0     1     1
  ## output name                                                                                                                                                                           
  set OutputArray PuppiParticles
  set OutputArrayTracks   puppiTracks
  set OutputArrayNeutrals puppiNeutrals
}

#####################                                                                                                                                                                      
## Make puppi jets ##                                                                                                                                                                      
#####################                                                                                                                                                                       

module FastJetFinder PuppiJetFinder {
  set InputArray RunPUPPI/PuppiParticles
  set OutputArray jets
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area                                               
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt                                                                                                  
  set JetAlgorithm 6
  set ParameterR   0.4
  set JetPTMin     10.
}

#####################################                                                                                                                                                      
## Maker rho corrections for Puppi ##                                                                                                                                                      
#####################################                                                                                                                                                       

module FastJetFinder PuppiRhoKt4 {
  set InputArray RunPUPPI/PuppiParticles
  ## compute rho clustering the event                                                                                                                                                      
  set ComputeRho     true
  ## not compute rho using the grid                                                                                                                                                        
  set ComputeRhoGrid false
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area                                               
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt                                                                                                  
  set JetAlgorithm 4
  set ParameterR   0.4
  set GhostEtaMax  5.0
  set RhoEtaMax    5.0
  add RhoEtaRange 0.0 2.5
  add RhoEtaRange 2.5 5.0
  set JetPTMin     0.0
}

module FastJetFinder PuppiRhoGridFastJet {
  ## take as input the particle flow particles                                                                                                                                             
  set InputArray RunPUPPI/PuppiParticles
  ## compute rho clustering the event                                                                                                                                                      
  set ComputeRho     false
  ## not compute rho using the grid                                                                                                                                                        
  set ComputeRhoGrid true
  set RhoOutputArray rho
  ## area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area                                              
  set AreaAlgorithm 1
  ## Clustering and Ghost parameter                                                                                                                                                        
  set GhostEtaMax 5.0
  set RhoEtaMax   5.0
  ## eta bins for rho evaluation                                                                                                                                                           
  add RhoEtaRange 0.0 2.5
  add RhoEtaRange 2.5 5.0
  set JetPTMin 0.0
}

########################                                                                                                                                                                   
## Correct puppi jets ##                                                                                                                                                                   
########################                                                                                                                                                                    
module JetPileUpSubtractor PuppiJetPileUpSubtractor { ## make the rho correction                                                                                                           
  set JetInputArray PuppiJetFinder/jets
  ## take the Rho from cluster the event with kt jets (decide to use this or the median grid or the safeAreaSbtraction)                                                                    
  set RhoInputArray PuppiRhoKt4/rho
  set OutputArray jets
  set doSafe4VAreaSubtraction false
  set JetPTMin 10.0
}



module JetPileUpSubtractor PuppiJetPileUpSubtractorGrid { ## make the rho correction                                                                                                       
  set JetInputArray PuppiJetFinder/jets
  ## take the Rho from cluster the event with kt jets (decide to use this or the median grid or the safeareasbtraction)                                                                    
  set RhoInputArray PuppiRhoGridFastJet/rho
  set OutputArray jets
  set doSafe4VAreaSubtraction false
  set JetPTMin 10.0
}

module JetPileUpSubtractor PuppiJetPileUpSubtractor4VArea { ## make the rho correction                                                                                                     
  set JetInputArray PuppiJetFinder/jets
  ## not used when doSafe4VAreaSubtraction is true                                                                                                                                         
  set RhoInputArray PuppiRhoGridFastJet/rho
  set OutputArray   jets
  set doSafe4VAreaSubtraction true
  set JetPTMin      10.0
  ## use this info only if doSafe4VAreaSubtraction is set to true                                                                                                                          
  set InputArray    RunPUPPI/PuppiParticles
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area                                               
  set AreaAlgorithmRho 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt                                                                                                  
  set JetAlgorithmRho 4
  set ParameterRRho   0.4
  set GhostEtaMaxRho  5.0
  set RhoEtaMaxRho    5.0

  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area                                               
  set AreaAlgorithm   1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt                                                                                                  
  set JetAlgorithm    6
  set ParameterR      0.4

 ## eta bins for rho evaluation                                                                                                                                                            
  add RhoEtaRange 0.0 2.5
  add RhoEtaRange 2.5 3.0
  add RhoEtaRange 3.0 10.0
}

#############################                                                                                                                                                            
## Jet Flavour Association ##                                                                                                                                                           
#############################                                                                                                                                                           

module JetFlavourAssociation  PuppiJetFlavourAssociation {

  set PartonInputArray    Delphes/partons
  set ParticleInputArray  Delphes/allParticles
  set LHEPartonInputArray Delphes/LHEParticles
  set JetInputArray       PuppiJetPileUpSubtractor/jets

  set DeltaR 0.4
  set PartonPTMin 0.5
  set PartonEtaMax 4.0

}


#####################                                                                                                                                                                      
# Missing ET merger #                                                                                                                                                                      #####################                                                                                                                                                                       
module Merger MissingET {
  add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
}


###############                                                                                                                                                                            
## PUPPI MET ##                                                                                                                                                                            
###############                                                                                                                                                                             
module Merger PuppiMissingET {
  add InputArray RunPUPPI/PuppiParticles
  set MomentumOutputArray momentum
}



###################
# Photon efficiency
###################

module Efficiency PhotonEfficiency {
  set InputArray Calorimeter/photons
  set OutputArray photons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
  # efficiency formula for photons
  set EfficiencyFormula {                                        (pt <= 10.0) * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 10.0)  * (0.9635) + \
                         (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 10.0)  * (0.9624) + \
                         (abs(eta) > 4.0)                                   * (0.00)}
}

####################                                                                                                                                                                       
# Photon isolation #                                                                                                                                                                       
####################                                                                                                                                                                        
module Isolation PhotonIsolation {
  # particle for which calculate the isolation                                                                                                                                             
  set CandidateInputArray        PhotonEfficiency/photons
  # neutral and charged particles for the whole event (no CHS applied)                                                                                                                     
  set NeutralIsolationInputArray Calorimeter/eflowTowers
  set ChargedIsolationInputArray Calorimeter/eflowTracks
  # select a rho for the isolation                                                                                                                                                         
  set RhoInputArray RhoKt4/rho
  # output array                                                                                                                                                                           
  set OutputArray photons
  # isolation cone                                                                                                                                                                         
  set DeltaRMax 0.3
  # minimum pT                                                                                                                                                                             
  set PTMin 0.5
  # iso ratio to cut                                                                                                                                                                       
  set PTRatioMax 9999.
}

#####################
# Electron efficiency
#####################

module Efficiency ElectronEfficiency {
  set InputArray TrackPileUpSubtractor/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for electrons
    set EfficiencyFormula {                                      (pt <= 4.0)  * (0.00) + \
                         (abs(eta) <= 1.45 ) * (pt >  4.0 && pt <= 6.0)   * (0.50) + \
                         (abs(eta) <= 1.45 ) * (pt >  6.0 && pt <= 8.0)   * (0.70) + \
                         (abs(eta) <= 1.45 ) * (pt >  8.0 && pt <= 10.0)  * (0.85) + \
                         (abs(eta) <= 1.45 ) * (pt > 10.0 && pt <= 30.0)  * (0.94) + \                                                      
                         (abs(eta) <= 1.45 ) * (pt > 30.0 && pt <= 50.0)  * (0.97) + \                          
                         (abs(eta) <= 1.45 ) * (pt > 50.0 && pt <= 70.0)  * (0.98) + \          
                         (abs(eta) <= 1.45 ) * (pt > 70.0 )  * (1.0) + \                                                                                                                 
                         (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt >  4.0 && pt <= 10.0)   * (0.35) + \
                         (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt > 10.0 && pt <= 30.0)   * (0.40) + \   
                         (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt > 30.0 && pt <= 70.0)   * (0.45) + \                                 
                         (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt > 70.0 )  * (0.55) + \    
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt >  4.0 && pt <= 10.0)  * (0.75) + \
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 10.0 && pt <= 30.0)  * (0.85) + \                                                      
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 30.0 && pt <= 50.0)  * (0.95) + \                          
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 50.0 && pt <= 70.0)  * (0.95) + \          
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 70.0 )  * (1.0) + \   
                         (abs(eta) >= 2.0 && abs(eta) <= 4.0 ) * (pt >  4.0 && pt <= 10.0)  * (0.65) + \
                         (abs(eta) >= 2.0 && abs(eta) <= 4.0 ) * (pt > 10.0 && pt <= 30.0)  * (0.75) + \                                                      
                         (abs(eta) >= 2.0 && abs(eta) <= 4.0 ) * (pt > 30.0 && pt <= 50.0)  * (0.90) + \                          
                         (abs(eta) >= 2.0 && abs(eta) <= 4.0 ) * (pt > 50.0 && pt <= 70.0)  * (0.90) + \          
                         (abs(eta) >= 2.0 && abs(eta) <= 4.0 ) * (pt > 70.0 )  * (0.90) + \                                                                                               
                         (abs(eta) > 4.0)                              * (0.00)}

}

######################                                                                                                                                                                     
# Electron isolation #                                                                                                                                                                     
######################                                                                                                                                                                      
module Isolation ElectronIsolation {
  set CandidateInputArray        ElectronEfficiency/electrons
  set NeutralIsolationInputArray Calorimeter/eflowTowers
  set ChargedIsolationInputArray Calorimeter/eflowTracks
  set RhoInputArray RhoKt4/rho
  set OutputArray electrons
  set DeltaRMax 0.3
  set PTMin 0.5
  set PTRatioMax 9999.
}


#################
# Muon efficiency
#################

module Efficiency MuonEfficiency {
  set InputArray TrackPileUpSubtractor/muons
  set OutputArray muons

  # set EfficiencyFormula {efficiency as a function of eta and pt}
  # efficiency formula for muons
    set EfficiencyFormula {                                    (pt <= 2.0)  * (0.00) + \  
                         (abs(eta) <= 4.00) * (pt >  2.0 && pt <= 3.0)  * (0.51) + \
                         (abs(eta) <= 4.00) * (pt >  3.0 && pt <= 4.0)  * (0.85) + \ 
                         (abs(eta) <= 4.00) * (pt >  4.0 && pt <= 11.0) * (0.93) + \               
                         (abs(eta) <= 4.00) * (pt >  11. && pt <= 50.)  * (0.96) + \   
                         (abs(eta) <= 4.00) * (pt >  50. && pt <= 70.)  * (0.98) + \                      
                         (abs(eta) <= 4.00) * (pt > 70.0 )  * (1.00) + \   
                         (abs(eta) > 4.00)  * (0.00)}

}

##################                                                                                                                                                                         
# Muon isolation #                                                                                                                                                                         
##################                                                                                                                                                                          
module Isolation MuonIsolation {
  set CandidateInputArray MuonEfficiency/muons
  set NeutralIsolationInputArray Calorimeter/eflowTowers
  set ChargedIsolationInputArray Calorimeter/eflowTracks
  set RhoInputArray RhoKt4/rho
  set OutputArray muons
  set DeltaRMax 0.3
  set PTMin 0.5
  set PTRatioMax 9999.
}


##########################                                                                                                                                                                 
# Hadronicn Tau tagging ##                                                                                                                                                                 
##########################                                                                                                                                                                  
module TauTagging TauTagging {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray   Delphes/partons
  set JetInputArray      JetPileUpSubtractor/jets
  set DeltaR 0.4
  set TauPTMin 1.0
  set TauEtaMax 4.0
    # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}                                                                                               
  # default efficiency formula (misidentification rate)                                                                                                                                    
  add EfficiencyFormula {0} {(abs(eta)<1.8)*0.006+(abs(eta)>1.8)*0.015}
  # efficiency formula for tau-jets                                                                                                                                                        
  add EfficiencyFormula {15} {(abs(eta)<2.4)*(0.65)+(abs(eta)>2.4)*(0.65)}
}


#####################
#### BTagging #######
#####################
module BTagging BTagging {

  set JetInputArray JetPileUpSubtractor/jets

  add EfficiencyFormulaLoose {0} {0.02}

  add EfficiencyFormulaLoose {4} { (pt <= 15.0) * (0.000) + \
                                   (abs(eta) <= 1.2) * (pt > 15.0) * (0.29*tanh(pt*0.0183 - 0.2196)) + \
                                   (abs(eta) > 1.2 && abs(eta) <= 4.0) * (pt > 15.0) * (0.29*tanh(pt*0.00997 - 0.143)) + \
                                   (abs(eta) > 4.0) * (0.000)
  }

  add EfficiencyFormulaLoose {5} { (pt <= 15.0) * (0.000) + \
                                   (abs(eta) <= 1.2) * (pt > 15.0) * (0.629858*tanh(pt*0.0166188 + 0.300119)) + \
                                   (abs(eta) > 1.2 && abs(eta) <= 4.0) * (pt > 15.0) * (0.584522*tanh(pt*0.0144387 + 0.397034)) + \
                                   (abs(eta) > 4.0) * (0.000)
  }

  add EfficiencyFormulaMedium {0} {0.001}
  # efficiency formula for c-jets (misidentification rate)                                                                                                                                  
  add EfficiencyFormulaMedium {4} { (pt <= 15.0) * (0.000) + \
                                    (abs(eta) <= 1.2) * (pt > 15.0) * (0.1873*tanh(pt*0.0183 - 0.2196)) + \
                                    (abs(eta) > 1.2 && abs(eta) <= 4.0) * (pt > 15.0) * (0.1898*tanh(pt*0.00997 - 0.143)) + \
                                    (abs(eta) > 4.0) * (0.000)
  }
  # efficiency formula for b-jets                                                                                                                                                           
  add EfficiencyFormulaMedium {5} { (pt <= 15.0) * (0.000) + \
                                    (abs(eta) <= 1.2) * (pt > 15.0) * (0.629858*tanh(pt*0.0166188 + 0.300119)) + \
                                    (abs(eta) > 1.2 && abs(eta) <= 4.0) * (pt > 15.0) * (0.584522*tanh(pt*0.0144387 + 0.397034)) + \
                                    (abs(eta) > 4.0) * (0.000)
  }

}

module BTagging PuppiBTagging {

  set JetInputArray PuppiJetPileUpSubtractor/jets

  add EfficiencyFormulaLoose {0} {0.02}

  add EfficiencyFormulaLoose {4} { (pt <= 15.0) * (0.000) + \
                                   (abs(eta) <= 1.2) * (pt > 15.0) * (0.29*tanh(pt*0.0183 - 0.2196)) + \
                                   (abs(eta) > 1.2 && abs(eta) <= 4.0) * (pt > 15.0) * (0.29*tanh(pt*0.00997 - 0.143)) + \
                                   (abs(eta) > 4.0) * (0.000)
  }

  add EfficiencyFormulaLoose {5} { (pt <= 15.0) * (0.000) + \
                                   (abs(eta) <= 1.2) * (pt > 15.0) * (0.629858*tanh(pt*0.0166188 + 0.300119)) + \
                                   (abs(eta) > 1.2 && abs(eta) <= 4.0) * (pt > 15.0) * (0.584522*tanh(pt*0.0144387 + 0.397034)) + \
                                   (abs(eta) > 4.0) * (0.000)
  }

  add EfficiencyFormulaMedium {0} {0.001}
  # efficiency formula for c-jets (misidentification rate)                                                                                                                                  
  add EfficiencyFormulaMedium {4} { (pt <= 15.0) * (0.000) + \
                                    (abs(eta) <= 1.2) * (pt > 15.0) * (0.1873*tanh(pt*0.0183 - 0.2196)) + \
                                    (abs(eta) > 1.2 && abs(eta) <= 4.0) * (pt > 15.0) * (0.1898*tanh(pt*0.00997 - 0.143)) + \
                                    (abs(eta) > 4.0) * (0.000)
  }
  # efficiency formula for b-jets                                                                                                                                                           
  add EfficiencyFormulaMedium {5} { (pt <= 15.0) * (0.000) + \
                                    (abs(eta) <= 1.2) * (pt > 15.0) * (0.629858*tanh(pt*0.0166188 + 0.300119)) + \
                                    (abs(eta) > 1.2 && abs(eta) <= 4.0) * (pt > 15.0) * (0.584522*tanh(pt*0.0144387 + 0.397034)) + \
                                    (abs(eta) > 4.0) * (0.000)
  }
}


###################                                                                                                                                                                      
## PileUp Jet ID ##                                                                                                                                                                        ###################                                                                                                                                                                     

module PileUpJetID PileUpJetID {
  set JetInputArray     JetPileUpSubtractor/jets
  set TrackInputArray   Calorimeter/eflowTracks
  set NeutralInputArray Calorimeter/eflowTowers
  set PVInputArray      ModifyBeamSpot/PV

  set OutputArray   jets

  set ParameterR 0.4
  set JetPTMin   10.0
  set UseConstituents 0

  add Cones  0.1 0.2 0.3 0.4 0.5 0.6 0.7

  add Pt010_Tight_betaStar   0.15 0.15 999. 999.
  add Pt1020_Tight_betaStar  0.15 0.15 999. 999.
  add Pt2030_Tight_betaStar  0.15 0.15 999. 999.
  add Pt3050_Tight_betaStar  0.15 0.15 999. 999.

  add Pt010_Tight_RMS  0.06 0.07 0.04 0.05
  add Pt1020_Tight_RMS 0.06 0.07 0.04 0.05
  add Pt2030_Tight_RMS 0.05 0.07 0.03 0.045
  add Pt3050_Tight_RMS 0.05 0.06 0.03 0.04

  add Pt010_Medium_betaStar  0.2 0.3 999. 999.
  add Pt1020_Medium_betaStar 0.2 0.3 999. 999.
  add Pt2030_Medium_betaStar 0.2 0.3 999. 999.
  add Pt3050_Medium_betaStar 0.2 0.3 999. 999.

  add Pt010_Medium_betaStar  0.2 0.3 999. 999.
  add Pt1020_Medium_betaStar 0.2 0.3 999. 999.
  add Pt2030_Medium_betaStar 0.2 0.3 999. 999.
  add Pt3050_Medium_betaStar 0.2 0.3 999. 999.

  add Pt010_Medium_RMS    0.06 0.03 0.03 0.04
  add Pt1020_Medium_RMS   0.06 0.03 0.03 0.04
  add Pt2030_Medium_RMS   0.06 0.03 0.03 0.04
  add Pt3050_Medium_RMS   0.06 0.03 0.03 0.04

  add Pt010_Loose_betaStar  0.2 0.3 999. 999
  add Pt1020_Loose_betaStar 0.2 0.3 999. 999
  add Pt2030_Loose_betaStar 0.2 0.3 999. 999
  add Pt3050_Loose_betaStar 0.2 0.3 999. 999

  add Pt010_Loose_RMS   0.06 0.05 0.05 0.07
  add Pt1020_Loose_RMS  0.06 0.05 0.05 0.07
  add Pt2030_Loose_RMS  0.06 0.05 0.05 0.055
  add Pt3050_Loose_RMS  0.06 0.05 0.05 0.055

}

module PileUpJetID PuppiPileUpJetID {
  set JetInputArray     PuppiJetPileUpSubtractor/jets
  set TrackInputArray   RunPUPPI/puppiTracks
  set NeutralInputArray RunPUPPI/puppiNeutrals
  set PVInputArray      ModifyBeamSpot/PV

  set OutputArray   jets

  set ParameterR 0.4
  set JetPTMin   10.0
  set UseConstituents 1

  add Cones  0.1 0.2 0.3 0.4 0.5 0.6 0.7

  add Pt010_Tight_betaStar   0.15 0.15 999. 999.
  add Pt1020_Tight_betaStar  0.15 0.15 999. 999.
  add Pt2030_Tight_betaStar  0.15 0.15 999. 999.
  add Pt3050_Tight_betaStar  0.15 0.15 999. 999.

  add Pt010_Tight_RMS  0.06 0.07 0.04 0.04
  add Pt1020_Tight_RMS 0.06 0.07 0.04 0.04
  add Pt2030_Tight_RMS 0.05 0.07 0.03 0.04
  add Pt3050_Tight_RMS 0.05 0.06 0.03 0.04

  add Pt010_Medium_betaStar  0.2 0.3 999. 999.
  add Pt1020_Medium_betaStar 0.2 0.3 999. 999.
  add Pt2030_Medium_betaStar 0.2 0.3 999. 999.
  add Pt3050_Medium_betaStar 0.2 0.3 999. 999.

  add Pt010_Medium_RMS    0.06 0.03 0.03 0.05
  add Pt1020_Medium_RMS   0.06 0.03 0.03 0.05
  add Pt2030_Medium_RMS   0.06 0.03 0.03 0.045
  add Pt3050_Medium_RMS   0.06 0.03 0.03 0.04

  add Pt010_Loose_betaStar  0.2 0.3 999. 999
  add Pt1020_Loose_betaStar 0.2 0.3 999. 999
  add Pt2030_Loose_betaStar 0.2 0.3 999. 999
  add Pt3050_Loose_betaStar 0.2 0.3 999. 999

  add Pt010_Loose_RMS   0.06 0.05 0.05 0.07
  add Pt1020_Loose_RMS  0.06 0.05 0.05 0.07
  add Pt2030_Loose_RMS  0.06 0.05 0.05 0.055
  add Pt3050_Loose_RMS  0.06 0.05 0.05 0.055

}

########################                                                                                                                                                                   ## Constituent filter ##                                                                                                                                                                   ########################                                                                                                                                                                  

module ConstituentFilter ConstituentFilter {

  set ConEMin 0.

  add JetInputArray GenJetFinderNoNu/jets
  add JetInputArray PileUpJetID/jets

  add ConstituentInputArray Delphes/stableParticles stableParticles
  add ConstituentInputArray TrackPileUpSubtractor/eflowTracks eflowTracks
  add ConstituentInputArray Calorimeter/eflowTowers eflowTowers
  add ConstituentInputArray MuonMomentumSmearing/muons muons

}

module ConstituentFilter ConstituentFilterPUPPI {

  set ConEMin 0.

  add JetInputArray GenJetFinderNoNu/jets
  add JetInputArray PileUpJetIDPUPPI/jets

  add ConstituentInputArray Delphes/stableParticles stableParticles
  add ConstituentInputArray TrackPileUpSubtractor/eflowTracks eflowTracks
  add ConstituentInputArray Calorimeter/eflowTowers eflowTowers
  add ConstituentInputArray MuonMomentumSmearing/muons muons

}


##################                                                                                                                                                                         
# ROOT tree writer                                                                                                                                                                         
##################                                                                                                                                                          

module TreeWriter TreeWriter {
  ## branch notation : <particle collection> <branch name> <type of object in classes/DelphesClass.h<                                                                                     
 
  ## input status 1 particle from Pythia8                                                                                                                                                  
  #add Branch Delphes/stableParticles GenParticles GenParticle                                                                                                                            
 
  ## input partons                                                                                                                                                                         
  #add Branch Delphes/partons GenParton GenParticle                                                                                                                                       
 
  ## LHE particles                                                                                                                                                                         
  add Branch Delphes/LHEParticles LHEParticles LHEParticle

  ## NPU after Pileup Merging                                                                                                                                                              
  add Branch PileUpMerger/NPU NPU ScalarHT

  ## gen particles after vertex smearing                                                                                                                                                   
  #add Branch GenBeamSpotFilter/beamSpotParticles GenBeamSpotParticles GenParticle                                                                                                
 
  ## particle after B field propagation                                                                                                                                                    
  #add Branch ParticlePropagator/stableParticles particlePropagator GenParticle                                                                                                            
  #add Branch ParticlePropagator/electrons       electronPropagator GenParticle                                                                                                            
  #add Branch ParticlePropagator/muons           muonPropagator GenParticle                                                                                                                

  ## after Pt filter: all delphes particles, not only status 1                                                                                                                             
  add Branch StatusPid/filteredParticles GenParticles GenParticle                                                                                                                  
 
  ## track collection after: charged hadrons smearing and track eff, electron smearing and track eff                                                                                       
  #add Branch TrackMerger/tracks trackCollectionNoMU Track                                                                                                                          
 
  ## output of the calorimeter simulation                                                                                                                                                  
  #add Branch Calorimeter/towers caloTowers Tower                                                                                                                                          
  #add Branch Calorimeter/photons RawPhotons Photon                                                                                                                                        
  #add Branch Calorimeter/eflowTracks eflowTracks Track                                                                                                                                    
  #add Branch Calorimeter/eflowTracks eflowTowers Tower                                                                                                                                
 
  ## tracks after CHS                                                                                                                                                                      
  #add Branch TrackPileUpSubtractor/eflowTracks trackCollectionCHS Track                                                                                                               
 
  ## eflow output                                                                                                                                                                          
  #add Branch EFlowMerger/eflow eflowCandidates Track                                                                                                                            
  
  ## Rho values                                                                                                                                                                            
  add Branch GlobalRhoKt4/rho GlobalRhoKt4 Rho
  add Branch GlobalRhoGridFastJet/rho GlobalRhoGridFastJet Rho
  add Branch RhoKt4/rho RhoKt4 Rho
  add Branch RhoGridFastJet/rho RhoGridFastJet Rho

  ## Standard Jets                                                                                                                                                                         
  #add Branch FastJetFinder/jets RawJet Jet                                                                                                                                               
  #add Branch GenJetFinder/jets GenJetWithNu Jet                                                                                                                                           
  add Branch GenJetFinderNoNu/jets GenJet Jet
  #add Branch JetPileUpSubtractor/jets Jet Jet                                                                                                                                             
  add Branch TrackJetFinder/jets TrackJet Jet
  #add Branch JetPileUpSubtractorGrid/jets Jet Jet                                                                                                                                         
  #add Branch JetPileUpSubtractor4VArea/jets Jet4VArea Jet                                                                                                                                 
  add Branch PileUpJetID/jets JetPUID Jet

  ## PUPPI                                                                                                                                                                                 
  #add Branch RunPUPPI/PuppiParticles puppiParticles GenParticle
  add Branch PuppiRhoKt4/rho         PuppiRhoKt4 Rho
  add Branch PuppiRhoGridFastJet/rho PuppiRhoGridFastJet Rho
  #add Branch PuppiJetFinder/jets     RawPuppiJet Jet                                                                                                                                      
  #add Branch PuppiJetPileUpSubtractor/jets PuppiJet Jet                                                                                                                                   
  #add Branch PuppiJetPileUpSubtractorGrid/jets PuppiJetGrid Jet                                                                                                                           
  #add Branch PuppiJetPileUpSubtractor4VArea/jets PuppiJet4VArea Jet                                                                                                                       
  add Branch PuppiPileUpJetID/jets PuppiJetPUID Jet

   ## MET                                                                                                                                                                          
  add Branch GenMissingET/momentum GenMissingET MissingET
  add Branch MissingET/momentum MissingET MissingET
  add Branch PuppiMissingET/momentum PuppiMissingET MissingET
  ## photons and leptons                                                                                                                                                                   
  add Branch ElectronIsolation/electrons Electron Electron
  add Branch PhotonIsolation/photons Photon Photon
  add Branch MuonIsolation/muons Muon Muon

  set fOffsetFromModifyBeamSpot 0
}

#####################################################                                                                                                                                      
# Find uniquely identified photons/electrons/tau/jets                                                                                                                                      
#####################################################                                                                                                                                 
#module UniqueObjectFinder UniqueObjectFinderGJ {                                                                                                                                          
#   add InputArray PhotonIsolation/photons photons                                                                                                                                         
#   add InputArray JetPileUpSubtractor/jets jets                                                                                                                                           
#}                                                                                                                                                                                   
#module UniqueObjectFinder UniqueObjectFinderEJ {                                                                                                                                          
#   add InputArray ElectronIsolation/electrons electrons                                                                                                                                   
#   add InputArray UniqueObjectFinderGJ/jets jets                                                                                                                                          
#}                                                                                                                                                                                      
#module UniqueObjectFinder UniqueObjectFinderMJ {                                                                                                                                          
#   add InputArray MuonIsolation/muons muons                                                                                                                                               
#   add InputArray UniqueObjectFinderEJ/jets jets                                                                                                                                          
#}                                                                                                                                                                                       
