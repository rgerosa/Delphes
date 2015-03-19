#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {

  PileUpMerger 

  ModifyBeamSpot

  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  ChargedHadronMomentumSmearing
  ElectronEnergySmearing
  MuonMomentumSmearing

  TrackMerger

  Calorimeter

  TrackMergerWithMuon
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

  RunPUPPI
  PuppiRhoKt4
  PuppiRhoGrid
  PuppiJetFinder

  PuppiJetPileUpSubtractor
  PuppiJetFlavourAssociation

  PuppiBTagging

  PuppiPileUpJetID

  PhotonEfficiency
  PhotonIsolation 

  ElectronEfficiency 
  ElectronIsolation 
 
  MuonEfficiency
  MuonIsolation  

  GenMissingET
  MissingET
  PuppiMissingET

  GenScalarHT
  ScalarHT
  PuppiScalarHT

  TreeWriter

}

### remove some modules

# GenBeamSpotFilter
# StatusPid
# JetPileUpSubtractorGrid
# JetPileUpSubtractor4VArea
# PuppiJetPileUpSubtractorGrid
# PuppiJetPileUpSubtractor4VArea

#### remove the module which do the filter of jet constituent
# ConstituentFilter
# PuppiConstituentFilter

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
 set PileUpFile MB_1.mb
 #average expected pile up <poissonian generation>
 set MeanPileUp 140
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
  ## particles after propagation
  set InputArray  ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons
  # tracking efficiency formula for charged hadrons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) + \
                                           (abs(eta) <= 1.0) * (pt > 0.1   && pt <= 1.0)   * (0.71) + \
                                           (abs(eta) <= 1.0) * (pt > 1.0)                  * (0.81) + \
                         (abs(eta) > 1.0 && abs(eta) <= 1.8) * (pt > 0.1   && pt <= 1.0)   * (0.51) + \
                         (abs(eta) > 1.0 && abs(eta) <= 1.8) * (pt > 1.0)                  * (0.58) + \
			 (abs(eta) > 1.8 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.62) + \
			 (abs(eta) > 1.8 && abs(eta) <= 2.5) * (pt > 1.0)                  * (0.71) + \    
                         (abs(eta) > 2.5)                                                  * (0.00)
  }

}


#####################################
# Electron tracking efficiency - ID
####################################

module Efficiency ElectronTrackingEfficiency {
  set InputArray  ParticlePropagator/electrons
  set OutputArray electrons
  # tracking efficiency formula for electrons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.85) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e2) * (0.97) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0e2)                * (0.99) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.85) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e2) * (0.90) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e2)                * (0.95) + \
			 (abs(eta) > 2.5)                                                  * (0.00)

  }
}

##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons
  # tracking efficiency formula for muons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.998) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0)                  * (0.9998) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.98) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0)                  * (0.98) + \
			 (abs(eta) > 2.5)                                                  * (0.00)
  }
    
}

########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  ## hadrons after having applied the tracking efficiency
  set InputArray  ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons
  # resolution formula for charged hadrons
  set ResolutionFormula {                  (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e1) * (0.013) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.02) + \
                                           (abs(eta) <= 1.5) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e1) * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.04) + \
			 (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 2.0e2)                * (0.05)
  }
}

#################################
# Energy resolution for electrons
#################################

module EnergySmearing ElectronEnergySmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons
  # set ResolutionFormula {resolution formula as a function of eta and energy}
  set ResolutionFormula {   (abs(eta) <= 2.5) * (energy > 0.1   && energy <= 2.5e1) * (energy*0.015) + \
  			    (abs(eta) <= 2.5) * (energy > 2.5e1)                    * sqrt(energy^2*0.005^2 + energy*0.027^2 + 0.15^2) + \
                            (abs(eta) > 2.5 && abs(eta) <= 3.0)                                       * sqrt(energy^2*0.005^2 + energy*0.027^2 + 0.15^2) + \
 			    (abs(eta) > 3.0 && abs(eta) <= 5.0)                                       * sqrt(energy^2*0.08^2 + energy*1.97^2)
  } 

}

###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons
  # resolution formula for muons
  set ResolutionFormula {                  (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e1) * (0.012) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.015) + \
                                           (abs(eta) <= 1.5) * (pt > 2.0e2)                * (0.03) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e1) * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.025) + \
			 (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 2.0e2)                * (0.03)
  }

}

#################################################################
# Track merger : merge two collection of object in a output one #
#################################################################

module Merger TrackMerger {
  ## take smeared charged hadron and electrons as only tracks to take into account in the calorimeter simulation
  add InputArray ChargedHadronMomentumSmearing/chargedHadrons
  add InputArray ElectronEnergySmearing/electrons
  set OutputArray tracks
}

##############################################################
# Calorimeter : emulate calorimiter answer making caloTowers #
##############################################################

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

  ### energy deposition for each particle type
  # default energy fractions {abs(PDG code)} {Fecal Fhcal}
  add EnergyFraction {0} {0.0 1.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11}  {1.0 0.0}
  add EnergyFraction {22}  {1.0 0.0}
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
  set ECalResolutionFormula { (abs(eta) <= 1.497) * sqrt(energy^2*0.007^2 + energy*0.029^2 + 1.01^2) + \
                              (abs(eta) > 1.497 && abs(eta)<=3.0) * sqrt(energy^2*0.02^2 + energy*0.087^2 + 1.95^2) + \
			      (abs(eta) > 3.0 && abs(eta) <= 5.0) * sqrt(energy^2*0.29^2 + energy*0.86^2 + 191^2)
  }


  # set HCalResolutionFormula {resolution formula as a function of eta and energy}
  set HCalResolutionFormula {  (abs(eta) <= 1.7) * (energy*0.132 - sqrt(energy)*0.285 + 10) + \
                               (abs(eta) > 1.7 && abs(eta)<=2.1) * (energy*0.0737 - sqrt(energy)*0.0343 + 7.3) + \
                               (abs(eta) > 2.1 && abs(eta)<=2.3) * (energy*0.239 + sqrt(energy)*1.95 + 19.1) + \
                               (abs(eta) > 2.3 && abs(eta) <= 5.0) * (energy*0.0732 + sqrt(energy)*14.7+42.8)
  }

}

####################################################################
## Track pile-up subtractor: apply CHS on top of track collection ##
####################################################################

module TrackPileUpSubtractor TrackPileUpSubtractor {
  ## take tracks from calorimeter module, smeared electrons and smeared muon. Take the PV from the ModifyBeamSpot --> pileup subtraction or CHS
  ## Is not useful to run this module on top of NoPU collections
  add InputArray Calorimeter/eflowTracks    eflowTracks
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
  set InputArray   EFlowMerger/eflow
  set OutputArray  jets
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm  6
  set ParameterR    0.4
  set JetPTMin      10.0
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
  set JetPTMin      1.0
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
  set JetPTMin     15.0
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
  set JetPTMin 20.0
}

module JetPileUpSubtractor JetPileUpSubtractorGrid { ## make the rho correction 
  ## input jets
  set JetInputArray FastJetFinder/jets
  ## take the Rho from cluster the event with kt jets (decide to use this or the median grid or the safeareasbtraction)
  set RhoInputArray RhoGridFastJet/rho
  ## output jets
  set OutputArray jets
  set doSafe4VAreaSubtraction false
  set JetPTMin 20.0
}

module JetPileUpSubtractor JetPileUpSubtractor4VArea { ## make the rho correction using safe 4V subtraction
  ## input jets
  set JetInputArray FastJetFinder/jets
  ## not used when doSafe4VAreaSubtraction is true
  set RhoInputArray RhoGridFastJet/rho
  ## output jets
  set OutputArray jets
  set JetPTMin 20.0

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
  set DeltaR        0.4
  set PartonPTMin   0.5
  set PartonEtaMax  2.5
 
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
  set JetPTMin   15.0
}

### -sum of all particles after filtering neutrinos
module Merger GenMissingET {
  add InputArray NeutrinoFilter/stableParticles
  set MomentumOutputArray momentum
}

###########################
### Run the puppi code  ###
###########################
module Merger TrackMergerWithMuon {
  ## take smeared charged hadron and electrons as only tracks to take into account in the calorimeter simulation
  add InputArray Calorimeter/eflowTracks    
  add InputArray MuonMomentumSmearing/muons 
  set OutputArray tracks
}


module RunPUPPI RunPUPPI {
  ## input information
  set TrackInputArray   TrackMergerWithMuon/tracks
  set NeutralInputArray Calorimeter/eflowTowers
  set PVInputArray      ModifyBeamSpot/PV
  ## min puppi weight and use dZ vertex option
  set MinPuppiWeight    0.10
  set UseExp            false
  ## define puppi algorithm parameters (more than one for the same eta region is possible) 
  add EtaMinBin           0.    2.5    2.5    3.0   3.0      
  add EtaMaxBin           2.5   3.0    3.0    10.0  10.0
  add PtMinBin            0.    0.5    0.5    0.5   0.5    
  add ConeSizeBin         0.2   0.2    0.2    0.2   0.2
  add RMSPtMinBin         0.1   0.5    0.5    0.5   0.5
  add RMSScaleFactorBin   1.0   1.0    1.0    1.0   1.0
  add NeutralMinEBin      0.2   1.0    1.0    1.5   1.5
  add NeutralPtSlope      0.02  0.02   0.02   0.02  0.02
  add ApplyCHS            true  true   true   true  true
  add UseCharged          true  false  false  false false
  add ApplyLowPUCorr      true  true   true   true  true
  add MetricId            5     5      0      5     0
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
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 4
  set ParameterR   0.4
  set GhostEtaMax  5.0
  set RhoEtaMax    5.0  
  add RhoEtaRange  0.0 2.5
  add RhoEtaRange  2.5 4.0
  add RhoEtaRange  4.0 5.0
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
  add RhoEtaRange 2.5 4.0
  add RhoEtaRange 4.0 5.0
  set JetPTMin 0.0
}

module FastJetFinder PuppiRhoGrid {
  ## take as input the particle flow particles
  set InputArray RunPUPPI/PuppiParticles
  ## compute rho clustering the event  
  set ComputeRho     false
  ## not compute rho using the grid
  set ComputeRhoGridParticles true
  set RhoOutputArray rho
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
  set RhoInputArray PuppiRhoGrid/rho
  set OutputArray jets
  set doSafe4VAreaSubtraction false
  set JetPTMin 20.0
}



module JetPileUpSubtractor PuppiJetPileUpSubtractorGrid { ## make the rho correction 
  set JetInputArray PuppiJetFinder/jets
  ## take the Rho from cluster the event with kt jets (decide to use this or the median grid or the safeareasbtraction)
  set RhoInputArray PuppiRhoGridFastJet/rho
  set OutputArray jets
  set doSafe4VAreaSubtraction false
  set JetPTMin 20.0
}

module JetPileUpSubtractor PuppiJetPileUpSubtractor4VArea { ## make the rho correction 
  set JetInputArray PuppiJetFinder/jets
  ## not used when doSafe4VAreaSubtraction is true
  set RhoInputArray PuppiRhoGridFastJet/rho
  set OutputArray   jets
  set doSafe4VAreaSubtraction true
  set JetPTMin      20.0
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

  set DeltaR       0.4
  set PartonPTMin  0.5
  set PartonEtaMax 2.5
 
}

#####################
# Missing ET merger #
#####################

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

#####################
# Photon efficiency #
#####################

module Efficiency PhotonEfficiency {
  ## input particles
  set InputArray Calorimeter/photons
  ## output particles
  set OutputArray photons
  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
  # efficiency formula for photons
  set EfficiencyFormula {  (pt <= 20.0) * (0.00) + \
                           (abs(eta) <= 1.5) * (pt > 20 && pt <= 30) * (0.24) + \
                           (abs(eta) <= 1.5) * (pt > 30 && pt <= 40) * (0.51) + \
                           (abs(eta) <= 1.5) * (pt > 40 && pt <= 50) * (0.73) + \
		           (abs(eta) <= 1.5) * (pt > 50 && pt <= 60) * (0.86) + \
		           (abs(eta) <= 1.5) * (pt > 60 && pt <= 70) * (0.91) + \
		           (abs(eta) <= 1.5) * (pt > 70 && pt <= 80) * (0.94) + \
		           (abs(eta) <= 1.5) * (pt > 80 && pt <= 90) * (0.96) + \
		           (abs(eta) <= 1.5) * (pt > 90 && pt <= 100) * (0.97) + \
		           (abs(eta) <= 1.5) * (pt > 100 && pt <= 110) * (0.97) + \
		           (abs(eta) <= 1.5) * (pt > 110 && pt <= 120) * (0.98) + \
		           (abs(eta) > 1.5 && abs(eta)<=2.5) * (pt > 20 && pt <= 30) * (0.14) + \
		           (abs(eta) > 1.5 && abs(eta)<=2.5) * (pt > 30 && pt <= 40) * (0.32) + \
		           (abs(eta) > 1.5 && abs(eta)<=2.5) * (pt > 40 && pt <= 50) * (0.50) + \
		           (abs(eta) > 1.5 && abs(eta)<=2.5) * (pt > 50 && pt <= 60) * (0.68) + \
		           (abs(eta) > 1.5 && abs(eta)<=2.5) * (pt > 60 && pt <= 70) * (0.78) + \
		           (abs(eta) > 1.5 && abs(eta)<=2.5) * (pt > 70 && pt <= 80) * (0.84) + \
		           (abs(eta) > 1.5 && abs(eta)<=2.5) * (pt > 80 && pt <= 90) * (0.89) + \
		           (abs(eta) > 1.5 && abs(eta)<=2.5) * (pt > 90 && pt <= 100) * (0.92) + \
		           (abs(eta) > 1.5 && abs(eta)<=2.5) * (pt > 100 && pt <= 110) * (0.93) + \
                           (abs(eta) > 2.5)                                   * (0.00)
      
  }
}

####################
# Photon isolation #
####################

module Isolation PhotonIsolation {
  # particle for which calculate the isolation
  set CandidateInputArray        PhotonEfficiency/photons 
  # neutral and charged particles for the whole event (no CHS applied)
  set NeutralIsolationInputArray Calorimeter/eflowTowers
  set ChargedIsolationInputArray TrackMergerWithMuon/tracks
  # select a rho for the isolation
  set RhoInputArray RhoKt4/rho
  # output array
  set OutputArray photons
  # isolation cone
  set DeltaRMax 0.3
  # minimum pT  
  set PTMin     0.5
  # iso ratio to cut
  set PTRatioMax 9999.
}


#######################
# Electron efficiency #
#######################

module Efficiency ElectronEfficiency {
  set InputArray TrackPileUpSubtractor/electrons
  set OutputArray electrons
  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
  # efficiency formula for electrons
  set EfficiencyFormula {(pt <= 4.0)  * (0.00) + \
                         (abs(eta) <= 1.45 ) * (pt >  4.0 && pt <= 6.0)   * (0.5*0.50) + \
                         (abs(eta) <= 1.45 ) * (pt >  6.0 && pt <= 8.0)   * (0.5*0.70) + \
                         (abs(eta) <= 1.45 ) * (pt >  8.0 && pt <= 10.0)  * (0.7*0.85) + \
                         (abs(eta) <= 1.45 ) * (pt > 10.0 && pt <= 30.0)  * (0.90*0.94) + \                                                      
                         (abs(eta) <= 1.45 ) * (pt > 30.0 && pt <= 50.0)  * (0.95*0.97) + \                          
                         (abs(eta) <= 1.45 ) * (pt > 50.0 && pt <= 70.0)  * (0.95*0.98) + \          
                         (abs(eta) <= 1.45 ) * (pt > 70.0 )  * (1.0) + \                                                                                                
                         (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt >  4.0 && pt <= 10.0)   * (0.5*0.35) + \
                         (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt > 10.0 && pt <= 30.0)   * (0.5*0.40) + \   
                         (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt > 30.0 && pt <= 70.0)   * (0.8*0.45) + \                                 
                         (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt > 70.0 )  * (0.8*0.45) + \    
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt >  4.0 && pt <= 10.0)  * (0.7*0.75) + \
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 10.0 && pt <= 30.0)  * (0.80*0.85) + \                                                      
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 30.0 && pt <= 50.0)  * (0.85*0.95) + \                          
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 50.0 && pt <= 70.0)  * (0.85*0.95) + \          
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 70.0 )  * (0.85*1.0) + \   
                         (abs(eta) >= 2.0 && abs(eta) <= 2.5 ) * (pt >  4.0 && pt <= 10.0)  * (0.7*0.65) + \
                         (abs(eta) >= 2.0 && abs(eta) <= 2.5 ) * (pt > 10.0 && pt <= 30.0)  * (0.7*0.75) + \                                                      
                         (abs(eta) >= 2.0 && abs(eta) <= 2.5 ) * (pt > 30.0 && pt <= 50.0)  * (0.8*0.85) + \                          
                         (abs(eta) >= 2.0 && abs(eta) <= 2.5 ) * (pt > 50.0 && pt <= 70.0)  * (0.8*0.85) + \          
                         (abs(eta) >= 2.0 && abs(eta) <= 2.5 ) * (pt > 70.0 )  * (0.8*0.85) + \                                                                         
                         (abs(eta) > 2.5)                              * (0.00)
  }
}

######################
# Electron isolation #
######################

module Isolation ElectronIsolation {
  set CandidateInputArray        ElectronEfficiency/electrons
  set NeutralIsolationInputArray Calorimeter/eflowTowers
  set ChargedIsolationInputArray TrackMergerWithMuon/tracks
  set RhoInputArray RhoKt4/rho
  set OutputArray electrons
  set DeltaRMax 0.3
  set PTMin 0.5
  set PTRatioMax 9999.
}

###################
# Muon efficiency #
###################

module Efficiency MuonEfficiency {
  set InputArray TrackPileUpSubtractor/muons
  set OutputArray muons
  # set EfficiencyFormula {efficiency as a function of eta and pt}
  # efficiency formula for muons
  set EfficiencyFormula {
                                  (pt <= 10.0)  * (0.00) + \  
                         (abs(eta)<=0.1)*(pt>10)*(0.89865) + \
                         (abs(eta)>0.1 && abs(eta)<=0.2)*(pt>10)*(0.894596) + \
			 (abs(eta)>0.2 && abs(eta)<=0.3)*(pt>10)*(0.764087) + \
                         (abs(eta)>0.3 && abs(eta)<=0.4)*(pt>10)*(0.881295) + \
                         (abs(eta)>0.4 && abs(eta)<=0.5)*(pt>10)*(0.913192) + \
                         (abs(eta)>0.5 && abs(eta)<=0.6)*(pt>10)*(0.897579) + \
                         (abs(eta)>0.6 && abs(eta)<=0.7)*(pt>10)*(0.894978) + \
                         (abs(eta)>0.7 && abs(eta)<=0.8)*(pt>10)*(0.878466) + \
                         (abs(eta)>0.8 && abs(eta)<=0.9)*(pt>10)*(0.831849) + \
                         (abs(eta)>0.9 && abs(eta)<=1.0)*(pt>10)*(0.806424) + \
                         (abs(eta)>1.0 && abs(eta)<=1.1)*(pt>10)*(0.756892) + \      
                         (abs(eta)>1.1 && abs(eta)<=1.2)*(pt>10)*(0.728583) + \
                         (abs(eta)>1.2 && abs(eta)<=1.3)*(pt>10)*(0.773855) + \
                         (abs(eta)>1.3 && abs(eta)<=1.4)*(pt>10)*(0.776296) + \
                         (abs(eta)>1.4 && abs(eta)<=1.5)*(pt>10)*(0.769977) + \
                         (abs(eta)>1.5 && abs(eta)<=1.6)*(pt>10)*(0.838174) + \
                         (abs(eta)>1.6 && abs(eta)<=1.7)*(pt>10)*(0.854358) + \
                         (abs(eta)>1.7 && abs(eta)<=1.8)*(pt>10)*(0.8565) + \
                         (abs(eta)>1.8 && abs(eta)<=1.9)*(pt>10)*(0.857182) + \
                         (abs(eta)>1.9 && abs(eta)<=2.0)*(pt>10)*(0.85591) + \
                         (abs(eta)>2.0 && abs(eta)<=2.1)*(pt>10)*(0.844826) + \
                         (abs(eta)>2.1 && abs(eta)<=2.2)*(pt>10)*(0.81742) + \
                         (abs(eta)>2.2 && abs(eta)<=2.3)*(pt>10)*(0.825831) + \
                         (abs(eta)>2.3 && abs(eta)<=2.4)*(pt>10)*(0.774208) + \
                         (abs(eta) > 2.40)  * (0.00)
  }
}

##################
# Muon isolation #
##################

module Isolation MuonIsolation {
  set CandidateInputArray MuonEfficiency/muons
  set NeutralIsolationInputArray Calorimeter/eflowTowers
  set ChargedIsolationInputArray TrackMergerWithMuon/tracks 
  set RhoInputArray RhoKt4/rho
  set OutputArray muons
  set DeltaRMax 0.3
  set PTMin 0.5
  set PTRatioMax 9999.
}


#############
# b-tagging #
#############
module BTagging BTagging {

  set JetInputArray JetPileUpSubtractor/jets

  add EfficiencyFormulaLoose {0} {
      (pt <= 20.0) * (0.000) + \
			     (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.0965) + \
			     (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.105) + \
			     (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.0762) + \
			     (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.0851) + \
			     (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.0739) + \
			     (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.0808) + \
			     (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.0863) + \
			     (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.0824) + \
			     (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.0891) + \
			     (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.0977) + \
			     (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.1028) + \
			     (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.105) + \
			     (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.1116) + \
			     (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.1225) + \
			     (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.1384) + \
			     (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.1535) + \
			     (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.1693) + \
			     (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.1869) + \
			     (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.2111) + \
			     (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.2063) + \
			     (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.2132) + \
			     (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.216) + \
			     (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.2273) + \
			     (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.2686) + \
			     (abs(eta) <= 1.8) * (pt > 2000.0) * (0.3134) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.054932) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.078226) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.059492) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.07064) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.071246) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.081144) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.088663) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.080107) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.087845) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.099813) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.103151) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.101119) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.109951) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.120709) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.1346) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.1524) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.165067) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.108622) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.124293) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.0823) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.086487) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.09222) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.000) + \
			     (abs(eta) > 2.4) * (0.000)
  }

  add EfficiencyFormulaLoose {4} {             (pt <= 20.0) * (0.000) + \
						   (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.307) + \
						   (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.355) + \
						   (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.315) + \
                                                (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.332) + \
                                                (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.313) + \
                                                (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.322) + \
                                                (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.331) + \
                                                (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.312) + \
                                                (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.319) + \
                                                (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.329) + \
                                                (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.317) + \
                                                (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.301) + \
                                                (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.306) + \
                                                (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.309) + \
                                                (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.309) + \
                                                (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.309) + \
                                                (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.313) + \
                                                (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.308) + \
                                                (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.321) + \
                                                (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.287) + \
                                                (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.295) + \
                                                (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.278) + \
                                                (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.293) + \                                    
                                                (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.351) + \
                                                (abs(eta) <= 1.8) * (pt > 2000.0) * (0.388) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.15416) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.20465) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.17009) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.18172) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.19284) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.19356) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.20196) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.18933) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.19708) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.20503) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.20163) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.18223) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.18792) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.19688) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.21584) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.22609) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.24573) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.15426) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.17006) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.14041) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.10447) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.15677) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.) + \
                                                (abs(eta) > 2.4) * (0.000)
  }

   # efficiency formula for b-jets
  add EfficiencyFormulaLoose {5} {                   (pt <= 20.0) * (0.000) + \
                                                (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.634) + \
                                                (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.723) + \
                                                (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.721) + \
                                                (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.747) + \
                                                (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.745) + \
                                                (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.755) + \
                                                (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.762) + \
                                                (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.762) + \
                                                (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.753) + \
                                                (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.75) + \
                                                (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.738) + \
                                                (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.723) + \
                                                (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.714) + \
                                                (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.691) + \
                                                (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.669) + \
                                                (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.646) + \
                                                (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.625) + \
                                                (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.614) + \
                                                (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.585) + \
                                                (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.519) + \
                                                (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.494) + \
                                                (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.453) + \
                                                (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.438) + \
                                                (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.486) + \
                                                (abs(eta) <= 1.8) * (pt > 2000.0) * (0.541) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.4585) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.5768) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.577) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.6064) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.6202) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.6085) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.6178) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.5966) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.587) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.5785) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.5605) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.5103) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.5111) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.4889) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.4697) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.4361) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.4178) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.3698) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.3255) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.2703) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.2767) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.2941) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) + \
                                                (abs(eta) > 2.4) * (0.000)
  }
 
  add EfficiencyFormulaMedium {0} {                                (pt <= 20.0) * (0.000) + \
                                                (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.00469) + \
                                                (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.00691) + \
                                                (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.00519) + \
                                                (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.00633) + \
                                                (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.00574) + \
                                                (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.00656) + \
                                                (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.0073) + \
                                                (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.00606) + \
                                                (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.00722) + \
                                                (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.00837) + \
                                                (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.01031) + \
                                                (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.01224) + \
                                                (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.01351) + \
                                                (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.01542) + \
                                                (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.01796) + \
                                                (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.02099) + \
                                                (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.0246) + \
                                                (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.01638) + \
                                                (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.01989) + \
                                                (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.01629) + \
                                                (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.01773) + \
                                                (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.01977) + \
                                                (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.02372) + \
                                                (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.0323) + \
                                                (abs(eta) <= 1.8) * (pt > 2000.0) * (0.04635) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.04635) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.004389) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.004706) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.00583) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.004895) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.006023) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.006487) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.005549) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.006939) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.008245) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.009879) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.011744) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.012714) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.014575) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.01848) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.022346) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.024952) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.007563) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.010131) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.004863) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.006965) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.007071) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.000) + \
                                                (abs(eta) > 2.4) * (0.000)
  }

 # efficiency formula for c-jets (misidentification rate)
  add EfficiencyFormulaMedium {4} {             (pt <= 20.0) * (0.000) + \
                                                (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.0534) + \
                                                (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.0677) + \
                                                (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.0559) + \
                                                (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.0616) + \
                                                (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.0603) + \
                                                (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.0641) + \
                                                (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.0647) + \
                                                (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.064) + \
                                                (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.066) + \
                                                (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.0666) + \
                                                (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.0679) + \
                                                (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.0701) + \
                                                (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.0664) + \
                                                (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.0688) + \
                                                (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.0671) + \
                                                (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.0654) + \
                                                (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.0651) + \
                                                (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.0452) + \
                                                (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.0484) + \
                                                (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.0346) + \
                                                (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.0357) + \
                                                (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.035) + \
                                                (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.0425) + \                                    
                                                (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.0635) + \
                                                (abs(eta) <= 1.8) * (pt > 2000.0) * (0.0951) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.0124) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.01787) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.01962) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.01831) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.01842) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.0224) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.0198) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.02005) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.02146) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.02519) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.02979) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.03011) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.03065) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.0338) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.03664) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.04036) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.04268) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.0142) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.00971) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.00759) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.00746) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.00423) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) + \
                                                (abs(eta) > 2.4) * (0.000)
  }

 # efficiency formula for b-jets
  add EfficiencyFormulaMedium {5} {                       (pt <= 20.0) * (0.000) + \
                                                (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.3392) + \
                                                (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.4447) + \
                                                (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.4628) + \
                                                (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.489) + \
                                                (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.5029) + \
                                                (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.5074) + \
                                                (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.5154) + \
                                                (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.5077) + \
                                                (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.5028) + \
                                                (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.4922) + \
                                                (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.4739) + \
                                                (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.4623) + \
                                                (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.4415) + \
                                                (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.4134) + \
                                                (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.3822) + \
                                                (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.351) + \
                                                (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.3212) + \
                                                (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.2507) + \
                                                (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.2098) + \
                                                (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.154) + \
                                                (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.1472) + \
                                                (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.136) + \
                                                (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.142) + \
                                                (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.1915) + \
                                                (abs(eta) <= 1.8) * (pt > 2000.0) * (0.2249) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.1792) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.2611) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.2846) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.2907) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.2949) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.2875) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.2812) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.2927) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.2668) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.2832) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.2488) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.2297) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.2106) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.1991) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.1764) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.1779) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.1569) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.0812) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.0634) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.0444) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.0625) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.0661) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) + \
                                                (abs(eta) > 2.4) * (0.000)
  }

  add EfficiencyFormulaTight {0} {
                             (pt <= 20.0) * (0.000) + \
                                                (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.0002) + \
                                                (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.00027) + \
                                                (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.000278) + \
                                                (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.000343) + \
                                                (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.000343) + \
                                                (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.00045) + \
                                                (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.000469) + \
                                                (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.000449) + \
                                                (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.000556) + \
                                                (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.000731) + \
                                                (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.000929) + \
                                                (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.001311) + \
                                                (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.00152) + \
                                                (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.001657) + \
                                                (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.002124) + \
                                                (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.00254) + \
                                                (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.00292) + \
                                                (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.00116) + \
                                                (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.001368) + \
                                                (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.001201) + \
                                                (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.001249) + \
                                                (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.001548) + \
                                                (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.001898) + \
                                                (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.003125) + \
                                                (abs(eta) <= 1.8) * (pt > 2000.0) * (0.004864) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.000476) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.000538) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.000468) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.000687) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.000624) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.00072) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.0008) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.000572) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.000843) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.00101) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.000999) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.000763) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.001088) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.001204) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.001871) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.00216) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.003148) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.003421) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.004692) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.005582) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.005732) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.007186) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.000) + \
                                                (abs(eta) > 2.4) * (0.000)
  }

 # efficiency formula for c-jets (misidentification rate)
  add EfficiencyFormulaTight {4} { (pt <= 20.0) * (0.000) + \
                                                (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.00329) + \
                                                (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.00403) + \
                                                (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.00373) + \
                                                (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.00437) + \
                                                (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.00525) + \
                                                (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.0049) + \
                                                (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.00506) + \
                                                (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.00559) + \
                                                (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.00605) + \
                                                (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.0069) + \
                                                (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.00725) + \
                                                (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.00805) + \
                                                (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.00741) + \
                                                (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.00763) + \
                                                (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.00872) + \
                                                (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.00731) + \
                                                (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.00773) + \
                                                (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.00383) + \
                                                (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.00377) + \
                                                (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.00239) + \
                                                (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.00264) + \
                                                (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.00266) + \
                                                (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.00362)+ \                                    
                                                (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.00498) + \
                                                (abs(eta) <= 1.8) * (pt > 2000.0) * (0.01455) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.00387) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.00553) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.00654) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.00657) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.00629) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.00595) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.00533) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.00361) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.00416) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.00658) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.0044) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.0036) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.00154) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.0028) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.00296) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.00352) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.00731) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.0044) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.01068) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.01138) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.00746) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.00847) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) + \
                                                (abs(eta) > 2.4) * (0.000)
  }

  # efficiency formula for b-jets
  add EfficiencyFormulaTight {5} { (pt <= 20.0) * (0.000) + \
                                                (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.1371) + \
                                                (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.1973) + \
                                                (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.2189) + \
                                                (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.231) + \
                                                (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.2494) + \
                                                (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.2514) + \
                                                (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.2529) + \
                                                (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.2482) + \
                                                (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.2464) + \
                                                (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.2328) + \
                                                (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.212) + \
                                                (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.1854) + \
                                                (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.1706) + \
                                                (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.1559) + \
                                                (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.1361) + \
                                                (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.1203) + \
                                                (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.1065) + \
                                                (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.0534) + \
                                                (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.0396) + \
                                                (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.0277) + \
                                                (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.0303) + \
                                                (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.0288) + \
                                                (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.0335) + \
                                                (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.0445) + \
                                                (abs(eta) <= 1.8) * (pt > 2000.0) * (0.0645) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.0804) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.1354) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.1715) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.182) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.1832) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.1818) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.1648) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.1621) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.1414) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.1446) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.1069) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.079) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.0736) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.0626) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.0484) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.0459) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.0384) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.0319) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.0401) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.037) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.0446) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.0661) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) + \
                                                (abs(eta) > 2.4) * (0.000)
  }

}







module BTagging PuppiBTagging {

  set JetInputArray PuppiJetPileUpSubtractor/jets

  add EfficiencyFormulaLoose {0} {
      (pt <= 20.0) * (0.000) + \
			     (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.0965) + \
			     (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.105) + \
			     (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.0762) + \
			     (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.0851) + \
			     (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.0739) + \
			     (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.0808) + \
			     (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.0863) + \
			     (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.0824) + \
			     (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.0891) + \
			     (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.0977) + \
			     (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.1028) + \
			     (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.105) + \
			     (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.1116) + \
			     (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.1225) + \
			     (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.1384) + \
			     (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.1535) + \
			     (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.1693) + \
			     (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.1869) + \
			     (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.2111) + \
			     (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.2063) + \
			     (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.2132) + \
			     (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.216) + \
			     (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.2273) + \
			     (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.2686) + \
			     (abs(eta) <= 1.8) * (pt > 2000.0) * (0.3134) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.054932) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.078226) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.059492) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.07064) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.071246) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.081144) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.088663) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.080107) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.087845) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.099813) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.103151) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.101119) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.109951) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.120709) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.1346) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.1524) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.165067) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.108622) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.124293) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.0823) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.086487) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.09222) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
			     (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.000) + \
			     (abs(eta) > 2.4) * (0.000)
  }

  add EfficiencyFormulaLoose {4} {             (pt <= 20.0) * (0.000) + \
						   (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.307) + \
						   (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.355) + \
						   (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.315) + \
                                                (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.332) + \
                                                (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.313) + \
                                                (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.322) + \
                                                (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.331) + \
                                                (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.312) + \
                                                (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.319) + \
                                                (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.329) + \
                                                (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.317) + \
                                                (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.301) + \
                                                (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.306) + \
                                                (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.309) + \
                                                (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.309) + \
                                                (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.309) + \
                                                (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.313) + \
                                                (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.308) + \
                                                (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.321) + \
                                                (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.287) + \
                                                (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.295) + \
                                                (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.278) + \
                                                (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.293) + \                                    
                                                (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.351) + \
                                                (abs(eta) <= 1.8) * (pt > 2000.0) * (0.388) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.15416) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.20465) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.17009) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.18172) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.19284) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.19356) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.20196) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.18933) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.19708) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.20503) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.20163) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.18223) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.18792) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.19688) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.21584) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.22609) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.24573) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.15426) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.17006) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.14041) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.10447) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.15677) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.) + \
                                                (abs(eta) > 2.4) * (0.000)
  }

   # efficiency formula for b-jets
  add EfficiencyFormulaLoose {5} {                   (pt <= 20.0) * (0.000) + \
                                                (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.634) + \
                                                (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.723) + \
                                                (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.721) + \
                                                (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.747) + \
                                                (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.745) + \
                                                (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.755) + \
                                                (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.762) + \
                                                (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.762) + \
                                                (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.753) + \
                                                (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.75) + \
                                                (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.738) + \
                                                (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.723) + \
                                                (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.714) + \
                                                (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.691) + \
                                                (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.669) + \
                                                (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.646) + \
                                                (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.625) + \
                                                (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.614) + \
                                                (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.585) + \
                                                (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.519) + \
                                                (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.494) + \
                                                (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.453) + \
                                                (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.438) + \
                                                (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.486) + \
                                                (abs(eta) <= 1.8) * (pt > 2000.0) * (0.541) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.4585) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.5768) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.577) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.6064) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.6202) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.6085) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.6178) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.5966) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.587) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.5785) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.5605) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.5103) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.5111) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.4889) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.4697) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.4361) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.4178) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.3698) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.3255) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.2703) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.2767) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.2941) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) + \
                                                (abs(eta) > 2.4) * (0.000)
  }
 
  add EfficiencyFormulaMedium {0} {                                (pt <= 20.0) * (0.000) + \
                                                (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.00469) + \
                                                (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.00691) + \
                                                (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.00519) + \
                                                (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.00633) + \
                                                (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.00574) + \
                                                (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.00656) + \
                                                (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.0073) + \
                                                (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.00606) + \
                                                (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.00722) + \
                                                (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.00837) + \
                                                (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.01031) + \
                                                (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.01224) + \
                                                (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.01351) + \
                                                (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.01542) + \
                                                (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.01796) + \
                                                (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.02099) + \
                                                (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.0246) + \
                                                (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.01638) + \
                                                (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.01989) + \
                                                (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.01629) + \
                                                (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.01773) + \
                                                (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.01977) + \
                                                (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.02372) + \
                                                (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.0323) + \
                                                (abs(eta) <= 1.8) * (pt > 2000.0) * (0.04635) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.04635) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.004389) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.004706) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.00583) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.004895) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.006023) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.006487) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.005549) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.006939) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.008245) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.009879) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.011744) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.012714) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.014575) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.01848) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.022346) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.024952) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.007563) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.010131) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.004863) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.006965) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.007071) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.000) + \
                                                (abs(eta) > 2.4) * (0.000)
  }

 # efficiency formula for c-jets (misidentification rate)
  add EfficiencyFormulaMedium {4} {             (pt <= 20.0) * (0.000) + \
                                                (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.0534) + \
                                                (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.0677) + \
                                                (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.0559) + \
                                                (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.0616) + \
                                                (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.0603) + \
                                                (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.0641) + \
                                                (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.0647) + \
                                                (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.064) + \
                                                (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.066) + \
                                                (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.0666) + \
                                                (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.0679) + \
                                                (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.0701) + \
                                                (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.0664) + \
                                                (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.0688) + \
                                                (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.0671) + \
                                                (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.0654) + \
                                                (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.0651) + \
                                                (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.0452) + \
                                                (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.0484) + \
                                                (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.0346) + \
                                                (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.0357) + \
                                                (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.035) + \
                                                (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.0425) + \                                    
                                                (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.0635) + \
                                                (abs(eta) <= 1.8) * (pt > 2000.0) * (0.0951) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.0124) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.01787) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.01962) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.01831) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.01842) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.0224) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.0198) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.02005) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.02146) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.02519) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.02979) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.03011) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.03065) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.0338) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.03664) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.04036) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.04268) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.0142) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.00971) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.00759) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.00746) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.00423) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) + \
                                                (abs(eta) > 2.4) * (0.000)
  }

 # efficiency formula for b-jets
  add EfficiencyFormulaMedium {5} {                       (pt <= 20.0) * (0.000) + \
                                                (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.3392) + \
                                                (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.4447) + \
                                                (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.4628) + \
                                                (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.489) + \
                                                (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.5029) + \
                                                (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.5074) + \
                                                (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.5154) + \
                                                (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.5077) + \
                                                (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.5028) + \
                                                (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.4922) + \
                                                (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.4739) + \
                                                (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.4623) + \
                                                (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.4415) + \
                                                (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.4134) + \
                                                (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.3822) + \
                                                (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.351) + \
                                                (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.3212) + \
                                                (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.2507) + \
                                                (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.2098) + \
                                                (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.154) + \
                                                (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.1472) + \
                                                (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.136) + \
                                                (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.142) + \
                                                (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.1915) + \
                                                (abs(eta) <= 1.8) * (pt > 2000.0) * (0.2249) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.1792) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.2611) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.2846) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.2907) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.2949) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.2875) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.2812) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.2927) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.2668) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.2832) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.2488) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.2297) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.2106) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.1991) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.1764) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.1779) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.1569) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.0812) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.0634) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.0444) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.0625) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.0661) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) + \
                                                (abs(eta) > 2.4) * (0.000)
  }

  add EfficiencyFormulaTight {0} {
                             (pt <= 20.0) * (0.000) + \
                                                (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.0002) + \
                                                (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.00027) + \
                                                (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.000278) + \
                                                (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.000343) + \
                                                (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.000343) + \
                                                (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.00045) + \
                                                (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.000469) + \
                                                (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.000449) + \
                                                (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.000556) + \
                                                (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.000731) + \
                                                (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.000929) + \
                                                (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.001311) + \
                                                (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.00152) + \
                                                (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.001657) + \
                                                (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.002124) + \
                                                (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.00254) + \
                                                (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.00292) + \
                                                (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.00116) + \
                                                (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.001368) + \
                                                (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.001201) + \
                                                (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.001249) + \
                                                (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.001548) + \
                                                (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.001898) + \
                                                (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.003125) + \
                                                (abs(eta) <= 1.8) * (pt > 2000.0) * (0.004864) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.000476) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.000538) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.000468) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.000687) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.000624) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.00072) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.0008) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.000572) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.000843) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.00101) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.000999) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.000763) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.001088) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.001204) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.001871) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.00216) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.003148) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.003421) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.004692) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.005582) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.005732) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.007186) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.000) + \
                                                (abs(eta) > 2.4) * (0.000)
  }

 # efficiency formula for c-jets (misidentification rate)
  add EfficiencyFormulaTight {4} { (pt <= 20.0) * (0.000) + \
                                                (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.00329) + \
                                                (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.00403) + \
                                                (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.00373) + \
                                                (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.00437) + \
                                                (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.00525) + \
                                                (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.0049) + \
                                                (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.00506) + \
                                                (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.00559) + \
                                                (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.00605) + \
                                                (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.0069) + \
                                                (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.00725) + \
                                                (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.00805) + \
                                                (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.00741) + \
                                                (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.00763) + \
                                                (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.00872) + \
                                                (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.00731) + \
                                                (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.00773) + \
                                                (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.00383) + \
                                                (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.00377) + \
                                                (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.00239) + \
                                                (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.00264) + \
                                                (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.00266) + \
                                                (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.00362)+ \                                    
                                                (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.00498) + \
                                                (abs(eta) <= 1.8) * (pt > 2000.0) * (0.01455) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.00387) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.00553) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.00654) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.00657) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.00629) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.00595) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.00533) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.00361) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.00416) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.00658) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.0044) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.0036) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.00154) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.0028) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.00296) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.00352) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.00731) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.0044) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.01068) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.01138) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.00746) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.00847) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) + \
                                                (abs(eta) > 2.4) * (0.000)
  }

  # efficiency formula for b-jets
  add EfficiencyFormulaTight {5} { (pt <= 20.0) * (0.000) + \
                                                (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.1371) + \
                                                (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.1973) + \
                                                (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.2189) + \
                                                (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.231) + \
                                                (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.2494) + \
                                                (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.2514) + \
                                                (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.2529) + \
                                                (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.2482) + \
                                                (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.2464) + \
                                                (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.2328) + \
                                                (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.212) + \
                                                (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.1854) + \
                                                (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.1706) + \
                                                (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.1559) + \
                                                (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.1361) + \
                                                (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.1203) + \
                                                (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.1065) + \
                                                (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.0534) + \
                                                (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.0396) + \
                                                (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.0277) + \
                                                (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.0303) + \
                                                (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.0288) + \
                                                (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.0335) + \
                                                (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.0445) + \
                                                (abs(eta) <= 1.8) * (pt > 2000.0) * (0.0645) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.0804) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.1354) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.1715) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.182) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.1832) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.1818) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.1648) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.1621) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.1414) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.1446) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.1069) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.079) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.0736) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.0626) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.0484) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.0459) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.0384) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.0319) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.0401) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.037) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.0446) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.0661) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
                                                (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) + \
                                                (abs(eta) > 2.4) * (0.000)
  }

}



###################
## PileUp Jet ID ##
###################

module PileUpJetID PileUpJetID {
  set JetInputArray     JetPileUpSubtractor/jets
  set TrackInputArray   TrackMergerWithMuon/tracks
  set NeutralInputArray Calorimeter/eflowTowers
  set PVInputArray      ModifyBeamSpot/PV

  set OutputArray   jets

  set ParameterR      0.4  
  set JetPTMin        20.0
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

  set ParameterR      0.4  
  set JetPTMin        20.0
  set UseConstituents 0
 
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


module Merger PileUpJetIDMissingET {
  add InputArray TrackPileUpSubtractor/eflowTracks
  add InputArray MuonMomentumSmearing/muons
  add InputArray PileUpJetID/eflowTowers
  set MomentumOutputArray momentum
}

########################
## Constituent filter ##
########################

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

##############
## Scalr HT ##
##############

module Merger ScalarHT {
 # add InputArray InputArray
 add InputArray  EFlowMerger/eflow
 set EnergyOutputArray energy
}

module Merger GenScalarHT {
 # add InputArray InputArray
 add InputArray NeutrinoFilter/stableParticles  
 set EnergyOutputArray energy
}

module Merger PuppiScalarHT {
 # add InputArray InputArray
 add InputArray RunPUPPI/PuppiParticles
 set EnergyOutputArray energy
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
  #add Branch StatusPid/filteredParticles GenParticles GenParticle
 
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
  add Branch PuppiRhoGrid/rho PuppiRhoGrid Rho
  #add Branch PuppiJetFinder/jets     RawPuppiJet Jet
  #add Branch PuppiJetPileUpSubtractor/jets PuppiJet Jet
  #add Branch PuppiJetPileUpSubtractorGrid/jets PuppiJetGrid Jet
  #add Branch PuppiJetPileUpSubtractor4VArea/jets PuppiJet4VArea Jet
  add Branch PuppiPileUpJetID/jets PuppiJetPUID Jet

  ## MET
  add Branch GenMissingET/momentum GenMissingET MissingET
  add Branch MissingET/momentum MissingET MissingET
  add Branch PuppiMissingET/momentum PuppiMissingET MissingET

  ## HT
  add Branch ScalarHT/energy HT ScalarHT
  add Branch GenScalarHT/energy GenHT ScalarHT
  add Branch PuppiScalarHT/energy PuppiHT ScalarHT

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


