#ifndef DelphesClasses_h
#define DelphesClasses_h

/**
 *  Definition of classes to be stored in the root tree.
 *  Function CompareXYZ sorts objects by the variable XYZ that MUST be
 *  present in the data members of the root tree class of the branch.
 *  $Date: 2008-06-04 13:57:24 $
 *  $Revision: 1.1 $
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

// Dependencies (#includes)

#include "TRef.h"
#include "TObject.h"
#include "TRefArray.h"
#include "TLorentzVector.h"

#include <vector>

#include "classes/SortableObject.h"

namespace fastjet {
  class PseudoJet;
}


class DelphesFactory;

//---------------------------------------------------------------------------

class Event: public TObject{

 public:

  Long64_t Number; // event number
  Float_t ReadTime;
  Float_t ProcTime;  
  ClassDef(Event, 1)
};

//---------------------------------------------------------------------------

class LHCOEvent: public Event{
 public:
  Int_t Trigger; // trigger word
  ClassDef(LHCOEvent, 1)
};

//---------------------------------------------------------------------------

class LHEFrwgt : public TObject {
  
 public:
  
  std::vector<int> opNum;
  std::vector<double> opVal;
  ClassDef(LHEFrwgt, 1)

};

//---------------------------------------------------------------------------

class LHEFEvent: public Event {
  
  public:

  Int_t   ProcessID; // subprocess code for the event | hepup.IDPRUP
  Float_t Weight; // weight for the event | hepup.XWGTUP
  Float_t ScalePDF; // scale in GeV used in the calculation of the PDFs in the event | hepup.SCALUP
  Float_t AlphaQED; // value of the QED coupling used in the event | hepup.AQEDUP
  Float_t AlphaQCD; // value of the QCD coupling used in the event | hepup.AQCDUP
  std::vector<Double_t> lheWeights;
  ClassDef(LHEFEvent, 3)
};

//---------------------------------------------------------------------------

class HepMCEvent: public Event{

  public:

  Int_t ProcessID; // unique signal process id | signal_process_id()
  Int_t MPI; // number of multi parton interactions | mpi () 

  Float_t Weight; // weight for the event

  Float_t Scale; // energy scale, see hep-ph/0109068 | event_scale()
  Float_t AlphaQED; // QED coupling, see hep-ph/0109068 | alphaQED()
  Float_t AlphaQCD; // QCD coupling, see hep-ph/0109068 | alphaQCD()

  Int_t ID1; // flavour code of first parton | pdf_info()->id1()
  Int_t ID2; // flavour code of second parton | pdf_info()->id2()     

  Float_t X1; // fraction of beam momentum carried by first parton ("beam side") | pdf_info()->x1()
  Float_t X2; // fraction of beam momentum carried by second parton ("target side") | pdf_info()->x2()

  Float_t ScalePDF; // Q-scale used in evaluation of PDF's (in GeV) | pdf_info()->scalePDF()

  Float_t PDF1; // PDF (id1, x1, Q) | pdf_info()->pdf1()
  Float_t PDF2; // PDF (id2, x2, Q) | pdf_info()->pdf2()

  ClassDef(HepMCEvent, 2)
};


//---------------------------------------------------------------------------
class GenParticle: public SortableObject{
  public:
  Int_t PID; // particle HEP ID number | hepevt.idhep[number]

  Int_t Status; // particle status | hepevt.isthep[number]
  Int_t IsPU; // 0 or 1 for particles from pile-up interactions
                  // Reflects decision ncoprorating z vertex resolution  
  
  Int_t M1; // particle 1st mother | hepevt.jmohep[number][0] - 1
  Int_t M2; // particle 2nd mother | hepevt.jmohep[number][1] - 1

  Int_t D1; // particle 1st daughter | hepevt.jdahep[number][0] - 1
  Int_t D2; // particle last daughter | hepevt.jdahep[number][1] - 1

  Int_t Charge; // particle charge

  Float_t Mass; // particle mass

  Float_t E; // particle energy | hepevt.phep[number][3]
  Float_t Px; // particle momentum vector (x component) | hepevt.phep[number][0]
  Float_t Py; // particle momentum vector (y component) | hepevt.phep[number][1]
  Float_t Pz; // particle momentum vector (z component) | hepevt.phep[number][2]

  Float_t PT; // particle transverse momentum
  Float_t Eta; // particle pseudorapidity
  Float_t Phi; // particle azimuthal angle

  Float_t Rapidity; // particle rapidity

  Float_t T; // particle vertex position (t component) | hepevt.vhep[number][3]
  Float_t X; // particle vertex position (x component) | hepevt.vhep[number][0]
  Float_t Y; // particle vertex position (y component) | hepevt.vhep[number][1]
  Float_t Z; // particle vertex position (z component) | hepevt.vhep[number][2]

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }
  
  TLorentzVector P4();

  ClassDef(GenParticle, 1)
};


//---------------------------------------------------------------------------
class LHEParticle: public SortableObject{
  public:

  Int_t PID;    // particle HEP ID number | hepevt.idhep[number]
  Int_t Status; // particle status | hepevt.isthep[number]
  Int_t IsPU;   // 0 or 1 for particles from pile-up interactions // Reflects decision ncoprorating z vertex resolution  
  
  Int_t M1; // particle 1st mother | hepevt.jmohep[number][0] - 1
  Int_t M2; // particle 2nd mother | hepevt.jmohep[number][1] - 1

  Int_t D1; // particle 1st daughter | hepevt.jdahep[number][0] - 1
  Int_t D2; // particle last daughter | hepevt.jdahep[number][1] - 1

  Int_t Spin ;

  Int_t Charge; // particle charge
  Float_t Mass; // particle mass

  Float_t E; // particle energy | hepevt.phep[number][3]
  Float_t Px; // particle momentum vector (x component) | hepevt.phep[number][0]
  Float_t Py; // particle momentum vector (y component) | hepevt.phep[number][1]
  Float_t Pz; // particle momentum vector (z component) | hepevt.phep[number][2]

  Float_t PT; // particle transverse momentum
  Float_t Eta; // particle pseudorapidity
  Float_t Phi; // particle azimuthal angle

  Float_t Rapidity; // particle rapidity

  Float_t T; // particle vertex position (t component) | hepevt.vhep[number][3]
  Float_t X; // particle vertex position (x component) | hepevt.vhep[number][0]
  Float_t Y; // particle vertex position (y component) | hepevt.vhep[number][1]
  Float_t Z; // particle vertex position (z component) | hepevt.vhep[number][2]

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }
  
  TLorentzVector P4();

  ClassDef(LHEParticle, 1)
};

//---------------------------------------------------------------------------

class MissingET: public TObject{
 public:
  Float_t MET; // mising transverse energy
  Float_t Phi; // mising energy azimuthal angle

  ClassDef(MissingET, 1)
};

//---------------------------------------------------------------------------

class ScalarHT: public TObject{

 public:
  Float_t HT; // scalar sum of transverse momenta

  ClassDef(ScalarHT, 1)
};

//---------------------------------------------------------------------------
class Rho: public TObject{
 public:
  Float_t Rho; // rho energy density
  Float_t Edges[2]; // pseudorapidity range edges
  
  ClassDef(Rho, 1)
};

//---------------------------------------------------------------------------

class Weight: public TObject {
public:
  Float_t Weight; // weight for the event

  ClassDef(Weight, 1)
};

//---------------------------------------------------------------------------

class Photon: public SortableObject {

 public:

  Float_t PT; // photon transverse momentum
  Float_t Eta; // photon pseudorapidity
  Float_t Phi; // photon azimuthal angle
  Float_t E; // photon energy

  Int_t Status; // particle status | hepevt.isthep[number]                                                                                                                                 
  Int_t IsRecoPU; // Assigned by charged hadron subtractor where applicable                                                                                                                
  Int_t IsPU;
  Int_t IsFakeObject;
  
  Float_t EhadOverEem; // ratio of the hadronic versus electromagnetic energy deposited in the calorimeter

  Float_t IsolationVarDBeta;
  Float_t IsolationVarRhoCorr;
  Float_t TrackIsolationVar;
  Float_t chargedHadronEnergy;
  Float_t neutralEnergy;
  Float_t chargedPUEnergy;
  Float_t allParticleEnergy;

  Float_t TOuter; // time at calorimeter face

  TRefArray Particles; // references to generated particles

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4();

  ClassDef(Photon, 3)
};

//---------------------------------------------------------------------------

class Electron: public SortableObject {
 public:

  Float_t PT; // electron transverse momentum
  Float_t Eta; // electron pseudorapidity
  Float_t Phi; // electron azimuthal angle

  Int_t Status; // particle status | hepevt.isthep[number]                                                                                                                                 
  Int_t IsRecoPU; // Assigned by charged hadron subtractor where applicable                                                                                                                
  Int_t IsPU;
  Int_t IsFakeObject;


  Int_t Charge; // electron charge
  Float_t EhadOverEem; // ratio of the hadronic versus electromagnetic energy deposited in the calorimeter

  Float_t IsolationVarDBeta;
  Float_t IsolationVarRhoCorr;
  Float_t TrackIsolationVar;
  Float_t chargedHadronEnergy;
  Float_t neutralEnergy;
  Float_t chargedPUEnergy;
  Float_t allParticleEnergy;


  Float_t TOuter; // Time at calorimeter face

  TRef Particle; // reference to generated particle

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4();

  ClassDef(Electron,3)
};

//---------------------------------------------------------------------------

class Muon: public SortableObject {

 public:

  Float_t PT; // muon transverse momentum
  Float_t Eta; // muon pseudorapidity
  Float_t Phi; // muon azimuthal angle

  Int_t Status; // particle status | hepevt.isthep[number]                                                                                                                                 
  Int_t IsRecoPU; // Assigned by charged hadron subtractor where applicable                                                                                                                
  Int_t IsPU;
  Int_t IsFakeObject;

  Int_t   Charge; // muon charge
  Float_t IsolationVarDBeta;
  Float_t IsolationVarRhoCorr;
  Float_t TrackIsolationVar;
  Float_t chargedHadronEnergy;
  Float_t neutralEnergy;
  Float_t chargedPUEnergy;
  Float_t allParticleEnergy;


  TRef Particle; // reference to generated particle

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4();

  ClassDef(Muon, 2)
};


//---------------------------------------------------------------------------

class IsoTrack: public SortableObject {
 
 public:

  Float_t PT; // IsoTrack transverse momentum
  Float_t Eta; // IsoTrack pseudorapidity
  Float_t Phi; // IsoTrack azimuthal angle

  Int_t   Charge; // IsoTrack charge
  Float_t IsolationVar;
  Float_t TrackIsolationVar;
  Float_t chargedHadronEnergy;
  Float_t neutralEnergy;
  Float_t chargedPUEnergy;
  Float_t allParticleEnergy;

  Int_t   IsEMCand; // Whether it is a Electron/Muon candidate

  TRef Particle; // reference to generated particle

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4();

  ClassDef(IsoTrack, 3)
};

//---------------------------------------------------------------------------

class Jet: public SortableObject {
 
 public:

  Float_t PT; // jet transverse momentum
  Float_t Eta; // jet pseudorapidity
  Float_t Phi; // jet azimuthal angle
  Float_t Mass; // jet invariant mass

  Float_t DeltaEta;  // jet radius in pseudorapidity
  Float_t DeltaPhi;  // jet radius in azimuthal angle

  Float_t Tau1;
  Float_t Tau2;
  Float_t Tau3;

  int NSubJetsTrimmed;

  float TrimmedMass;
  float TrimmedPt;
  float TrimmedEta;
  float TrimmedPhi;

  float TrimmedMassSub1;
  float TrimmedPtSub1;
  float TrimmedEtaSub1;
  float TrimmedPhiSub1;

  float TrimmedMassSub2;
  float TrimmedPtSub2;
  float TrimmedEtaSub2;
  float TrimmedPhiSub2;

  float TrimmedMassSub3;
  float TrimmedPtSub3;
  float TrimmedEtaSub3;
  float TrimmedPhiSub3;

  int NSubJetsPruned;

  float PrunedMass;
  float PrunedPt;
  float PrunedEta;
  float PrunedPhi;

  float PrunedMassSub1;
  float PrunedPtSub1;
  float PrunedEtaSub1;
  float PrunedPhiSub1;

  float PrunedMassSub2;
  float PrunedPtSub2;
  float PrunedEtaSub2;
  float PrunedPhiSub2;

  float PrunedMassSub3;
  float PrunedPtSub3;
  float PrunedEtaSub3;
  float PrunedPhiSub3;

  int NSubJetsSoftDrop;

  float SoftDropMass;
  float SoftDropPt;
  float SoftDropEta;
  float SoftDropPhi;

  float SoftDropMassSub1;
  float SoftDropPtSub1;
  float SoftDropEtaSub1;
  float SoftDropPhiSub1;

  float SoftDropMassSub2;
  float SoftDropPtSub2;
  float SoftDropEtaSub2;
  float SoftDropPhiSub2;

  float SoftDropMassSub3;
  float SoftDropPtSub3;
  float SoftDropEtaSub3;
  float SoftDropPhiSub3;

  Float_t AreaX,AreaY,AreaZ,AreaT; // 4-vector area

  UInt_t BTagAlgo;
  UInt_t BTagDefault;
  UInt_t BTagPhysics;
  UInt_t BTagNearest2;
  UInt_t BTagNearest3;
  UInt_t BTagHeaviest;
  UInt_t BTagHighestPt;

  UInt_t flavourAlgo;
  UInt_t flavourDefault;
  UInt_t flavourPhysics;
  UInt_t flavourNearest2;
  UInt_t flavourNearest3;
  UInt_t flavourHeaviest;
  UInt_t flavourHighestPt;
  
  UInt_t TauTag; // 0 or 1 for a jet that has been tagged as a tau

  Int_t Charge; // tau charge

  Float_t EhadOverEem; // ratio of the hadronic versus electromagnetic energy deposited in the calorimeter

  // pileup jet ID
  float dRMean  ;
  float dR2Mean ;
  float ptD     ;
  float sumPt   ;
  float sumPt2  ;

  float dRMeanEm ;
  float ptDNe    ;
  float sumPtNe  ;
  float nNeutral ;
  float neuEMfrac ;
  float dRMeanNeut ;
  float neuHadfrac ;

  float dRMeanCh ;
  float ptDCh    ;
  float sumPtCh  ;
  float nCharged ;

  float chgEMfrac  ;
  float chgHadfrac ;

  float betaClassic ;
  float betaClassicStar ;
  float beta     ;
  float betaStar ;
  float constituents ;

  float dZ ;
  float d0 ;

  float etaW ;
  float phiW ;
  float jetW ;

  float majW ; 
  float minW ; 
  float dRLeadCent ;
  float dRLead2nd ;

  float ptMean ; 
  float ptRMS ; 
  float pt2A ; 
  float sumChPt ; 
  float sumNePt ;
  float axis2;

  float leadFrac ; 
  float secondFrac ; 
  float thirdFrac ; 
  float fourthFrac ;
  float leadChFrac ; 
  float secondChFrac ; 
  float thirdChFrac ; 
  float fourthChFrac ;
  float leadEmFrac ; 
  float secondEmFrac ; 
  float thirdEmFrac ; 
  float fourthEmFrac ;
  float leadNeutFrac ; 
  float secondNeutFrac ; 
  float thirdNeutFrac ; 
  float fourthNeutFrac ;

  int pileupIDFlagCutBased;

  std::vector<float> FracPt;
  std::vector<float> emFracPt;
  std::vector<float> neutFracPt;
  std::vector<float> chFracPt;

  TRefArray Constituents; // references to constituents
  TRefArray Particles; // references to generated particles

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4();
  TLorentzVector AreaP4();

  ClassDef(Jet,4)
};

//---------------------------------------------------------------------------

class Track: public SortableObject  {

 public:  

  Int_t PID; // HEP ID number
  Int_t Status; // particle status | hepevt.isthep[number]                                                                                             
  Int_t IsRecoPU; // Assigned by charged hadron subtractor where applicable
  Int_t IsPU;
  Int_t Charge; // track charge

  Float_t PT; // track transverse momentum
  Float_t Eta; // track pseudorapidity
  Float_t Phi; // track azimuthal angle
  Float_t Mass;

  Float_t EtaOuter; // track pseudorapidity at the tracker edge
  Float_t PhiOuter; // track azimuthal angle at the tracker edge

  Float_t X; // track vertex position (x component)
  Float_t Y; // track vertex position (y component)
  Float_t Z; // track vertex position (z component)
  Float_t T; // track vertex time: I think this is cheating, but I guess no more than (X,Y,Z)

  Float_t XOuter; // track position (x component) at the tracker edge
  Float_t YOuter; // track position (y component) at the tracker edge
  Float_t ZOuter; // track position (z component) at the tracker edge

  Float_t TOuter; // track time at tracker edge

  TRef Particle; // reference to generated particle

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4();

  ClassDef(Track, 2)
};

//---------------------------------------------------------------------------

class Tower: public SortableObject  {

 public:
  Float_t ET; // calorimeter tower transverse energy
  Float_t Eta; // calorimeter tower pseudorapidity
  Float_t Phi; // calorimeter tower azimuthal angle

  Float_t E; // calorimeter tower energy

  Float_t Eem; // calorimeter tower electromagnetic energy
  Float_t Ehad; // calorimeter tower hadronic energy

  Float_t Edges[4]; // calorimeter tower edges

  Float_t TOuter; // ecal deposit time, averaged by sqrt(EM energy) over all particles, not smeared
  Float_t nTimes;

  TRefArray Particles; // references to generated particles

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4();

  ClassDef(Tower, 2)
};

//---------------------------------------------------------------------------

class Candidate: public SortableObject  {
  friend class DelphesFactory;

 public:

  Candidate();

  // Momentum and position
  TLorentzVector Momentum, Position, Area;

  // Generator level info
  Int_t PID;
  Int_t Status;
  Int_t M1, M2, D1, D2;
  Int_t Charge;
  Int_t Spin;
  Float_t Mass;
  
  // Isolatiom info
  Float_t IsolationVarDBeta;
  Float_t IsolationVarRhoCorr;
  Float_t TrackIsolationVar;
  Float_t chargedHadronEnergy;
  Float_t neutralEnergy;
  Float_t chargedPUEnergy;
  Float_t allParticleEnergy;

  // Flag info
  Int_t IsPU;
  Int_t IsRecoPU; // Assigned by charged hadron subtractor where applicable Reflects decision ncoprorating z vertex resolution
  Int_t IsConstituent;
  Int_t IsEMCand; // Whether it is a Electron/Muon candidate
  Int_t IsFakeObject;

  // BTagging info
  UInt_t BTagAlgo;
  UInt_t BTagDefault;
  UInt_t BTagPhysics;
  UInt_t BTagNearest2;
  UInt_t BTagNearest3;
  UInt_t BTagHeaviest;
  UInt_t BTagHighestPt;

  UInt_t flavourAlgo;
  UInt_t flavourDefault;
  UInt_t flavourPhysics;
  UInt_t flavourNearest2;
  UInt_t flavourNearest3;
  UInt_t flavourHeaviest;
  UInt_t flavourHighestPt;
  

  // Tau tag info
  UInt_t TauTag;

  // Calorimeter info
  Float_t Eem;
  Float_t Ehad;

  // Substructure info
  Float_t Tau1; 
  Float_t Tau2; 
  Float_t Tau3; 

  int NSubJetsTrimmed;

  float TrimmedMass;
  float TrimmedPt;
  float TrimmedEta;
  float TrimmedPhi;

  float TrimmedMassSub1;
  float TrimmedPtSub1;
  float TrimmedEtaSub1;
  float TrimmedPhiSub1;

  float TrimmedMassSub2;
  float TrimmedPtSub2;
  float TrimmedEtaSub2;
  float TrimmedPhiSub2;

  float TrimmedMassSub3;
  float TrimmedPtSub3;
  float TrimmedEtaSub3;
  float TrimmedPhiSub3;

  int NSubJetsPruned;

  float PrunedMass;
  float PrunedPt;
  float PrunedEta;
  float PrunedPhi;

  float PrunedMassSub1;
  float PrunedPtSub1;
  float PrunedEtaSub1;
  float PrunedPhiSub1;

  float PrunedMassSub2;
  float PrunedPtSub2;
  float PrunedEtaSub2;
  float PrunedPhiSub2;

  float PrunedMassSub3;
  float PrunedPtSub3;
  float PrunedEtaSub3;
  float PrunedPhiSub3;

  int NSubJetsSoftDrop;

  float SoftDropMass;
  float SoftDropPt;
  float SoftDropEta;
  float SoftDropPhi;

  float SoftDropMassSub1;
  float SoftDropPtSub1;
  float SoftDropEtaSub1;
  float SoftDropPhiSub1;

  float SoftDropMassSub2;
  float SoftDropPtSub2;
  float SoftDropEtaSub2;
  float SoftDropPhiSub2;

  float SoftDropMassSub3;
  float SoftDropPtSub3;
  float SoftDropEtaSub3;
  float SoftDropPhiSub3;

  Int_t   nTimes;

  // pileup jet ID
  float dRMean  ;
  float dR2Mean ;
  float ptD     ;
  float sumPt   ;
  float sumPt2  ;

  float dRMeanEm ;
  float ptDNe    ;
  float sumPtNe  ;
  float nNeutral ;
  float neuEMfrac ;
  float dRMeanNeut ;
  float neuHadfrac ;

  float dRMeanCh ;
  float ptDCh    ;
  float sumPtCh  ;
  float nCharged ;

  float chgEMfrac  ;
  float chgHadfrac ;

  float betaClassic ;
  float betaClassicStar ;
  float beta     ;
  float betaStar ;
  float constituents ;

  float dZ ;
  float d0 ;

  float etaW ;
  float phiW ;
  float jetW ;

  float majW ; 
  float minW ; 
  float dRLeadCent ;
  float dRLead2nd ;

  float ptMean ; 
  float ptRMS ; 
  float pt2A ; 
  float sumChPt ; 
  float sumNePt ;
  float axis2;

  float leadFrac ; 
  float secondFrac ; 
  float thirdFrac ; 
  float fourthFrac ;
  float leadChFrac ; 
  float secondChFrac ; 
  float thirdChFrac ; 
  float fourthChFrac ;
  float leadEmFrac ; 
  float secondEmFrac ; 
  float thirdEmFrac ; 
  float fourthEmFrac ;
  float leadNeutFrac ; 
  float secondNeutFrac ; 
  float thirdNeutFrac ; 
  float fourthNeutFrac ;

  int pileupIDFlagCutBased;  

  Float_t Edges[4];
  Float_t DeltaEta;
  Float_t DeltaPhi;

  std::vector<float> FracPt; 
  std::vector<float> emFracPt; 
  std::vector<float> neutFracPt;
  std::vector<float> chFracPt; 

  // Towers
  std::vector<std::pair<float,float> > ecal_E_t; 

  // helper methods
  static CompBase *fgCompare; 
  const  CompBase *GetCompare() const { return fgCompare; }

  void       AddCandidate(Candidate *object);
  TObjArray *GetCandidates();

  Bool_t Overlaps(const Candidate *object) const;

  virtual void     Copy(TObject &object) const;
  virtual TObject *Clone(const char *newname = "") const;
  virtual void     Clear(Option_t* option = ""); 

 private:

  DelphesFactory *fFactory; //!
  TObjArray *fArray; //!  
  void SetFactory(DelphesFactory *factory) { fFactory = factory; }

  ClassDef(Candidate, 1)
};

#endif // DelphesClasses_h


