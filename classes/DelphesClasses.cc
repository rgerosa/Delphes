/**
 *  Definition of classes to be stored in the root tree.
 *  Function CompareXYZ sorts objects by the variable XYZ that MUST be
 *  present in the data members of the root tree class of the branch.
 *  $Date: 2008-06-04 13:57:24 $
 *  $Revision: 1.1 $
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"

#include "classes/DelphesFactory.h"
#include "classes/SortableObject.h"

CompBase *GenParticle::fgCompare = 0;
CompBase *Photon::fgCompare   = CompPT<Photon>::Instance();
CompBase *Electron::fgCompare = CompPT<Electron>::Instance();
CompBase *Muon::fgCompare     = CompPT<Muon>::Instance();
CompBase *IsoTrack::fgCompare = CompPT<IsoTrack>::Instance();
CompBase *Jet::fgCompare      = CompPT<Jet>::Instance();
CompBase *Track::fgCompare    = CompPT<Track>::Instance();
CompBase *Tower::fgCompare    = CompE<Tower>::Instance();
CompBase *Candidate::fgCompare = CompMomentumPt<Candidate>::Instance();

//------------------------------------------------------------------------------
TLorentzVector GenParticle::P4(){
  TLorentzVector vec;
  vec.SetPxPyPzE(Px, Py, Pz, E);
  return vec;
}

//------------------------------------------------------------------------------
TLorentzVector Photon::P4(){
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------
TLorentzVector Electron::P4(){
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------
TLorentzVector Muon::P4(){
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------
TLorentzVector IsoTrack::P4(){
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------
TLorentzVector Jet::P4(){
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, Mass);
  return vec;
}

TLorentzVector Jet::AreaP4(){
  TLorentzVector vec(AreaX,AreaY,AreaZ,AreaT);
  return vec;
}

//------------------------------------------------------------------------------
TLorentzVector Track::P4(){
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------
TLorentzVector Tower::P4(){
  TLorentzVector vec;
  vec.SetPtEtaPhiM(ET, Eta, Phi, 0.0);
  return vec;
}


//------------------------------------------------------------------------------
Candidate::Candidate() :
  PID(0), Status(0), M1(-1), M2(-1), D1(-1), D2(-1),
  IsolationVar(0), TrackIsolationVar(0),
  Charge(0), Mass(0.0),
  IsPU(0), IsRecoPU(0), IsEMCand(0), IsConstituent(0),
  BTag(0), TauTag(0), Eem(0.0), Ehad(0.0),
  Tau1(0), Tau2(0), Tau3(0),
  NSubJetsTrimmed(0),
  TrimmedMass(0), TrimmedPt(0), TrimmedEta(0), TrimmedPhi(0),
  TrimmedMassSub1(0),TrimmedPtSub1(0),TrimmedEtaSub1(0),TrimmedPhiSub1(0),
  TrimmedMassSub2(0),TrimmedPtSub2(0),TrimmedEtaSub2(0),TrimmedPhiSub2(0),
  TrimmedMassSub3(0),TrimmedPtSub3(0),TrimmedEtaSub3(0),TrimmedPhiSub3(0),
  PrunedMass(0),PrunedPt(0),PrunedEta(0),PrunedPhi(0),
  PrunedMassSub1(0),PrunedPtSub1(0),PrunedEtaSub1(0),PrunedPhiSub1(0),
  PrunedMassSub2(0),PrunedPtSub2(0),PrunedEtaSub2(0),PrunedPhiSub2(0),
  PrunedMassSub3(0),PrunedPtSub3(0),PrunedEtaSub3(0),PrunedPhiSub3(0),
  NSubJetsSoftDrop(0),
  SoftDropMass(0),SoftDropPt(0),SoftDropEta(0),SoftDropPhi(0),
  SoftDropMassSub1(0),SoftDropPtSub1(0),SoftDropEtaSub1(0),SoftDropPhiSub1(0),
  SoftDropMassSub2(0),SoftDropPtSub2(0),SoftDropEtaSub2(0),SoftDropPhiSub2(0),
  SoftDropMassSub3(0),SoftDropPtSub3(0),SoftDropEtaSub3(0),SoftDropPhiSub3(0),
  Beta(-999.), BetaStar(-999.), 
  MeanSqDeltaR(-999.), PTD(-999.),
  NCharged(-999), NNeutrals(-999),
  DeltaEta(0.0), DeltaPhi(0.0),
  Momentum(0.0, 0.0, 0.0, 0.0),
  Position(0.0, 0.0, 0.0, 0.0),
  Area(0.0, 0.0, 0.0, 0.0),
  fFactory(0),
  fArray(0){
  Edges[0] = 0.0;
  Edges[1] = 0.0;
  Edges[2] = 0.0;
  Edges[3] = 0.0;
  for (int i = 0 ; i < 5 ; i++) {
    FracPt[i] = -999.;
  }
}

void Candidate::AddCandidate(Candidate *object){
  if(!fArray) fArray = fFactory->NewArray();
  fArray->Add(object);
}

TObjArray *Candidate::GetCandidates(){
  if(!fArray) fArray = fFactory->NewArray();
  return fArray;
}

Bool_t Candidate::Overlaps(const Candidate *object) const{
  const Candidate *candidate;

  if(object->GetUniqueID() == GetUniqueID()) return kTRUE;

  if(fArray){
    TIter it(fArray);
    while((candidate = static_cast<Candidate *>(it.Next()))){
      if(candidate->Overlaps(object)) return kTRUE;
    }
  }

  if(object->fArray){
    TIter it(object->fArray);
    while((candidate = static_cast<Candidate *>(it.Next()))){
      if(candidate->Overlaps(this)) return kTRUE;
    }
  }

  return kFALSE;
}


TObject *Candidate::Clone(const char *newname) const{
  Candidate *object = fFactory->NewCandidate();
  Copy(*object);
  return object;
}

void Candidate::Copy(TObject &obj) const{
  Candidate &object = static_cast<Candidate &>(obj);
  Candidate *candidate;

  object.PID = PID;
  object.Status = Status;
  object.IsolationVar = IsolationVar;
  object.TrackIsolationVar = TrackIsolationVar;
  object.M1 = M1;
  object.M2 = M2;
  object.D1 = D1;
  object.D2 = D2;
  object.Charge = Charge;
  object.Mass = Mass;
  object.IsPU = IsPU;
  object.IsRecoPU = IsRecoPU;
  object.IsEMCand = IsEMCand;
  object.IsConstituent = IsConstituent;
  object.BTag = BTag;
  object.TauTag = TauTag;
  object.Eem = Eem;
  object.Ehad = Ehad;
  object.Edges[0] = Edges[0];
  object.Edges[1] = Edges[1];
  object.Edges[2] = Edges[2];
  object.Edges[3] = Edges[3];
  object.DeltaEta = DeltaEta;
  object.DeltaPhi = DeltaPhi;
  object.Momentum = Momentum;
  object.Position = Position;
  object.Area = Area;
  object.Beta = Beta;
  object.BetaStar = BetaStar;
  object.MeanSqDeltaR = MeanSqDeltaR;
  object.PTD = PTD;
  object.NCharged = NCharged; 
  object.NNeutrals = NNeutrals;

  object.Tau1 = Tau1 ;
  object.Tau2 = Tau2 ;
  object.Tau3 = Tau3 ;

  object.NSubJetsTrimmed = NSubJetsTrimmed;

  object.TrimmedMass = TrimmedMass; 
  object.TrimmedPt   = TrimmedPt;
  object.TrimmedEta  = TrimmedEta;
  object.TrimmedPhi  = TrimmedPhi;

  object.TrimmedMassSub1 = TrimmedMassSub1;
  object.TrimmedPtSub1   = TrimmedPtSub1;
  object.TrimmedEtaSub1  = TrimmedEtaSub1;
  object.TrimmedPhiSub1  = TrimmedPhiSub1;

  object.TrimmedMassSub2 = TrimmedMassSub2;
  object.TrimmedPtSub2   = TrimmedPtSub2;
  object.TrimmedEtaSub2  = TrimmedEtaSub2;
  object.TrimmedPhiSub2  = TrimmedPhiSub2;

  object.TrimmedMassSub3 = TrimmedMassSub3;
  object.TrimmedPtSub3   = TrimmedPtSub3;
  object.TrimmedEtaSub3  = TrimmedEtaSub3;
  object.TrimmedPhiSub3  = TrimmedPhiSub3;

  object.NSubJetsPruned = NSubJetsPruned;

  object.PrunedMass = PrunedMass; 
  object.PrunedPt   = PrunedPt;
  object.PrunedEta  = PrunedEta;
  object.PrunedPhi  = PrunedPhi;

  object.PrunedMassSub1 = PrunedMassSub1;
  object.PrunedPtSub1   = PrunedPtSub1;
  object.PrunedEtaSub1  = PrunedEtaSub1;
  object.PrunedPhiSub1  = PrunedPhiSub1;

  object.PrunedMassSub2 = PrunedMassSub2;
  object.PrunedPtSub2   = PrunedPtSub2;
  object.PrunedEtaSub2  = PrunedEtaSub2;
  object.PrunedPhiSub2  = PrunedPhiSub2;

  object.PrunedMassSub3 = PrunedMassSub3;
  object.PrunedPtSub3   = PrunedPtSub3;
  object.PrunedEtaSub3  = PrunedEtaSub3;
  object.PrunedPhiSub3  = PrunedPhiSub3;

  object.NSubJetsSoftDrop = NSubJetsSoftDrop;

  object.SoftDropMass = SoftDropMass; 
  object.SoftDropPt   = SoftDropPt;
  object.SoftDropEta  = SoftDropEta;
  object.SoftDropPhi  = SoftDropPhi;

  object.SoftDropMassSub1 = SoftDropMassSub1;
  object.SoftDropPtSub1   = SoftDropPtSub1;
  object.SoftDropEtaSub1  = SoftDropEtaSub1;
  object.SoftDropPhiSub1  = SoftDropPhiSub1;

  object.SoftDropMassSub2 = SoftDropMassSub2;
  object.SoftDropPtSub2   = SoftDropPtSub2;
  object.SoftDropEtaSub2  = SoftDropEtaSub2;
  object.SoftDropPhiSub2  = SoftDropPhiSub2;

  object.SoftDropMassSub3 = SoftDropMassSub3;
  object.SoftDropPtSub3   = SoftDropPtSub3;
  object.SoftDropEtaSub3  = SoftDropEtaSub3;
  object.SoftDropPhiSub3  = SoftDropPhiSub3;

  for (int i = 0 ; i < 5 ; i++) {
    object.FracPt[i] = FracPt[i];
  }

  // Copy cluster timing info
  copy(ecal_E_t.begin(),ecal_E_t.end(),back_inserter(object.ecal_E_t));

  object.fFactory = fFactory;
  object.fArray = 0;

  if(fArray && fArray->GetEntriesFast() > 0){
    TIter itArray(fArray);
    TObjArray *array = object.GetCandidates();
    while((candidate = static_cast<Candidate *>(itArray.Next()))){
      array->Add(candidate);
    }
  }
}

void Candidate::Clear(Option_t* option){
  SetUniqueID(0);
  ResetBit(kIsReferenced);
  PID = 0;
  Status = 0;
  IsolationVar =0.;
  TrackIsolationVar =0.;
  M1 = -1; M2 = -1; D1 = -1; D2 = -1;
  Charge = 0;
  Mass = 0.0;
  IsPU = 0;
  IsRecoPU = 0;
  IsConstituent = 0;
  IsEMCand = 0;
  BTag = 0;
  TauTag = 0;
  Eem = 0.0;
  Ehad = 0.0;
  Edges[0] = 0.0;
  Edges[1] = 0.0;
  Edges[2] = 0.0;
  Edges[3] = 0.0;
  DeltaEta = 0.0;
  DeltaPhi = 0.0;
  Momentum.SetXYZT(0.0, 0.0, 0.0, 0.0);
  Position.SetXYZT(0.0, 0.0, 0.0, 0.0);
  Area.SetXYZT(0.0, 0.0, 0.0, 0.0);
  Beta = -999.;
  BetaStar = -999.;
  MeanSqDeltaR = -999.;
  PTD = -999.;
  NCharged = -1;
  NNeutrals = -1;
  fArray = 0;
  ecal_E_t.clear();

  Tau1 = -999;
  Tau2 = -999;
  Tau3 = -999;

  NSubJetsTrimmed=-999;
  NSubJetsPruned=-999;
  NSubJetsSoftDrop=-999;

  TrimmedMass = -999;
  TrimmedEta = -999;
  TrimmedPhi = -999;
  TrimmedPt  = -999;

  TrimmedMassSub1 = -999;
  TrimmedEtaSub1 = -999;
  TrimmedPhiSub1 = -999;
  TrimmedPtSub1  = -999;

  TrimmedMassSub2 = -999;
  TrimmedEtaSub2 = -999;
  TrimmedPhiSub2 = -999;
  TrimmedPtSub2  = -999;

  TrimmedMassSub3 = -999;
  TrimmedEtaSub3 = -999;
  TrimmedPhiSub3 = -999;
  TrimmedPtSub3  = -999;

  PrunedMass = -999;
  PrunedEta = -999;
  PrunedPhi = -999;
  PrunedPt  = -999;

  PrunedMassSub1 = -999;
  PrunedEtaSub1 = -999;
  PrunedPhiSub1 = -999;
  PrunedPtSub1  = -999;

  PrunedMassSub2 = -999;
  PrunedEtaSub2 = -999;
  PrunedPhiSub2 = -999;
  PrunedPtSub2  = -999;

  PrunedMassSub3 = -999;
  PrunedEtaSub3 = -999;
  PrunedPhiSub3 = -999;
  PrunedPtSub3  = -999;

  SoftDropMass = -999;
  SoftDropEta = -999;
  SoftDropPhi = -999;
  SoftDropPt  = -999;

  SoftDropMassSub1 = -999;
  SoftDropEtaSub1 = -999;
  SoftDropPhiSub1 = -999;
  SoftDropPtSub1  = -999;

  SoftDropMassSub2 = -999;
  SoftDropEtaSub2 = -999;
  SoftDropPhiSub2 = -999;
  SoftDropPtSub2  = -999;

  SoftDropMassSub3 = -999;
  SoftDropEtaSub3 = -999;
  SoftDropPhiSub3 = -999;
  SoftDropPtSub3  = -999;

}

