#ifndef PileUpJetID_h
#define PileUpJetID_h

#include "classes/DelphesModule.h"
#include "TMatrixDSym.h"

#include <deque>

class TObjArray;
class DelphesFormula;

class PileUpJetID: public DelphesModule {

 public:

  PileUpJetID();
  ~PileUpJetID();

  void Init();
  void Process();
  void Finish();
  void assign(const std::vector<float> &, float &, float &, float &, float &);
  std::pair<int,int> getJetIdKey(float jetPt, float jetEta);
  int computeCutIDflag(float,float,float,float,float);

 private:

  Double_t fJetPTMin;
  Double_t fParameterR;

  Int_t    fUseConstituents; 

  std::vector<double> fCones ;

  std::vector<double> Pt010_Tight_betaStar ;
  std::vector<double> Pt1020_Tight_betaStar ;
  std::vector<double> Pt2030_Tight_betaStar ;
  std::vector<double> Pt3050_Tight_betaStar ;

  std::vector<double> Pt010_Medium_betaStar ;
  std::vector<double> Pt1020_Medium_betaStar ;
  std::vector<double> Pt2030_Medium_betaStar ;
  std::vector<double> Pt3050_Medium_betaStar ;

  std::vector<double> Pt010_Loose_betaStar ;
  std::vector<double> Pt1020_Loose_betaStar ;
  std::vector<double> Pt2030_Loose_betaStar ;
  std::vector<double> Pt3050_Loose_betaStar;

  std::vector<double> Pt010_Tight_RMS ;
  std::vector<double> Pt1020_Tight_RMS ;
  std::vector<double> Pt2030_Tight_RMS ;
  std::vector<double> Pt3050_Tight_RMS ;

  std::vector<double> Pt010_Medium_RMS ;
  std::vector<double> Pt1020_Medium_RMS ;
  std::vector<double> Pt2030_Medium_RMS ;
  std::vector<double> Pt3050_Medium_RMS ;

  std::vector<double> Pt010_Loose_RMS ;
  std::vector<double> Pt1020_Loose_RMS ;
  std::vector<double> Pt2030_Loose_RMS ;
  std::vector<double> Pt3050_Loose_RMS;

  Float_t pileUpIDCut_betaStar [3][4][4];
  Float_t pileUpIDCut_RMS [3][4][4];
 
  TIterator *fItJetInputArray; //!
  TIterator *fItTrackInputArray; // SCZ
  TIterator *fItNeutralInputArray; // SCZ
  TIterator *fItPVInputArray; //!                                                                                                                                                     

  const TObjArray *fJetInputArray; //!
  const TObjArray *fTrackInputArray; // SCZ
  const TObjArray *fNeutralInputArray; 
  const TObjArray *fPVInputArray; //!                                                                                                                                                 

  TObjArray *fOutputArray; //!
  TObjArray *fNeutralsInPassingJets; // SCZ

  ClassDef(PileUpJetID, 2)
};

#endif
