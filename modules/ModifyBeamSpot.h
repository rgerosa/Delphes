#ifndef ModifyBeamSpot_h
#define ModifyBeamSpot_h

/** \class ModifyBeamSpot
 *
 *  \author S. Zenz
 *
 */

#include "classes/DelphesModule.h"
#include "modules/simpleVariableCollector.h"

class TIterator;
class TObjArray;
class DelphesFormula;

class ModifyBeamSpot: public DelphesModule
{
public:

  ModifyBeamSpot();
  ~ModifyBeamSpot();

  void Init();
  void Process();
  void Finish();

private:

  DelphesFormula *fFormula; //!

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!
  
  TObjArray *fOutputArray; //!

  Double_t fOutputBSX, fOutputBSY, fOutputBSZ ;
  Double_t fZVertexSpread;
  Double_t fTVertexSpread;
  Double_t currentZ, currentT;

  // Store Z of PV
  TObjArray *fPVOutputArray; //!

  simpleVariableCollector fDebugOutputCollector ;    
  int fEventCounter ;                                                                                                                                             

  ClassDef(ModifyBeamSpot, 1)
};

#endif
