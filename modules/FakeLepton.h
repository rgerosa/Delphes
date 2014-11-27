#ifndef FakeLepton_h
#define FakeLepton_h

/** \class FakeLepton
 *  Selects candidates from the InputArray according to the efficiency formula.
 *  $Date: 2013-02-12 14:57:44 +0100 (Tue, 12 Feb 2013) $
 *  $Revision: 905 $
 *  \author P. Demin - UCL, Louvain-la-Neuve
*/

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;
class DelphesFormula;

class FakeLepton: public DelphesModule {

public:

  FakeLepton();
  ~FakeLepton();

  void Init();
  void Process();
  void Finish();

private:

  DelphesFormula *fFormula; //!

  TIterator *fItInputFakeArray; //!

  const TObjArray *fInputFakeArray; //!

  TObjArray *fInputLeptonArray; //!

  ClassDef(FakeLepton,1)
};

#endif
