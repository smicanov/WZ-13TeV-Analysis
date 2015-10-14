#ifndef Event_h
#define Event_h

#define BASECLASS EventTree_V07_04_09_01

#include "EventTree_V07_04_09_01.h"
#include "Particles.h"

#include <vector>
#include <iostream>


enum FinalState
{
  undefined,    // 0
  WEleFakeEle,  // 1
  WEleFakeMu,   // 2
  WMuFakeEle,   // 3
  WMuFakeMu     // 4
};


enum SelectionLevel
{
  Undefined,            // 0
  FailsPreselection,    // 1
  Preselection,         // 2
  SSSelection,          // 3
  OFSelection,          // 4
  FullSelection         // 5
};


class Event : public BASECLASS
{

  friend class WJetsSelection;

public:

  Event(TTree* tree);

  bool PassesPreselection();
  bool PassesSSSelection();
  bool PassesOFSelection();
  bool PassesFullSelection();
  bool PassesTight();
  bool PassesTrigger();

  FinalState GetFinalState() { return fFinalState; }
  SelectionLevel GetSelectionLevel() { return fSelectionLevel; }

  void Read();
  void Dump(std::ostream& out, int verbosity = 0);


protected:

  void Clear();

  FinalState fFinalState;
  SelectionLevel fSelectionLevel;

  vector<Lepton*> fLeptons;
  vector<Jet*> fJets;
  vector<unsigned int> fGoodJetsIndex;
  pair<int, int> fCandidateLeptonIndex;  // (W lepton, fake (FOMLoose) lepton)

};

#endif

