#ifndef Event_h
#define Event_h

#define BASECLASS EventTree_V07_04_09_01

#include "EventTree_V07_04_09_01.h"
#include "Particles.h"

#include <vector>
#include <iostream>


enum FinalState
{
  undefined,   // 0
  WEleFakeMu,  // 1
  WMuFakeEle   // 2
};


enum SelectionLevel
{
  Undefined,            // 0
  FailsPreselection,    // 1
  Preselection,         // 2
  SSOFSelection,        // 3
  FullSelection         // 4
};


class Event : public BASECLASS
{

  friend class WJetsSelection;

public:

  Event(TTree* tree);

  bool PassesPreselection();
  bool PassesSSOFSelection();
  bool PassesFullSelection();

  FinalState GetFinalState() { return fFinalState; }
  SelectionLevel GetSelectionLevel() { return fSelectionLevel; }

  void Read();
  void Dump(std::ostream& out, int verbosity=0);


protected:

  void Clear();

  FinalState fFinalState;
  SelectionLevel fSelectionLevel;

  vector<Lepton*> fLeptons;
  pair<int, int> fCandidateLeptonIndex;  // (W lepton, fake (FOMLoose) lepton)

};

#endif

