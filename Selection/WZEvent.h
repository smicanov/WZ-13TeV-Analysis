#ifndef WZEvent_h
#define WZEvent_h

#define WZBASECLASS EventTree_V07_04_14_00

#include "EventTree_V07_04_14_00.h"
#include "Particles.h"

#include <vector>


enum FinalState
{
  undefined,   // 0
  eee,         // 1
  eem,         // 2
  mme,         // 3
  mmm          // 4
};


enum SelectionLevel
{
  Undefined,            // 0
  FailsPreselection,    // 1
  Preselection,         // 2
  ZSelection,           // 3
  WSelection,           // 4
  FullSelection         // 5
};


enum SelectionType
{
  Nominal,      // 0
  MatrixMethod  // 1
};


class WZEvent : public WZBASECLASS
{

  friend class WZSelection;
  friend class WZJets;

public:

  WZEvent(TTree* tree);

  vector<bool> GetHLT25ns() { return fHLT25ns; }

  bool PassesTriggerAll();
  bool PassesTriggerDoubleEG();
  bool PassesTriggerDoubleMuon();
  bool PassesTriggerMuonEG();

  bool PassesPreselection(SelectionType type = Nominal);
  bool PassesZSelection(SelectionType type = Nominal);
  bool PassesWSelection(SelectionType type = Nominal);
  bool PassesFullSelection(SelectionType type = Nominal);

  FinalState GetFinalState() { return fFinalState; }
  SelectionLevel GetSelectionLevel() { return fSelectionLevel; }

  Lepton* GetWLepton() { return fLeptons.at(fWLeptonIndex); }
  pair<Lepton*, Lepton*> GetZLeptons();

  void Read();

  void Dump(std::ostream& out, int verbosity = 0);


protected:

  void Clear();

  vector<bool> fHLT25ns;
// 0 - HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v        - bit  7
// 1 - HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v              - bit 20
// 2 - HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v            - bit 21
// 3 - HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v   - bit 50
// 4 - HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v  - bit 51

  FinalState fFinalState;
  SelectionLevel fSelectionLevel;

  vector<Lepton*> fLeptons;
  vector<unsigned int> fTightLeptonsIndex;
  pair<int, int> fZLeptonsIndex;
  int fWLeptonIndex;

};

#endif

