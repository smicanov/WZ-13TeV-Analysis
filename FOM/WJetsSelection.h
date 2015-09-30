#ifndef WJetsSelection_h
#define WJetsSelection_h

#include "Generic.h"


class WJetsSelection : public Generic
{

public:

  WJetsSelection(Event* event, TFile* outputFile);

  unsigned int GetNAnalyzed() { return nAnalyzedEvents; }
  unsigned int GetNSelected() { return nSelectedEvents; }

  void Init();
  void Analysis();
  void Finish();


protected:

  unsigned int nAnalyzedEvents;
  unsigned int nSelectedEvents;

  unsigned int yieldsByChannelPreselection[6];
  unsigned int yieldsByChannelSSSelection[6];
  unsigned int yieldsByChannelOFSelection[6];
  unsigned int yieldsByChannelFullSelection[6];

/*
  std::ofstream eventLists1[4];
  std::ofstream eventLists2[4];
  std::ofstream eventLists3[4];
  std::ofstream eventLists4[4];
*/

  TH1D* hWPt[6];
  TH1D* hWEta[6];
  TH1D* hWPhi[6];
  TH1D* hWRelIso[6];

  TH1D* hFakePt[6];
  TH1D* hFakeEta[6];
  TH1D* hFakePhi[6];
  TH1D* hFakeRelIso[6];

  TH1D* hDeltaPhiWFake[6];
  TH1D* hDeltaRWFake[6];

  TH1D* hMET[6];
  TH1D* hMETPhi[6];

  TH1D* hDeltaRWMET[6];
  TH1D* hDeltaPhiWMET[6];

  TH1D* hDeltaRFakeMET[6];
  TH1D* hDeltaPhiFakeMET[6];

  TH1D* hMt[6];
  TH1D* h2LMass[6];

  TH1D* hGoodJets[6];
  TH1D* hDRminJetWl[6];
  TH1D* hDRminJetFakel[6];
};

#endif

