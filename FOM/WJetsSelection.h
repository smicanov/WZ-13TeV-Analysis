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

  TH1D* hWMt[6];
  TH1D* hFakeMt[6];

  TH1D* hDeltaPhiWBosonWl[6];
  TH1D* hDeltaRWBosonWl[6];
  TH1D* hDeltaPhiWBosonFakel[6];
  TH1D* hDeltaRWBosonFakel[6];

  TH1D* hDeltaPhiFakeBosonWl[6];
  TH1D* hDeltaRFakeBosonWl[6];
  TH1D* hDeltaPhiFakeBosonFakel[6];
  TH1D* hDeltaRFakeBosonFakel[6];

  TH1D* h2LMass[6];

  TH1D* hGoodJets[6];
  TH1D* hDRminGoodJetWl[6];
  TH1D* hDRminGoodJetFakel[6];

  TH1D* hGoodJetsCut[6];
  TH1D* hDRminGoodJetCutWl[6];
  TH1D* hDRminGoodJetCutFakel[6];

  TH1D* hGoodJetsLeadCut[6];
  TH1D* hDRminGoodJetWlLeadCut[6];
  TH1D* hDRminGoodJetFakelLeadCut[6];

  TH1D* hGoodJetsCutLeadCut[6];
  TH1D* hDRminGoodJetCutWlLeadCut[6];
  TH1D* hDRminGoodJetCutFakelLeadCut[6];

  TH1D* hLeadJetPt[6];
  TH1D* hLeadJetEta[6];
  TH1D* hLeadJetPhi[6];
  TH1D* hLeadJetEt[6];

};

#endif

