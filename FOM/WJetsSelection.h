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

  unsigned int yieldsByChannelPreselection[4];
  unsigned int yieldsByChannelSSOFSelection[4];
  unsigned int yieldsByChannelFullSelection[4];

/*
  std::ofstream eventLists2[2];
  std::ofstream eventLists3[2];
  std::ofstream eventLists4[2];
*/

  TH1D* hWPt[4];
  TH1D* hFakePt[4];
  TH1D* hDeltaRWFake[4];

  TH1D* hMET[4];
  TH1D* hMt[4];
  TH1D* h2LMass[4];

};

#endif

