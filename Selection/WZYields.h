#ifndef WZYields_h
#define WZYields_h

#include "Generic.h"


class WZYields : public Generic
{

public:

  WZYields(WZEvent*, TFile* outputFile);

	unsigned int GetNAnalyzed() { return nAnalyzedEvents; }
	unsigned int GetNSelected() { return nSelectedEvents; }

  void Init();

  void Analysis();

  void Finish();


protected:

  long unsigned int minRun;
  long unsigned int maxRun;
  vector<long unsigned int> runNumber;

  unsigned int nAnalyzedEvents;
  unsigned int nSelectedEvents;

  unsigned int yieldsByChannelTrigger[6];
  unsigned int yieldsByChannelPreselection[6];
  unsigned int yieldsByChannelZSelection[6];
  unsigned int yieldsByChannelWSelection[6];
  unsigned int yieldsByChannelFullSelection[6];


  std::ofstream eventLists1[4];
  std::ofstream eventLists2[4];
  std::ofstream eventLists3[4];
  std::ofstream eventLists4[4];


  TH1D* hZmass[5];
  TH1D* hZpt[5];
  TH1D* hMET[5];
  TH1D* hMt[5];
  TH1D* hZl1pt[5];
  TH1D* hZl2pt[5];
  TH1D* hWlpt[5];
  TH1D* hNJets[5];
  TH1D* hNJetsNoMuIso[5];
  TH1D* hNJetsNoEleIso[5];
  TH1D* hNJetsNoIso[5];
  TH1D* h3LMass[5];
  TH1D* hDeltaR[5];
  TH1D* hDeltaRMin[5];

};

#endif

