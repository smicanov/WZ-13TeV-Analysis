#ifndef WZJets_h
#define WZJets_h

#include "Generic.h"


class WZJets : public Generic
{

public:

  WZJets(WZEvent*, TFile* outputFile);

	unsigned int GetNAnalyzed() { return nAnalyzedEvents; }
	unsigned int GetNSelected() { return nSelectedEvents; }
	unsigned int GetNGoodJets() { return nGoodJets; }

  void Init();

  void Analysis();

  void Finish();


protected:

  unsigned int nAnalyzedEvents;
  unsigned int nSelectedEvents;
  unsigned int nGoodJets;

  TH1D* hNJets;
  TH1D* hNGoodJets;

  TH1D* hNJetsWMu;
  TH1D* hNJetsWEle;
  TH1D* hNJetsWL;

  TH2D* hNJetsVsNMu;
  TH2D* hNJetsVsNEle;
  TH2D* hNJetsVsNL;

  TH1D* hDeltaRMu;
  TH1D* hDeltaREle;
  TH1D* hDeltaRL;

};

#endif
