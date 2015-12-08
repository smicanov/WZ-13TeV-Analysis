#ifndef WZControlRegionsYields_h
#define WZControlRegionsYields_h

#include "Generic.h"


class WZControlRegionsYields : public Generic
{

public:

  WZControlRegionsYields(WZEvent*, TFile* outputFile);

	unsigned int GetNAnalyzed() { return nAnalyzedEvents; }
	unsigned int GetNSelected() { return nSelectedEvents; }

  void Init();

  void Analysis();

  void Finish();


protected:

  unsigned int nAnalyzedEvents;
  unsigned int nSelectedEvents;

  unsigned int yieldsByChannelControlRegions[9][5];

};

#endif

