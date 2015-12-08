#include "WZControlRegionsYields.h"
#include "Constants.h"
#include "TLatex.h"

#include <cmath>
#include <sstream>
#include <string>
#include <boost/lexical_cast.hpp>


using namespace std;


WZControlRegionsYields::WZControlRegionsYields(WZEvent* e, TFile* outputFile)
 : Generic(e, outputFile)
{

}


void WZControlRegionsYields::Init()
{
  nAnalyzedEvents = 0;
  nSelectedEvents = 0;

  for (unsigned int i = 0; i < 9; i++) {
    for (unsigned int j = 0; j < 5; j++)
      yieldsByChannelControlRegions[i][j] = 0;
  }
}


void WZControlRegionsYields::Analysis()
{
  nAnalyzedEvents++;

  if (fWZEvent->PassesTriggerMuonEG() && fWZEvent->PassesFullSelection(MatrixMethod) && fWZEvent->PassesMETFilters())
    yieldsByChannelControlRegions[fWZEvent->GetControlRegion()][fWZEvent->GetFinalState()]++;

  if (fWZEvent->PassesTriggerMuonEG() && fWZEvent->PassesFullSelection(MatrixMethod) && !(fWZEvent->PassesMETFilters())) {
    cout << "Failed MET Filter !!!\t-\t" << fWZEvent->run << ":" << fWZEvent->lumis << ":" << fWZEvent->event << endl;
    nSelectedEvents++;
  }

}


void WZControlRegionsYields::Finish()
{
  cout << endl << "Done." << endl;

  cout << "Analyzed events : " << GetNAnalyzed() << "\n"
       << "Failed MET Filters : "<< GetNSelected() << "\n\n";


  cout << "CHANNEL\tNone\tPPP\tPPF\tPFP\tFPP\tFFP\tFPF\tPFF\tFFF" << "\n";
  for (int j = 0; j < 5; j++) {
    cout << j;
    for (int i = 0; i < 9; i++)
      cout << "\t" << yieldsByChannelControlRegions[i][j];
    cout << "\n";
  }

  cout << endl;
}

