#include "WJetsSelection.h"
#include "Constants.h"
#include "TLatex.h"

#include <cmath>
#include <sstream>
#include <string>
#include <boost/lexical_cast.hpp>


using namespace std;


WJetsSelection::WJetsSelection(Event* event, TFile* outputFile)
 : Generic(event, outputFile)
{

}


void WJetsSelection::Init()
{
  nAnalyzedEvents = 0;
  nSelectedEvents = 0;

  for (unsigned int i = 0; i <= 3; i++) {
    hWPt[i] = bookTH1D(("hWPt_" + boost::lexical_cast<string>(i)).c_str(),
                        "W lepton pt", 100, 0, 200);
    hFakePt[i] = bookTH1D(("hFakePt_" + boost::lexical_cast<string>(i)).c_str(),
                           "Fake lepton pt", 100, 0, 200);
    hDeltaRWFake[i] = bookTH1D(("hDeltaRWFake_" + boost::lexical_cast<string>(i)).c_str(),
                          "#DeltaR (W, Fake)", 100, 0, 5);

    hMET[i] = bookTH1D(("hMET_" + boost::lexical_cast<string>(i)).c_str(),
                       "Missing E_{trans}", 100, 0, 200);
    hMt[i] = bookTH1D(("hMt_" + boost::lexical_cast<string>(i)).c_str(),
                       "M_{trans}", 100, 0, 200);
    h2LMass[i] = bookTH1D(("h2LMass_" + boost::lexical_cast<string>(i)).c_str(),
                          "mass(2l)", 100, 0, 200);
  }


  for (int i = 0; i <= 3; i++) {
    yieldsByChannelPreselection[i] = 0;
    yieldsByChannelSSOFSelection[i] = 0;
    yieldsByChannelFullSelection[i] = 0;
  }

/*
// Setup selected event lists 
  for (int i = 1; i <= 2; i++) {
    ostringstream outputFileName2;
    outputFileName2 << "/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/fom/test/WJets_Preselection_" << i << ".txt";
    cout << "File name : " << outputFileName2.str() << endl;
    eventLists2[i-1].open(outputFileName2.str().c_str());

    ostringstream outputFileName3;
    outputFileName3 << "/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/fom/test/WJets_SSOFSelection_" << i << ".txt";
    cout << "File name : " << outputFileName3.str() << endl;
    eventLists3[i-1].open(outputFileName3.str().c_str());

    ostringstream outputFileName4;
    outputFileName4 << "/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/fom/test/WJets_FullSelection_" << i << ".txt";
    cout << "File name : " << outputFileName4.str() << endl;
    eventLists4[i-1].open(outputFileName4.str().c_str());
  }
*/
}


void WJetsSelection::Analysis()
{
  nAnalyzedEvents++;

  if (fEvent->PassesFullSelection()) {
    yieldsByChannelFullSelection[fEvent->GetFinalState()]++;
    yieldsByChannelFullSelection[3]++;
//    fEvent->DumpEvent(eventLists4[fEvent->GetFinalState()-1], 1);
  }

  if (fEvent->PassesSSOFSelection()) {
    yieldsByChannelSSOFSelection[fEvent->GetFinalState()]++;
    yieldsByChannelSSOFSelection[3]++;
//    fEvent->DumpEvent(eventLists3[fEvent->GetFinalState()-1], 1);
  }

  if (fEvent->PassesPreselection()) {
    yieldsByChannelPreselection[fEvent->GetFinalState()]++;
    yieldsByChannelPreselection[3]++;
//    fEvent->DumpEvent(eventLists2[fEvent->GetFinalState()-1], 1);
  }

  if (!(fEvent->PassesPreselection()))  return;
//  if (!(fEvent->PassesOSSFSelection()))  return;
//  if (!(fEvent->PassesFullSelection()))  return;

  nSelectedEvents++;

  Lepton* wLepton = fEvent->fLeptons.at(fEvent->fCandidateLeptonIndex.first);
  Lepton* fakeLepton = fEvent->fLeptons.at(fEvent->fCandidateLeptonIndex.second);

  const double wPt = wLepton->Pt();
  const double fakePt = fakeLepton->Pt();
  const double deltaRWFake = wLepton->DeltaR(*fakeLepton);

  const double met = fEvent->pfMET;
  const double phiMET= fEvent->pfMETPhi;
  const double pxMET = met * cos(phiMET);
  const double pyMET = met * sin(phiMET);
  TLorentzVector lMET(pxMET, pyMET, 0., met);
  const double mt = sqrt(2 * met * wPt * (1 - cos(wLepton->DeltaPhi(lMET))));
  const double mass2L = (*wLepton + *fakeLepton).M();

  hWPt[3]->Fill(wPt);
  hFakePt[3]->Fill(fakePt);
  hDeltaRWFake[3]->Fill(deltaRWFake);

  hMET[3]->Fill(met);
  hMt[3]->Fill(mt);
  h2LMass[3]->Fill(mass2L);

  hWPt[fEvent->GetFinalState()]->Fill(wPt);
  hFakePt[fEvent->GetFinalState()]->Fill(fakePt);
  hDeltaRWFake[fEvent->GetFinalState()]->Fill(deltaRWFake);

  hMET[fEvent->GetFinalState()]->Fill(met);
  hMt[fEvent->GetFinalState()]->Fill(mt);
  h2LMass[fEvent->GetFinalState()]->Fill(mass2L);
}


void WJetsSelection::Finish()
{
  cout << endl << "Done." << endl;

  cout << "Analyzed Events = " << GetNAnalyzed() << "\n"
       << "Selected events = " << GetNSelected() << "\n\n";

  cout << "CHANNEL \tPreselection \tSSOF Selection \tFull Selection" << "\n";
  for (int i = 0; i <= 3; i++) {
    cout << i << "\t" << yieldsByChannelPreselection[i]
              << "\t" << yieldsByChannelSSOFSelection[i]
              << "\t" << yieldsByChannelFullSelection[i] << "\n";
  }
  cout << endl;
}

