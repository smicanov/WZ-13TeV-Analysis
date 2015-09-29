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

  for (unsigned int i = 0; i < 6; i++) {
    hWPt[i] = bookTH1D(("hWPt_" + boost::lexical_cast<string>(i)).c_str(),
                       "W lepton p_{t}", 100, 0, 200);
    hWEta[i] = bookTH1D(("hWEta_" + boost::lexical_cast<string>(i)).c_str(),
                        "W lepton #eta", 50, -2.5, 2.5);
    hWPhi[i] = bookTH1D(("hWPhi_" + boost::lexical_cast<string>(i)).c_str(),
                        "W lepton #phi", 72, -3.1416, 3.1416);
    hWRelIso[i] = bookTH1D(("hWRelIso_" + boost::lexical_cast<string>(i)).c_str(),
                           "W lepton relIso", 100, 0, 0.2);

    hFakePt[i] = bookTH1D(("hFakePt_" + boost::lexical_cast<string>(i)).c_str(),
                          "Fake lepton p_{t}", 100, 0, 100);
    hFakeEta[i] = bookTH1D(("hFakeEta_" + boost::lexical_cast<string>(i)).c_str(),
                           "Fake lepton #eta", 50, -2.5, 2.5);
    hFakePhi[i] = bookTH1D(("hFakePhi_" + boost::lexical_cast<string>(i)).c_str(),
                           "Fake lepton #phi", 72, -3.1416, 3.1416);
    hFakeRelIso[i] = bookTH1D(("hFakeRelIso_" + boost::lexical_cast<string>(i)).c_str(),
                              "Fake lepton relIso", 100, 0, 3);


    hDeltaPhiWFake[i] = bookTH1D(("hDeltaPhiWFake_" + boost::lexical_cast<string>(i)).c_str(),
                                 "|#Delta#Phi (W, Fake)|", 72, 0, 3.1416);
    hDeltaRWFake[i] = bookTH1D(("hDeltaRWFake_" + boost::lexical_cast<string>(i)).c_str(),
                               "#DeltaR (W, Fake)", 100, 0, 5);

    hMET[i] = bookTH1D(("hMET_" + boost::lexical_cast<string>(i)).c_str(),
                       "Missing E_{t}", 100, 0, 200);
    hMETPhi[i] = bookTH1D(("hMETPhi_" + boost::lexical_cast<string>(i)).c_str(),
                          "Missing E_{t} #phi", 100, -3.1416, 3.1416);

    hDeltaPhiWMET[i] = bookTH1D(("hDeltaPhiWMET_" + boost::lexical_cast<string>(i)).c_str(),
                                 "|#Delta#Phi (W, miss E_{t})|", 72, 0, 3.1416);
    hDeltaRWMET[i] = bookTH1D(("hDeltaRWMET_" + boost::lexical_cast<string>(i)).c_str(),
                               "#DeltaR (W, miss E_{t})", 100, 0, 5);

    hDeltaPhiFakeMET[i] = bookTH1D(("hDeltaPhiFakeMET_" + boost::lexical_cast<string>(i)).c_str(),
                                 "|#Delta#Phi (Fake, miss E_{t})|", 72, 0, 3.1416);
    hDeltaRFakeMET[i] = bookTH1D(("hDeltaRFakeMET_" + boost::lexical_cast<string>(i)).c_str(),
                               "#DeltaR (Fake, miss E_{t})", 100, 0, 5);

    hMt[i] = bookTH1D(("hMt_" + boost::lexical_cast<string>(i)).c_str(),
                       "M_{t}", 100, 0, 200);
    h2LMass[i] = bookTH1D(("h2LMass_" + boost::lexical_cast<string>(i)).c_str(),
                          "mass(2l)", 100, 0, 200);
  }

  for (int i = 0; i < 6; i++) {
    yieldsByChannelPreselection[i] = 0;
    yieldsByChannelSSSelection[i] = 0;
    yieldsByChannelOFSelection[i] = 0;
    yieldsByChannelFullSelection[i] = 0;
  }

/*
// Setup selected event lists 
  for (int i = 1; i <= 4; i++) {
    ostringstream outputFileName1;
    outputFileName1 << "/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/fom/test_WEleTight_FOMLoose-NoRelIso/WJets_Preselection_" << i << ".txt";
    cout << "File name : " << outputFileName1.str() << endl;
    eventLists1[i-1].open(outputFileName1.str().c_str());

    ostringstream outputFileName2;
    outputFileName2 << "/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/fom/test_WEleTight_FOMLoose-NoRelIso/WJets_Preselection_" << i << ".txt";
    cout << "File name : " << outputFileName2.str() << endl;
    eventLists2[i-1].open(outputFileName2.str().c_str());

    ostringstream outputFileName3;
    outputFileName3 << "/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/fom/test_WEleTight_FOMLoose-NoRelIso/WJets_SSOFSelection_" << i << ".txt";
    cout << "File name : " << outputFileName3.str() << endl;
    eventLists3[i-1].open(outputFileName3.str().c_str());

    ostringstream outputFileName4;
    outputFileName4 << "/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/fom/test_WEleTight_FOMLoose-NoRelIso/WJets_FullSelection_" << i << ".txt";
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
    yieldsByChannelFullSelection[5]++;
//    fEvent->DumpEvent(eventLists4[fEvent->GetFinalState()-1], 1);
  }

  if (fEvent->PassesOFSelection()) {
    yieldsByChannelOFSelection[fEvent->GetFinalState()]++;
    yieldsByChannelOFSelection[5]++;
//    fEvent->DumpEvent(eventLists3[fEvent->GetFinalState()-1], 1);
  }

  if (fEvent->PassesSSSelection()) {
    yieldsByChannelSSSelection[fEvent->GetFinalState()]++;
    yieldsByChannelSSSelection[5]++;
//    fEvent->DumpEvent(eventLists2[fEvent->GetFinalState()-1], 1);
  }

  if (fEvent->PassesPreselection()) {
    yieldsByChannelPreselection[fEvent->GetFinalState()]++;
    yieldsByChannelPreselection[5]++;
//    fEvent->DumpEvent(eventLists1[fEvent->GetFinalState()-1], 1);
  }

  if (!(fEvent->PassesPreselection()))  return;
//  if (!(fEvent->PassesSSSelection()))  return;
//  if (!(fEvent->PassesFullSelection()))  return;

  nSelectedEvents++;

  Lepton* wLepton = fEvent->fLeptons.at(fEvent->fCandidateLeptonIndex.first);
  const double wPt = wLepton->Pt();
  const double wEta = wLepton->Eta();
  const double wPhi = wLepton->Phi();
  const double wRelIso = wLepton->GetRelIso();

  Lepton* fakeLepton = fEvent->fLeptons.at(fEvent->fCandidateLeptonIndex.second);
  const double fakePt = fakeLepton->Pt();
  const double fakeEta = fakeLepton->Eta();
  const double fakePhi = fakeLepton->Phi();
  const double fakeRelIso = fakeLepton->GetRelIso();

  const double deltaPhiWFake = abs(wLepton->DeltaPhi(*fakeLepton));
  const double deltaRWFake = wLepton->DeltaR(*fakeLepton);

  const double met = fEvent->pfMET;
  const double phiMET= fEvent->pfMETPhi;
  const double pxMET = met * cos(phiMET);
  const double pyMET = met * sin(phiMET);
  TLorentzVector lMET(pxMET, pyMET, 0., met);

  const double deltaPhiWMET = abs(wLepton->DeltaPhi(lMET));
  const double deltaRWMET = wLepton->DeltaR(lMET);

  const double deltaPhiFakeMET = abs(lMET.DeltaPhi(*fakeLepton));
  const double deltaRFakeMET = lMET.DeltaR(*fakeLepton);

  const double mt = sqrt(2 * met * wPt * (1 - cos(wLepton->DeltaPhi(lMET))));
  const double mass2L = (*wLepton + *fakeLepton).M();

  hWPt[5]->Fill(wPt);
  hWEta[5]->Fill(wEta);
  hWPhi[5]->Fill(wPhi);
  hWRelIso[5]->Fill(wRelIso);

  hFakePt[5]->Fill(fakePt);
  hFakeEta[5]->Fill(fakeEta);
  hFakePhi[5]->Fill(fakePhi);
  hFakeRelIso[5]->Fill(fakeRelIso);

  hDeltaPhiWFake[5]->Fill(deltaPhiWFake);
  hDeltaRWFake[5]->Fill(deltaRWFake);

  hMET[5]->Fill(met);
  hMETPhi[5]->Fill(phiMET);

  hDeltaPhiWMET[5]->Fill(deltaPhiWMET);
  hDeltaRWMET[5]->Fill(deltaRWMET);

  hDeltaPhiFakeMET[5]->Fill(deltaPhiFakeMET);
  hDeltaRFakeMET[5]->Fill(deltaRFakeMET);

  hMt[5]->Fill(mt);
  h2LMass[5]->Fill(mass2L);


  hWPt[fEvent->GetFinalState()]->Fill(wPt);
  hWEta[fEvent->GetFinalState()]->Fill(wEta);
  hWPhi[fEvent->GetFinalState()]->Fill(wPhi);
  hWRelIso[fEvent->GetFinalState()]->Fill(wRelIso);

  hFakePt[fEvent->GetFinalState()]->Fill(fakePt);
  hFakeEta[fEvent->GetFinalState()]->Fill(fakeEta);
  hFakePhi[fEvent->GetFinalState()]->Fill(fakePhi);
  hFakeRelIso[fEvent->GetFinalState()]->Fill(fakeRelIso);

  hDeltaPhiWFake[fEvent->GetFinalState()]->Fill(deltaPhiWFake);
  hDeltaRWFake[fEvent->GetFinalState()]->Fill(deltaRWFake);

  hMET[fEvent->GetFinalState()]->Fill(met);
  hMETPhi[fEvent->GetFinalState()]->Fill(phiMET);

  hDeltaPhiWMET[fEvent->GetFinalState()]->Fill(deltaPhiWMET);
  hDeltaRWMET[fEvent->GetFinalState()]->Fill(deltaRWMET);

  hDeltaPhiFakeMET[fEvent->GetFinalState()]->Fill(deltaPhiFakeMET);
  hDeltaRFakeMET[fEvent->GetFinalState()]->Fill(deltaRFakeMET);

  hMt[fEvent->GetFinalState()]->Fill(mt);
  h2LMass[fEvent->GetFinalState()]->Fill(mass2L);
}


void WJetsSelection::Finish()
{
  cout << endl << "Done." << endl;

  cout << "Analyzed events : " << GetNAnalyzed() << "\n"
       << "Selected events : " << GetNSelected() << "\n\n";

  cout << "CHANNEL\tPreselection\tSS Selection\tOF Selection\tFull Selection" << "\n";
  for (int i = 0; i < 6; i++) {
    cout << i << "\t" << yieldsByChannelPreselection[i]
              << "\t" << yieldsByChannelSSSelection[i]
              << "\t" << yieldsByChannelOFSelection[i]
              << "\t" << yieldsByChannelFullSelection[i] << "\n";
  }
  cout << endl;
}

