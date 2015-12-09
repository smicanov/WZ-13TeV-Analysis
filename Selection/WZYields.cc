#include "WZYields.h"
#include "Constants.h"
#include "TLatex.h"

#include <cmath>
#include <sstream>
#include <string>
#include <boost/lexical_cast.hpp>


using namespace std;


WZYields::WZYields(WZEvent* e, TFile* outputFile)
 : Generic(e, outputFile)
{

}


void WZYields::Init()
{
  nAnalyzedEvents = 0;
  nSelectedEvents = 0;

  for (unsigned int i = 0; i <= 4; i++) {
    hZmass[i] = bookTH1D(("hZmass_" + boost::lexical_cast<string>(i)).c_str(),
                         "Z mass", 64, 58, 122);
    hZpt[i] = bookTH1D(("hZpt_" + boost::lexical_cast<string>(i)).c_str(),
                       "Z pt", 100, 0, 200);
    hMET[i] = bookTH1D(("hMET_" + boost::lexical_cast<string>(i)).c_str(),
                       "Missing Et", 100, 0, 200);
    hMt[i] = bookTH1D(("hMt_" + boost::lexical_cast<string>(i)).c_str(),
                       "Mt", 100, 0, 200);
    hZl1pt[i] = bookTH1D(("hZl1pt_" + boost::lexical_cast<string>(i)).c_str(),
                         "Z Lepton1 pt", 100, 0, 200);
    hZl2pt[i] = bookTH1D(("hZl2pt_" + boost::lexical_cast<string>(i)).c_str(),
                         "Z Lepton2 pt", 100, 0, 200);
    hWlpt[i] = bookTH1D(("hWlpt_" + boost::lexical_cast<string>(i)).c_str(),
                        "W Lepton pt", 100, 0, 200);
    hNJets[i] = bookTH1D(("hNJets_" + boost::lexical_cast<string>(i)).c_str(),
                         "No. Jets", 7, -0.5, 6.5);
    hNJetsNoMuIso[i] = bookTH1D(("hNJetsNoMuIso_" + boost::lexical_cast<string>(i)).c_str(),
                                "No. Jets", 11, -0.5, 10.5);
    hNJetsNoEleIso[i] = bookTH1D(("hNJetsNoEleIso_" + boost::lexical_cast<string>(i)).c_str(),
                                 "No. Jets", 11, -0.5, 10.5);
    hNJetsNoIso[i] = bookTH1D(("hNJetsNoIso_" + boost::lexical_cast<string>(i)).c_str(),
                              "No. Jets", 11, -0.5, 10.5);
    h3LMass[i] = bookTH1D(("h3LMass_" + boost::lexical_cast<string>(i)).c_str(),
                          "3L Mass", 150, 50, 350);
    hDeltaR[i] = bookTH1D(("hDeltaR_" + boost::lexical_cast<string>(i)).c_str(),
                          "#deltaR (Wl, Zl)", 100, 0, 5);
    hDeltaRMin[i] = bookTH1D(("hDeltaRMin_" + boost::lexical_cast<string>(i)).c_str(),
                             "#deltaR_{min} (Wl, Zl)", 100, 0, 5);
  }

  for (int i = 0; i <= 5; i++) {
    yieldsByChannelPreselection[i] = 0;
    yieldsByChannelZSelection[i] = 0;
    yieldsByChannelWSelection[i] = 0;
    yieldsByChannelFullSelection[i] = 0;
    yieldsByChannelFullSelectionProblematic[i] = 0;
  }

/*
  // Setup selected event lists 
  for (int i = 1; i <= 4; i++) {

    ostringstream outputFileName1;
    outputFileName1 << "/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/yields/synchronization/data/WZ_FullSelection_" << i << "_test.txt";
    cout << "File name : " << outputFileName1.str() << endl;
    eventLists1[i-1].open(outputFileName1.str().c_str());

    ostringstream outputFileName2;
    outputFileName2 << "/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/yields/synchronization/data/WZ_WSelection_" << i << "_test.txt";
    cout << "File name : " << outputFileName2.str() << endl;
    eventLists2[i-1].open(outputFileName2.str().c_str());

    ostringstream outputFileName3;
    outputFileName3 << "/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/yields/synchronization/data/WZ_ZSelection_" << i << "_test.txt";
    cout << "File name : " << outputFileName3.str() << endl;
    eventLists3[i-1].open(outputFileName3.str().c_str());

    ostringstream outputFileName4;
    outputFileName4 << "/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/yields/synchronization/data/WZ_Preselection_" << i << "_test.txt";
    cout << "File name : " << outputFileName4.str() << endl;
    eventLists4[i-1].open(outputFileName4.str().c_str());
  }
*/

  minRun = maxRun = 0;
  runNumber.clear();
}


void WZYields::Analysis()
{
  nAnalyzedEvents++;

  if (fWZEvent->PassesTriggerMuonEG() && fWZEvent->PassesFullSelection() && !(fWZEvent->PassesMETFilters()))
    cout << "Selected event FAILED MET Filter !!!\t-\t" << fWZEvent->run << ":" << fWZEvent->lumis << ":" << fWZEvent->event << endl;

  bool test = false;
  for (vector<Lepton*>::const_iterator lIt = fWZEvent->fLeptons.begin();
       lIt != fWZEvent->fLeptons.end(); ++lIt)
    if ((*lIt)->IsTight() && !((*lIt)->IsLoose()))  test = true;
  if (test)  nSelectedEvents++;

  runNumber.push_back(fWZEvent->run);

  if (fWZEvent->PassesTriggerMuonEG() && fWZEvent->PassesFullSelection() && fWZEvent->PassesMETFilters()) {
    yieldsByChannelFullSelection[fWZEvent->GetFinalState()]++;
    yieldsByChannelFullSelection[5]++;
//    fWZEvent->Dump(eventLists1[fWZEvent->GetFinalState()-1], 10);
    if (!((fWZEvent->GetZLeptons().first)->IsLoose()) ||
        !((fWZEvent->GetZLeptons().second)->IsLoose()) ||
        !((fWZEvent->GetWLepton())->IsLoose())) {
      yieldsByChannelFullSelectionProblematic[fWZEvent->GetFinalState()]++;
      yieldsByChannelFullSelectionProblematic[5]++;
      cout << "Event info : " << fWZEvent->run << ":" << fWZEvent->lumis << ":" << fWZEvent->event << endl;
    }
  }

  if (fWZEvent->PassesTriggerMuonEG() && fWZEvent->PassesWSelection() && fWZEvent->PassesMETFilters()) {
    yieldsByChannelWSelection[fWZEvent->GetFinalState()]++;
    yieldsByChannelWSelection[5]++;
//    fWZEvent->Dump(eventLists2[fWZEvent->GetFinalState()-1], 10);
  }

  if (fWZEvent->PassesTriggerMuonEG() && fWZEvent->PassesZSelection() && fWZEvent->PassesMETFilters()){
    yieldsByChannelZSelection[fWZEvent->GetFinalState()]++;
    yieldsByChannelZSelection[5]++;
//    fWZEvent->Dump(eventLists3[fWZEvent->GetFinalState()-1], 7);
  }

  if (fWZEvent->PassesTriggerMuonEG() && fWZEvent->PassesPreselection() && fWZEvent->PassesMETFilters()) {
    yieldsByChannelPreselection[fWZEvent->GetFinalState()]++;
    yieldsByChannelPreselection[5]++;
//    fWZEvent->Dump(eventLists4[fWZEvent->GetFinalState()-1], 5);
  }

  if (fWZEvent->PassesTriggerMuonEG() && fWZEvent->PassesMETFilters()) {
    yieldsByChannelTrigger[fWZEvent->GetFinalState()]++;
    yieldsByChannelTrigger[5]++;
  }

  if (!(fWZEvent->PassesTriggerMuonEG()) || !(fWZEvent->PassesFullSelection()) || !(fWZEvent->PassesMETFilters()))  return;

  const double massZ = (*(fWZEvent->GetZLeptons().first) + *(fWZEvent->GetZLeptons().second)).M();
  const double ptZ = (*(fWZEvent->GetZLeptons().first) + *(fWZEvent->GetZLeptons().second)).Pt();
  const double ptZl1 = fWZEvent->GetZLeptons().first->Pt();
  const double ptZl2 = fWZEvent->GetZLeptons().second->Pt();
  const double ptWl = fWZEvent->GetWLepton()->Pt();
  const double met = fWZEvent->pfMET;
  const double phiMET= fWZEvent->pfMETPhi;
  const double pxMET = met * cos(phiMET);
  const double pyMET = met * sin(phiMET);
  TLorentzVector lMET(pxMET, pyMET, 0., met);
  const double mt = sqrt(2 * met * ptWl * (1 - cos(fWZEvent->GetWLepton()->DeltaPhi(lMET))));
  const double mass3L = (*(fWZEvent->GetZLeptons().first) + *(fWZEvent->GetZLeptons().second) +
                         *(fWZEvent->GetWLepton())).M();
  const double dRWlZl1 = fWZEvent->GetWLepton()->DeltaR(*(fWZEvent->GetZLeptons().first));
  const double dRWlZl2 = fWZEvent->GetWLepton()->DeltaR(*(fWZEvent->GetZLeptons().second));
  const double minDR = min(dRWlZl1, dRWlZl2);

// Counting accompanying jets
  unsigned int nSelectedJets = 0;
  unsigned int nSelectedJetsNoMuIso = 0;
  unsigned int nSelectedJetsNoEleIso = 0;
  unsigned int nSelectedJetsNoIso = 0;
  for (unsigned int i = 0; i < fWZEvent->jetPt->size(); i++) {

//  in ggNtuplizer V07-04-05 use explicit cuts
/*
    if (!(fWZEvent->jetNHF->at(i) < 0.99) || !(fWZEvent->jetNEF->at(i) < 0.99) ||
        !(fWZEvent->jetCEF->at(i) < 0.99) || !(fWZEvent->jetNConstituents->at(i) > 1) ||
        !(fWZEvent->jetCHF->at(i) > 0) || !(fWZEvent->jetNCH->at(i) > 0))
      continue;
*/

//  while in version V07-04-09 use bool for 'PURE09' and 'LOOSE' (defined with cuts above)
    for (vector<bool>::const_iterator bIt = fWZEvent->jetPFLooseId->begin();
         bIt != fWZEvent->jetPFLooseId->end(); ++bIt)
      if (!(*bIt))  continue;

    const double ptJet = fWZEvent->jetPt->at(i);
    const double etaJet = fWZEvent->jetEta->at(i);
    if (!(ptJet > JET_PTCUT) || !(abs(etaJet) < JET_ETACUT))  continue;

    nSelectedJetsNoIso++;
    const double phiJet = fWZEvent->jetPhi->at(i);
    const double eJet = fWZEvent->jetEn->at(i);  // present only in V07-04-09+, otherwise use M=0
    TLorentzVector lJet;
    lJet.SetPtEtaPhiE(ptJet, etaJet, phiJet, eJet);
//    lJet.SetPtEtaPhiM(ptJet, etaJet, phiJet, 0);
    const double deltaRJetWl = fWZEvent->GetWLepton()->DeltaR(lJet);
    const double deltaRJetZl1 = fWZEvent->GetZLeptons().first->DeltaR(lJet);
    const double deltaRJetZl2 = fWZEvent->GetZLeptons().second->DeltaR(lJet);

    if (deltaRJetWl > LEPTONJET_DELTARMIN &&
        deltaRJetZl1 > LEPTONJET_DELTARMIN && deltaRJetZl2 > LEPTONJET_DELTARMIN)
      nSelectedJets++;

    if (fWZEvent->GetFinalState() == mmm) {
      nSelectedJetsNoMuIso++;
      if (deltaRJetWl > ELEJET_DELTARMIN &&
          deltaRJetZl1 > ELEJET_DELTARMIN && deltaRJetZl2 > ELEJET_DELTARMIN)
        nSelectedJetsNoEleIso++;
    }
    else if (fWZEvent->GetFinalState() == mme) {
      if (deltaRJetWl > ELEJET_DELTARMIN)  nSelectedJetsNoMuIso++;
      else if (deltaRJetZl1 > ELEJET_DELTARMIN && deltaRJetZl2 > ELEJET_DELTARMIN)
        nSelectedJetsNoEleIso++;
    }
    else if (fWZEvent->GetFinalState() == eem) {
      if (deltaRJetWl > ELEJET_DELTARMIN)  nSelectedJetsNoEleIso++;
      else if (deltaRJetZl1 > ELEJET_DELTARMIN && deltaRJetZl2 > ELEJET_DELTARMIN)
        nSelectedJetsNoMuIso++;
    }
    else if (fWZEvent->GetFinalState() == eee) {
      nSelectedJetsNoEleIso++;
      if (deltaRJetWl > ELEJET_DELTARMIN &&
          deltaRJetZl1 > ELEJET_DELTARMIN && deltaRJetZl2 > ELEJET_DELTARMIN)
        nSelectedJetsNoMuIso++;
    }
    else  continue;
  }

  hZmass[4]->Fill(massZ);
  hZpt[4]->Fill(ptZ);
  hMET[4]->Fill(met);
  hMt[4]->Fill(mt);
  hZl1pt[4]->Fill(ptZl1);
  hZl2pt[4]->Fill(ptZl2);
  hWlpt[4]->Fill(ptWl);
  hNJets[4]->Fill(nSelectedJets);
  hNJetsNoMuIso[4]->Fill(nSelectedJetsNoMuIso);
  hNJetsNoEleIso[4]->Fill(nSelectedJetsNoEleIso);
  hNJetsNoIso[4]->Fill(nSelectedJetsNoIso);
  h3LMass[4]->Fill(mass3L);
  hDeltaR[4]->Fill(dRWlZl1);
  hDeltaR[4]->Fill(dRWlZl2);
  hDeltaRMin[4]->Fill(minDR);

  hZmass[fWZEvent->GetFinalState()-1]->Fill(massZ);
  hZpt[fWZEvent->GetFinalState()-1]->Fill(ptZ);
  hMET[fWZEvent->GetFinalState()-1]->Fill(met);
  hMt[fWZEvent->GetFinalState()-1]->Fill(mt);
  hZl1pt[fWZEvent->GetFinalState()-1]->Fill(ptZl1);
  hZl2pt[fWZEvent->GetFinalState()-1]->Fill(ptZl2);
  hWlpt[fWZEvent->GetFinalState()-1]->Fill(ptWl);
  hNJets[fWZEvent->GetFinalState()-1]->Fill(nSelectedJets);
  hNJetsNoMuIso[fWZEvent->GetFinalState()-1]->Fill(nSelectedJetsNoMuIso);
  hNJetsNoEleIso[fWZEvent->GetFinalState()-1]->Fill(nSelectedJetsNoEleIso);
  hNJetsNoIso[fWZEvent->GetFinalState()-1]->Fill(nSelectedJetsNoIso);
  h3LMass[fWZEvent->GetFinalState()-1]->Fill(mass3L);
  hDeltaR[fWZEvent->GetFinalState()-1]->Fill(dRWlZl1);
  hDeltaR[fWZEvent->GetFinalState()-1]->Fill(dRWlZl2);
  hDeltaRMin[fWZEvent->GetFinalState()-1]->Fill(minDR);
}


void WZYields::Finish()
{
  cout << endl << "Done." << endl;

  cout << "Analyzed events : " << GetNAnalyzed() << "\n"
       << "Problematic events : " << GetNSelected() << "\n\n";

  vector<long unsigned int>::iterator min = min_element(runNumber.begin(), runNumber.end());
  minRun = *min;
  vector<long unsigned int>::iterator max = max_element(runNumber.begin(), runNumber.end());
  maxRun = *max;

  cout << "Minimum Run # : " << minRun << "\n"
       << "Maximum Run # : " << maxRun << "\n\n";

  cout << "CHANNEL\tTriggers\tPreselection\tZ Selection\tW Selection\tFull Selection\t\tProblematic"<< "\n";
  for (int i = 0; i <= 5; i++) {
    cout << i << "\t" << yieldsByChannelTrigger[i]
              << "\t" << yieldsByChannelPreselection[i]
              << "\t" << yieldsByChannelZSelection[i]
              << "\t" << yieldsByChannelWSelection[i]
              << "\t" << yieldsByChannelFullSelection[i]
              << "\t\t" << yieldsByChannelFullSelectionProblematic[i] << "\n";
  }

  cout << endl;
}

