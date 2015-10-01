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
                              "Fake lepton relIso", 80, 0, 4);


    hDeltaPhiWFake[i] = bookTH1D(("hDeltaPhiWFake_" + boost::lexical_cast<string>(i)).c_str(),
                                 "|#Delta#Phi (W l, Fake l)|", 72, 0, 3.1416);
    hDeltaRWFake[i] = bookTH1D(("hDeltaRWFake_" + boost::lexical_cast<string>(i)).c_str(),
                               "#DeltaR (W l, Fake l)", 100, 0, 5);

    hMET[i] = bookTH1D(("hMET_" + boost::lexical_cast<string>(i)).c_str(),
                       "Missing E_{t}", 100, 0, 200);
    hMETPhi[i] = bookTH1D(("hMETPhi_" + boost::lexical_cast<string>(i)).c_str(),
                          "Missing E_{t} #phi", 100, -3.1416, 3.1416);

    hDeltaPhiWMET[i] = bookTH1D(("hDeltaPhiWMET_" + boost::lexical_cast<string>(i)).c_str(),
                                 "|#Delta#Phi (W l, miss E_{t})|", 72, 0, 3.1416);
    hDeltaRWMET[i] = bookTH1D(("hDeltaRWMET_" + boost::lexical_cast<string>(i)).c_str(),
                               "#DeltaR (W l, miss E_{t})", 100, 0, 5);

    hDeltaPhiFakeMET[i] = bookTH1D(("hDeltaPhiFakeMET_" + boost::lexical_cast<string>(i)).c_str(),
                                 "|#Delta#Phi (Fake l, miss E_{t})|", 72, 0, 3.1416);
    hDeltaRFakeMET[i] = bookTH1D(("hDeltaRFakeMET_" + boost::lexical_cast<string>(i)).c_str(),
                               "#DeltaR (Fake l, miss E_{t})", 100, 0, 5);

    hWMt[i] = bookTH1D(("hWMt_" + boost::lexical_cast<string>(i)).c_str(),
                       "W M_{t}", 100, 0, 200);
    hFakeMt[i] = bookTH1D(("hFakeMt_" + boost::lexical_cast<string>(i)).c_str(),
                       "Fake M_{t}", 100, 0, 200);

    hDeltaPhiWBosonWl[i] = bookTH1D(("hDeltaPhiWBosonWl_" + boost::lexical_cast<string>(i)).c_str(),
                                 "|#Delta#Phi (W boson, W lepton)|", 72, 0, 3.1416);
    hDeltaRWBosonWl[i] = bookTH1D(("hDeltaRWBosonWl_" + boost::lexical_cast<string>(i)).c_str(),
                               "#DeltaR (W boson, W lepton)", 100, 0, 5);
    hDeltaPhiWBosonFakel[i] = bookTH1D(("hDeltaPhiWBosonFakel_" + boost::lexical_cast<string>(i)).c_str(),
                                 "|#Delta#Phi (W boson, Fake lepton)|", 72, 0, 3.1416);
    hDeltaRWBosonFakel[i] = bookTH1D(("hDeltaRWBosonFakel_" + boost::lexical_cast<string>(i)).c_str(),
                               "#DeltaR (W boson, Fake lepton)", 100, 0, 5);

    hDeltaPhiFakeBosonWl[i] = bookTH1D(("hDeltaPhiFakeBosonWl_" + boost::lexical_cast<string>(i)).c_str(),
                                 "|#Delta#Phi (Fake W boson, W lepton)|", 72, 0, 3.1416);
    hDeltaRFakeBosonWl[i] = bookTH1D(("hDeltaRFakeBosonWl_" + boost::lexical_cast<string>(i)).c_str(),
                               "#DeltaR (Fake W boson, W lepton)", 100, 0, 5);
    hDeltaPhiFakeBosonFakel[i] = bookTH1D(("hDeltaPhiFakeBosonFakel_" + boost::lexical_cast<string>(i)).c_str(),
                                 "|#Delta#Phi (Fake W boson, Fake lepton)|", 72, 0, 3.1416);
    hDeltaRFakeBosonFakel[i] = bookTH1D(("hDeltaRFakeBosonFakel_" + boost::lexical_cast<string>(i)).c_str(),
                               "#DeltaR (Fake W boson, Fake lepton)", 100, 0, 5);

    h2LMass[i] = bookTH1D(("h2LMass_" + boost::lexical_cast<string>(i)).c_str(),
                          "mass(2l)", 100, 0, 200);

    hGoodJets15[i] = bookTH1D(("hGoodJets15_" + boost::lexical_cast<string>(i)).c_str(),
                            "Number of good Jets (p_{t}^{Good} > 15 GeV)", 11, -0.5, 10.5);
    hDRminGoodJet15Wl[i] = bookTH1D(("hDRminGoodJet15Wl_" + boost::lexical_cast<string>(i)).c_str(),
                              "#DeltaR_{min} (W l, good jet (15 GeV))", 70, 0, 7);
    hDRminGoodJet15Fakel[i] = bookTH1D(("hDRminGoodJet15Fakel_" + boost::lexical_cast<string>(i)).c_str(),
                                 "#DeltaR_{min} (Fake l, good jet (15 GeV))", 70, 0, 7);

    hGoodJetsLead20[i] = bookTH1D(("hGoodJetsLead20_" + boost::lexical_cast<string>(i)).c_str(),
                            "Number of good Jets (p_{t}^{Lead} > 20 GeV)", 11, -0.5, 10.5);
    hDRminGoodJetWlLead20[i] = bookTH1D(("hDRminGoodJetWlLead20_" + boost::lexical_cast<string>(i)).c_str(),
                              "#DeltaR_{min} (W l, good jet (20 GeV))", 70, 0, 7);
    hDRminGoodJetFakelLead20[i] = bookTH1D(("hDRminGoodJetFakelLead20_" + boost::lexical_cast<string>(i)).c_str(),
                                 "#DeltaR_{min} (Fake l, good jet (20 Gev))", 70, 0, 7);

    hLeadJetPt[i] = bookTH1D(("hLeadJetPt_" + boost::lexical_cast<string>(i)).c_str(),
                          "Lead jet p_{t}", 100, 0, 200);
    hLeadJetEta[i] = bookTH1D(("hLeadJetEta_" + boost::lexical_cast<string>(i)).c_str(),
                           "Lead jet #eta", 50, -2.5, 2.5);
    hLeadJetPhi[i] = bookTH1D(("hLeadJetPhi_" + boost::lexical_cast<string>(i)).c_str(),
                           "Lead jet #phi", 72, -3.1416, 3.1416);
    hLeadJetEt[i] = bookTH1D(("hLeadJetEt_" + boost::lexical_cast<string>(i)).c_str(),
                          "Lead jet E_{t}", 100, 0, 200);

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

  if (fEvent->PassesSSSelection() && fEvent->PassesTight()) {
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

//  if (!(fEvent->PassesPreselection()))  return;
  if (!(fEvent->PassesSSSelection()))  return;
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
  const TLorentzVector lMET(pxMET, pyMET, 0., met);

  const double deltaPhiWMET = abs(wLepton->DeltaPhi(lMET));
  const double deltaRWMET = wLepton->DeltaR(lMET);

  const double deltaPhiFakeMET = abs(lMET.DeltaPhi(*fakeLepton));
  const double deltaRFakeMET = lMET.DeltaR(*fakeLepton);

  const double wMt = sqrt(2 * met * wPt * (1 - cos(wLepton->DeltaPhi(lMET))));
  const double fakeMt = sqrt(2 * met * fakePt * (1 - cos(fakeLepton->DeltaPhi(lMET))));

  const TLorentzVector wBoson = lMET + *(wLepton);
  const double deltaPhiWBosonWl = abs(wLepton->DeltaPhi(wBoson));
  const double deltaRWBosonWl = wLepton->DeltaR(wBoson);
  const double deltaPhiWBosonFakel = abs(fakeLepton->DeltaPhi(wBoson));
  const double deltaRWBosonFakel = fakeLepton->DeltaR(wBoson);

  const TLorentzVector fakeBoson = lMET + *(fakeLepton);
  const double deltaPhiFakeBosonWl = abs(fakeBoson.DeltaPhi(*wLepton));
  const double deltaRFakeBosonWl = fakeBoson.DeltaR(*wLepton);
  const double deltaPhiFakeBosonFakel = abs(fakeBoson.DeltaPhi(*fakeLepton));
  const double deltaRFakeBosonFakel = fakeBoson.DeltaR(*fakeLepton);

  const double mass2L = (*wLepton + *fakeLepton).M();

  unsigned int nGoodJets = 0;
  vector<double> dRGoodJetWl;
  vector<double> dRGoodJetFakel;
  for (vector<unsigned int>::const_iterator jIt = fEvent->fGoodJetsIndex.begin();
       jIt != fEvent->fGoodJetsIndex.end(); ++jIt) {
    if (fEvent->fJets.at(*jIt)->Pt() > GOODJET_PTCUT) {
      nGoodJets++;
      const double dRjWl = fEvent->fJets.at(*jIt)->DeltaR(*wLepton);
      dRGoodJetWl.push_back(dRjWl);
      const double dRjFakel = fEvent->fJets.at(*jIt)->DeltaR(*fakeLepton);
      dRGoodJetFakel.push_back(dRjFakel);
    }
  }

  hGoodJets15[5]->Fill(nGoodJets);
  hGoodJets15[fEvent->GetFinalState()]->Fill(nGoodJets);

  if (!(dRGoodJetWl.empty())) {
    hDRminGoodJet15Wl[5]->Fill(*(min_element(dRGoodJetWl.begin(), dRGoodJetWl.end())));
    hDRminGoodJet15Wl[fEvent->GetFinalState()]->Fill(*(min_element(dRGoodJetWl.begin(), dRGoodJetWl.end())));
    hDRminGoodJet15Fakel[5]->Fill(*(min_element(dRGoodJetFakel.begin(), dRGoodJetFakel.end())));
    hDRminGoodJet15Fakel[fEvent->GetFinalState()]->Fill(*(min_element(dRGoodJetFakel.begin(), dRGoodJetFakel.end())));
  }

  if (!(fEvent->fGoodJetsIndex.empty())) {
    vector<Jet> jets;
    for (vector<unsigned int>::const_iterator jIt = fEvent->fGoodJetsIndex.begin();
         jIt != fEvent->fGoodJetsIndex.end(); ++jIt) {
      const Jet goodJet = *(fEvent->fJets.at(*jIt));
      jets.push_back(goodJet);
    }
    sort(jets.begin(), jets.end(), HigherPtJets);

    const double leadJetPt = jets.front().Pt();
    const double leadJetEta = jets.front().Eta();
    const double leadJetPhi = jets.front().Phi();
    const double leadJetEt = jets.front().Et();

    hLeadJetPt[5]->Fill(leadJetPt);
    hLeadJetPt[fEvent->GetFinalState()]->Fill(leadJetPt);
    hLeadJetEta[5]->Fill(leadJetEta);
    hLeadJetEta[fEvent->GetFinalState()]->Fill(leadJetEta);
    hLeadJetPhi[5]->Fill(leadJetPhi);
    hLeadJetPhi[fEvent->GetFinalState()]->Fill(leadJetPhi);
    hLeadJetEt[5]->Fill(leadJetEt);
    hLeadJetEt[fEvent->GetFinalState()]->Fill(leadJetEt);

    if (leadJetPt > LEADJET_PTCUT) {
      hGoodJetsLead20[5]->Fill(fEvent->fGoodJetsIndex.size());
      hGoodJetsLead20[fEvent->GetFinalState()]->Fill(fEvent->fGoodJetsIndex.size());

      vector<double> dRJetWl;
      vector<double> dRJetFakel;
      for (vector<Jet*>::const_iterator jIt = fEvent->fJets.begin(); jIt != fEvent->fJets.end(); ++jIt) {
        const double dRjWl = (*jIt)->DeltaR(*wLepton);
        dRJetWl.push_back(dRjWl);
        const double dRjFakel = (*jIt)->DeltaR(*fakeLepton);
        dRJetFakel.push_back(dRjFakel);
      }

      hDRminGoodJetWlLead20[5]->Fill(*(min_element(dRJetWl.begin(), dRJetWl.end())));
      hDRminGoodJetWlLead20[fEvent->GetFinalState()]->Fill(*(min_element(dRJetWl.begin(), dRJetWl.end())));
      hDRminGoodJetFakelLead20[5]->Fill(*(min_element(dRJetFakel.begin(), dRJetFakel.end())));
      hDRminGoodJetFakelLead20[fEvent->GetFinalState()]->Fill(*(min_element(dRJetFakel.begin(), dRJetFakel.end())));
    } else {
      hGoodJetsLead20[5]->Fill(0.);
      hGoodJetsLead20[fEvent->GetFinalState()]->Fill(0.);
    }
  } else {
    hGoodJetsLead20[5]->Fill(0.);
    hGoodJetsLead20[fEvent->GetFinalState()]->Fill(0.);
  }

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

  hWMt[5]->Fill(wMt);
  hFakeMt[5]->Fill(fakeMt);

  hDeltaPhiWBosonWl[5]->Fill(deltaPhiWBosonWl);
  hDeltaRWBosonWl[5]->Fill(deltaRWBosonWl);
  hDeltaPhiWBosonFakel[5]->Fill(deltaPhiWBosonFakel);
  hDeltaRWBosonFakel[5]->Fill(deltaRWBosonFakel);

  hDeltaPhiFakeBosonWl[5]->Fill(deltaPhiFakeBosonWl);
  hDeltaRFakeBosonWl[5]->Fill(deltaRFakeBosonWl);
  hDeltaPhiFakeBosonFakel[5]->Fill(deltaPhiFakeBosonFakel);
  hDeltaRFakeBosonFakel[5]->Fill(deltaRFakeBosonFakel);

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

  hWMt[fEvent->GetFinalState()]->Fill(wMt);
  hFakeMt[fEvent->GetFinalState()]->Fill(fakeMt);

  hDeltaPhiWBosonWl[fEvent->GetFinalState()]->Fill(deltaPhiWBosonWl);
  hDeltaRWBosonWl[fEvent->GetFinalState()]->Fill(deltaRWBosonWl);
  hDeltaPhiWBosonFakel[fEvent->GetFinalState()]->Fill(deltaPhiWBosonFakel);
  hDeltaRWBosonFakel[fEvent->GetFinalState()]->Fill(deltaRWBosonFakel);

  hDeltaPhiFakeBosonWl[fEvent->GetFinalState()]->Fill(deltaPhiFakeBosonWl);
  hDeltaRFakeBosonWl[fEvent->GetFinalState()]->Fill(deltaRFakeBosonWl);
  hDeltaPhiFakeBosonFakel[fEvent->GetFinalState()]->Fill(deltaPhiFakeBosonFakel);
  hDeltaRFakeBosonFakel[fEvent->GetFinalState()]->Fill(deltaRFakeBosonFakel);

  h2LMass[fEvent->GetFinalState()]->Fill(mass2L);
}


void WJetsSelection::Finish()
{
  cout << endl << "Done." << endl;

  cout << "Analyzed events : " << GetNAnalyzed() << "\n"
       << "Selected events : " << GetNSelected() << "\n\n";

  cout << "CHANNEL\tPreselection\tSS Selection\tTight Selection\tFull Selection" << "\n";
  for (int i = 0; i < 6; i++) {
    cout << i << "\t" << yieldsByChannelPreselection[i]
              << "\t" << yieldsByChannelSSSelection[i]
              << "\t" << yieldsByChannelOFSelection[i]
              << "\t" << yieldsByChannelFullSelection[i] << "\n";
  }
  cout << endl;
}

