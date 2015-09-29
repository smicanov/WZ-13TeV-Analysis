#include "Event.h"
#include "Constants.h"

#include <iomanip>


using namespace std;


Event::Event(TTree* tree) :
  BASECLASS(tree)
{
// Set static pointer to Particle class
  Particle::SetTree(this);
}


void
Event::Clear()
{

  fFinalState = undefined;
  fSelectionLevel = Undefined;

// Empty leptons
  for (vector<Lepton*>::iterator lIt = fLeptons.begin(); lIt != fLeptons.end(); ) {
    delete *lIt;  
    lIt = fLeptons.erase(lIt);
  }

  fCandidateLeptonIndex = make_pair(-1, -1);
}


void
Event::Read()
{
  Clear();

// Electrons
  for (unsigned int indexEle = 0; indexEle < eleCharge->size(); indexEle++) {
    Electron* ele = new Electron(indexEle, elePt->at(indexEle), eleEta->at(indexEle),
                                 elePhi->at(indexEle), eleCharge->at(indexEle));
    fLeptons.push_back(ele);
  }

// Muons
  for (unsigned int indexMu = 0; indexMu < muCharge->size(); indexMu++) {
    Muon* mu = new Muon(indexMu, muPt->at(indexMu), muEta->at(indexMu),
                        muPhi->at(indexMu), muCharge->at(indexMu));
    fLeptons.push_back(mu);
  }
}


bool
Event::PassesPreselection()
{
  if      (!(fSelectionLevel < Preselection))     return true;
  else if (fSelectionLevel == FailsPreselection)  return false;

  bool passed = false;
  fSelectionLevel = FailsPreselection;

  if((nEle + nMu) < N_LEPTONS)  return passed;

  unsigned int nFake = 0;
  vector<unsigned int> indexFake;
  for (vector<Lepton*>::iterator lIt = fLeptons.begin(); lIt != fLeptons.end(); ++lIt) {
    if ((*lIt)->PassesPtCut() && (*lIt)->PassesEtaCut() && (*lIt)->IsFOMLoose()) {
      nFake++;
      unsigned int index = distance(fLeptons.begin(), lIt);
      indexFake.push_back(index);
    }
  }

  if (nFake != indexFake.size()) {
    cout << "Error: Number of Fake candidates DIFFERS from number of Fake candidate indices !!!" << endl;
    cout << nFake << " != " << indexFake.size() << endl;
    return passed;
  }

// ### put the conditions on number of selected leptons below!!! ###
  if (indexFake.size() != N_LEPTONS)  return passed;
  else {
    if (fLeptons.at(indexFake.at(0))->IsTight() || fLeptons.at(indexFake.at(1))->IsTight())
      passed = true;
  }

  if (passed) {
    fSelectionLevel = Preselection;
    const unsigned int ind1 = indexFake.at(0);
    const unsigned int ind2 = indexFake.at(1);

    if (fLeptons.at(ind1)->IsTight() && !(fLeptons.at(ind2)->IsTight())) {
      fCandidateLeptonIndex.first = ind1;
      fCandidateLeptonIndex.second = ind2;
    } else if (!(fLeptons.at(ind1)->IsTight()) && fLeptons.at(ind2)->IsTight()) {
      fCandidateLeptonIndex.first = ind2;
      fCandidateLeptonIndex.second = ind1;
    } else {
      if (fLeptons.at(ind2)->Pt() > fLeptons.at(ind1)->Pt()) {
        fCandidateLeptonIndex.first = ind2;
        fCandidateLeptonIndex.second = ind1;
      } else {
        fCandidateLeptonIndex.first = ind1;
        fCandidateLeptonIndex.second = ind2;
      }
    }
  }

  if (fCandidateLeptonIndex.first != fCandidateLeptonIndex.second &&
      !(fCandidateLeptonIndex.first < 0) && !(fCandidateLeptonIndex.second < 0) &&
      fCandidateLeptonIndex.first < (int)fLeptons.size() &&
      fCandidateLeptonIndex.second < (int)fLeptons.size()) {
    const unsigned int wFlavor = fLeptons.at(fCandidateLeptonIndex.first)->GetPdgId();
    const unsigned int fakeFlavor = fLeptons.at(fCandidateLeptonIndex.second)->GetPdgId();
    if (wFlavor == 11) {
      if      (fakeFlavor == 11)  fFinalState = WEleFakeEle;
      else if (fakeFlavor == 13)  fFinalState = WEleFakeMu;
    } else if (wFlavor == 13) {
      if      (fakeFlavor == 11)  fFinalState = WMuFakeEle;
      else if (fakeFlavor == 13)  fFinalState = WMuFakeMu;
    }
  }

  return passed;
}

/*
bool
Event::PassesPreselection()
{
  if      (!(fSelectionLevel < Preselection))     return true;
  else if (fSelectionLevel == FailsPreselection)  return false;

  bool passed = false;
  fSelectionLevel = FailsPreselection;

  if((nEle + nMu) < N_LEPTONS)  return passed;

  unsigned int nWLeptons = 0;
  vector<unsigned int> indexWLeptons;
  for (vector<Lepton*>::iterator lIt = fLeptons.begin(); lIt != fLeptons.end(); ++lIt) {
    if ((*lIt)->PassesPtCut() && (*lIt)->PassesEtaCut() && (*lIt)->IsTight()) {
      nWLeptons++;
      unsigned int index = distance(fLeptons.begin(), lIt);
      indexWLeptons.push_back(index);
    }
  }

  if (nWLeptons != indexWLeptons.size()) {
    cout << "Error: Number of W Leptons DIFFERS from number of W Lepton indices !!!" << endl;
    cout << nWLeptons << " != " << indexWLeptons.size() << endl;
    return passed;
  }

  if (nWLeptons != 1)  return passed;
  else  fCandidateLeptonIndex.first = indexWLeptons.at(0);

  unsigned int nFakeLeptons = 0;
  vector<unsigned int> indexFakeLeptons;
  for (vector<Lepton*>::iterator lIt = fLeptons.begin(); lIt != fLeptons.end(); ++lIt) {
    if ((*lIt)->PassesPtCut() && (*lIt)->PassesEtaCut() && (*lIt)->IsFOMLoose()) {
      unsigned int index = distance(fLeptons.begin(), lIt);
      if ((int)index == fCandidateLeptonIndex.first)  continue;
      else {
        nFakeLeptons++;
        indexFakeLeptons.push_back(index);
      }
    }
  }

  if (nFakeLeptons != indexFakeLeptons.size()) {
    cout << "Error: Number of Fake Leptons DIFFERS from number of Fake Lepton indices !!!" << endl;
    cout << nFakeLeptons << " != " << indexFakeLeptons.size() << endl;
    return passed;
  }

  if (nFakeLeptons != 1)  return passed;
  else  fCandidateLeptonIndex.second = indexFakeLeptons.at(0);

// ### put the condition on number of selected leptons here!!! ###
  if (indexWLeptons.size() + indexFakeLeptons.size() == N_LEPTONS &&
      nWLeptons == nFakeLeptons &&
      fCandidateLeptonIndex.first != fCandidateLeptonIndex.second &&
      !(fCandidateLeptonIndex.first < 0) && !(fCandidateLeptonIndex.second < 0) &&
      fCandidateLeptonIndex.first < (int)fLeptons.size() &&
      fCandidateLeptonIndex.second < (int)fLeptons.size())
    passed = true;

  if (passed) {
    fSelectionLevel = Preselection;
    const unsigned int wFlavor = fLeptons.at(fCandidateLeptonIndex.first)->GetPdgId();
    const unsigned int fakeFlavor = fLeptons.at(fCandidateLeptonIndex.second)->GetPdgId();
    if (wFlavor == 11) {
      if      (fakeFlavor == 11)  fFinalState = WEleFakeEle;
      else if (fakeFlavor == 13)  fFinalState = WEleFakeMu;
    } else if (wFlavor == 13) {
      if      (fakeFlavor == 11)  fFinalState = WMuFakeEle;
      else if (fakeFlavor == 13)  fFinalState = WMuFakeMu;
    }
  }

  return passed;
}
*/

bool
Event::PassesSSSelection()
{
  if (!(fSelectionLevel < SSSelection))  return true;

  bool passed = false;
  if (!PassesPreselection())  return passed;

  const double wCharge = fLeptons.at(fCandidateLeptonIndex.first)->GetCharge();
  const double fakeCharge = fLeptons.at(fCandidateLeptonIndex.second)->GetCharge();
  if (wCharge != fakeCharge)  return passed;
  else  passed = true;

  if (passed)  fSelectionLevel = SSSelection;

  return passed;
}


bool
Event::PassesOFSelection()
{
  if (!(fSelectionLevel < OFSelection))  return true;

  bool passed = false;
  if (!PassesPreselection() || !PassesSSSelection())  return passed;

  const unsigned int wFlavor = fLeptons.at(fCandidateLeptonIndex.first)->GetPdgId();
  const unsigned int fakeFlavor = fLeptons.at(fCandidateLeptonIndex.second)->GetPdgId();
  if (wFlavor == fakeFlavor)  return passed;
  else  passed = true;

  if (passed)  fSelectionLevel = OFSelection;

  return passed;
}


bool
Event::PassesFullSelection()
{
  if (!(fSelectionLevel < FullSelection))  return true;

  bool passed = false;
  if (!PassesPreselection() || !PassesSSSelection() || !PassesOFSelection())  return passed;

  const double wPt = fLeptons.at(fCandidateLeptonIndex.first)->Pt();
  const double fakePt = fLeptons.at(fCandidateLeptonIndex.second)->Pt();

  const double pxMET = pfMET * cos(pfMETPhi);
  const double pyMET = pfMET * sin(pfMETPhi);
  TLorentzVector lMET(pxMET, pyMET, 0., pfMET);
  const double deltaPhiWMET = fLeptons.at(fCandidateLeptonIndex.first)->DeltaPhi(lMET);
  const double mt = sqrt(2 * pfMET * wPt * (1 - cos(deltaPhiWMET)));
  const double mass2L = (*(fLeptons.at(fCandidateLeptonIndex.first)) +
                         *(fLeptons.at(fCandidateLeptonIndex.second))).M();
  const double deltaRWFake =
    fLeptons.at(fCandidateLeptonIndex.first)->DeltaR(*(fLeptons.at(fCandidateLeptonIndex.second)));

  if (wPt > WLEPTON_PTCUT && fakePt > FAKELEPTON_PTCUT && deltaRWFake > WFAKE_DELTARCUT &&
      pfMET > METCUT && mt > MTCUT && mass2L > MASS2LCUT)
    passed = true;

  if (passed)  fSelectionLevel = FullSelection;

  return passed;
}


void
Event::Dump(ostream& out, int verbosity)
{
  out << run << ":" << lumis << ":" << event;

// Format for event listing:
// wl_pt:wl_eta:wl_phi:wl_iso:
// fakel_pt:fakel_eta:fakel_phi:fakel_iso:
// deltaR_wl_fakel:MET:MET_phi:Mt:mass(2l)

  if (verbosity > 0) {
    if (fCandidateLeptonIndex.first != fCandidateLeptonIndex.second &&
        !(fCandidateLeptonIndex.first < 0) && !(fCandidateLeptonIndex.second < 0) &&
        fCandidateLeptonIndex.first < (int)fLeptons.size() &&
        fCandidateLeptonIndex.second < (int)fLeptons.size()) {
      Lepton* wLepton = fLeptons.at(fCandidateLeptonIndex.first);
      Lepton* fakeLepton = fLeptons.at(fCandidateLeptonIndex.second);

      const double pxMET = pfMET * cos(pfMETPhi);
      const double pyMET = pfMET * sin(pfMETPhi);
      TLorentzVector lMET(pxMET, pyMET, 0., pfMET);
      const double mt = sqrt(2 * pfMET * wLepton->Pt() * (1 - cos(wLepton->DeltaPhi(lMET))));
      const double mass2L = (*wLepton + *fakeLepton).M();
      const double deltaRWFake = wLepton->DeltaR(*fakeLepton);

      out << fixed << setprecision(4)
          << ":" << wLepton->Pt()
          << ":" << wLepton->Eta()
          << ":" << wLepton->Phi()
          << ":" << wLepton->GetRelIso()
          << ":" << fakeLepton->Pt()
          << ":" << fakeLepton->Eta()
          << ":" << fakeLepton->Phi()
          << ":" << fakeLepton->GetRelIso()
          << ":" << deltaRWFake << ":" << pfMET << ":" << pfMETPhi << ":" << mt << ":" << mass2L;
    }
  }

  out << endl;
}

