#include "WZEvent.h"
#include "Constants.h"

#include <iomanip>


using namespace std;


WZEvent::WZEvent(TTree* tree) :
  WZBASECLASS(tree)
{
// Set static pointer to Particle class
  Particle::SetWZTree(this);
}


void WZEvent::Clear()
{
  fHLT25ns.clear();

  fFinalState = undefined;
  fSelectionLevel = Undefined;

// Empty leptons
  for (vector<Lepton*>::iterator lIt = fLeptons.begin(); lIt != fLeptons.end(); ) {
    delete *lIt;  
    lIt = fLeptons.erase(lIt);
  }

  fTightLeptonsIndex.clear();
  fZLeptonsIndex = make_pair(-1, -1);
  fWLeptonIndex = -1;
}


void WZEvent::Read()
{
  Clear();

// in ggNtuplizer versions V07-04-14-00
  if (HLTEleMuX) {
    const vector<unsigned int> hlt25nsBits { 7, 20, 21, 50, 51 };
        for (vector<unsigned int>::const_iterator bIt = hlt25nsBits.begin();
         bIt != hlt25nsBits.end(); ++bIt) {
      (HLTEleMuX>>(*bIt) & 1)  ?  fHLT25ns.push_back(true)  :  fHLT25ns.push_back(false);
    }
  }

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
WZEvent::PassesTrigger()
{
  bool passed = false;

  for (vector<bool>::const_iterator bIt = fHLT25ns.begin(); bIt != fHLT25ns.end(); ++bIt) {
    if (*bIt)  passed = true;
    else  continue;
  }

  return passed;
}


bool WZEvent::PassesPreselection(SelectionType type)
{
  if      (!(fSelectionLevel < Preselection))     return true;
  else if (fSelectionLevel == FailsPreselection)  return false;

  bool passed = false;
  fSelectionLevel = FailsPreselection;

  if((nEle + nMu) < 3)  return passed;

  unsigned int nEleTight = 0;
  unsigned int nMuTight = 0;

// ## Nominal Selection ##
  if (type == Nominal)
  {
    for (vector<Lepton*>::iterator lIt = fLeptons.begin(); lIt != fLeptons.end(); ++lIt) {
      if ((*lIt)->PassesPtCut() && (*lIt)->PassesEtaCut() && (*lIt)->IsTight()) {
        unsigned int index = distance(fLeptons.begin(), lIt);
        fTightLeptonsIndex.push_back(index);
        if ((*lIt)->GetPdgId() == 11)       nEleTight++;
        else if ((*lIt)->GetPdgId() == 13)  nMuTight++;
      }
    }

    if (nMuTight + nEleTight != fTightLeptonsIndex.size()) {
      cout << "ERROR: Number of Tight leptons DIFFERS from number of Tight lepton indices !!!" << endl;
      cout << nMuTight << " + " << nEleTight << " != " << fTightLeptonsIndex.size() << endl;
      return passed;
    }

    // ### put the condition on number of tight leptons here!!! ###
    if (fTightLeptonsIndex.size() == N_TIGHTLEPTONS)  passed = true;
  }

// ## Matrix Method Selection ##
  if (type == MatrixMethod)
  {
    for (vector<Lepton*>::iterator lIt = fLeptons.begin(); lIt != fLeptons.end(); ++lIt) {
      if ((*lIt)->PassesPtCut() && (*lIt)->PassesEtaCut() && (*lIt)->IsLoose()) {
        unsigned int index = distance(fLeptons.begin(), lIt);
        fTightLeptonsIndex.push_back(index);
        if ((*lIt)->GetPdgId() == 11)       nEleTight++;
        else if ((*lIt)->GetPdgId() == 13)  nMuTight++;
      }
    }

    if (nMuTight + nEleTight != fTightLeptonsIndex.size()) {
      cout << "ERROR: Number of Loose leptons DIFFERS from number of Loose lepton indices !!!" << endl;
      cout << nMuTight << " + " << nEleTight << " != " << fTightLeptonsIndex.size() << endl;
      return passed;
    }

    // ### put the condition on number of loose and tight leptons here!!! ###
    if (fTightLeptonsIndex.size() >= N_TIGHTLEPTONS) {
      unsigned int nTight = 0;
      for (vector<unsigned int>::const_iterator iIt = fTightLeptonsIndex.begin();
           iIt != fTightLeptonsIndex.end(); ++iIt)
        if (fLeptons.at(*iIt)->IsTight())  nTight++;
      if (nTight <= N_TIGHTLEPTONS)  passed = true;
    }
  }

  if (passed) {
    fSelectionLevel = Preselection;
    if      (nEleTight == 3 && nMuTight == 0)  fFinalState = eee;
    else if (nEleTight == 2 && nMuTight == 1)  fFinalState = eem;
    else if (nEleTight == 1 && nMuTight == 2)  fFinalState = mme;
    else if (nEleTight == 0 && nMuTight == 3)  fFinalState = mmm;
  }

  return passed;
}


bool WZEvent::PassesZSelection(SelectionType type)
{
  if (!(fSelectionLevel < ZSelection))  return true;

  bool passed = false;
  if (!PassesPreselection(type))  return passed;

  vector<pair<unsigned int, unsigned int> > iZl;
  vector<double> dMassZ;

  for (vector<unsigned int>::const_iterator iIt1 = fTightLeptonsIndex.begin();
       iIt1 != fTightLeptonsIndex.end(); ++iIt1) {
    for (vector<unsigned int>::const_iterator iIt2 = fTightLeptonsIndex.begin()+1;
         iIt2 != fTightLeptonsIndex.end(); ++iIt2) {
      unsigned int il1 = *iIt1;
      unsigned int il2 = *iIt2;

      if (abs(fLeptons.at(il1)->GetPdgId()) != abs(fLeptons.at(il2)->GetPdgId()) ||
          fLeptons.at(il1)->GetCharge() == fLeptons.at(il2)->GetCharge())
        continue;

      const double mZCand = (*(fLeptons.at(il1)) + *(fLeptons.at(il2))).M();
      if (mZCand < MZ_MIN || mZCand > MZ_MAX ||
          !(max(fLeptons.at(il1)->Pt(), fLeptons.at(il2)->Pt()) > ZLEADINGLEPTON_PTCUT))
        continue;
      else {
        iZl.push_back(make_pair(il1, il2));
        const double dMZ = abs(mZCand - PDG_ZMASS);
        dMassZ.push_back(dMZ);
      }
    }
  }

  if (iZl.size() == 0)  return passed;
  else {
    passed = true;
    vector<double>::iterator bestMin = min_element(dMassZ.begin(), dMassZ.end());
    unsigned int iMin = distance(dMassZ.begin(), bestMin);
    fZLeptonsIndex = iZl.at(iMin);
  }

  if (fZLeptonsIndex.first != fZLeptonsIndex.second) {
    if (fLeptons.at(fZLeptonsIndex.first)->Pt() < fLeptons.at(fZLeptonsIndex.second)->Pt())
      swap(fZLeptonsIndex.first, fZLeptonsIndex.second);
  }

  if (passed)  fSelectionLevel = ZSelection;

  return passed;
}


bool WZEvent::PassesWSelection(SelectionType type)
{
  if (!(fSelectionLevel < WSelection))  return true;

  bool passed = false;
  if (!PassesPreselection(type) || !PassesZSelection(type))  return passed;

  vector<unsigned int> iWl;
  vector<double> ptWl;

  for (vector<unsigned int>::const_iterator iIt = fTightLeptonsIndex.begin();
       iIt != fTightLeptonsIndex.end(); ++iIt) {
    if ((int)*iIt != fZLeptonsIndex.first && (int)*iIt != fZLeptonsIndex.second) {
      const double wlPt = fLeptons.at(*iIt)->Pt();
      const double deltaR1 = fLeptons.at(*iIt)->DeltaR(*(fLeptons.at(fZLeptonsIndex.first)));
      const double deltaR2 = fLeptons.at(*iIt)->DeltaR(*(fLeptons.at(fZLeptonsIndex.second)));

      if (wlPt > WLEPTON_PTCUT && deltaR1 > WZ_DELTARCUT && deltaR2 > WZ_DELTARCUT) {

// Requirement for W lepton to be TIGHT
/*
        if (fLeptons.at(*iIt)->GetPdgId() == 11) {
          const unsigned int index = fLeptons.at(*iIt)->GetIndex();
          if (!(eleIDbit->at(index)>>ELETIGHT_BIT&1))  continue;
        }
*/

        iWl.push_back(*iIt);
        ptWl.push_back(wlPt);
      } else continue;
    }
  }

  if (iWl.size() > 1)  cout << "WARNING: More than 1 W lepton - applying highest Pt criterium !!!" << endl;

  if (iWl.size() == 0)  return passed;
  else {
    passed = true;
    vector<double>::iterator bestMax = max_element(ptWl.begin(), ptWl.end());
    unsigned int iMax = distance(ptWl.begin(), bestMax);
    fWLeptonIndex = iWl.at(iMax);
  }

  if (passed)  fSelectionLevel = WSelection;

  return passed;
}


bool WZEvent::PassesFullSelection(SelectionType type)
{
  if (!(fSelectionLevel < FullSelection))  return true;

  bool passed = false;
  if (!PassesPreselection(type) || !PassesZSelection(type) || !PassesWSelection(type))  return passed;

  const double mass3L = (*(GetZLeptons().first) + *(GetZLeptons().second) + *(GetWLepton())).M();
  if (pfMET > METCUT && mass3L > MASS3LCUT &&
      fZLeptonsIndex.first != fZLeptonsIndex.second &&
      fZLeptonsIndex.first != fWLeptonIndex &&
      fZLeptonsIndex.second != fWLeptonIndex &&
      !(fZLeptonsIndex.first < 0) && fZLeptonsIndex.first < (int)fLeptons.size() &&
      !(fZLeptonsIndex.second < 0) && fZLeptonsIndex.second < (int)fLeptons.size() &&
      !(fWLeptonIndex < 0) && fWLeptonIndex < (int)fLeptons.size())
    passed = true;

  if (passed) {
    fSelectionLevel = FullSelection;
    unsigned int zFlavor = abs(fLeptons.at(fZLeptonsIndex.first)->GetPdgId());
    unsigned int wFlavor = abs(fLeptons.at(fWLeptonIndex)->GetPdgId());

    if (zFlavor == 11) {
      if      (wFlavor == 11)  fFinalState = eee;
      else if (wFlavor == 13)  fFinalState = eem;
    } else if (zFlavor == 13) {
      if      (wFlavor == 11)  fFinalState = mme;
      else if (wFlavor == 13)  fFinalState = mmm;
    }
  }

  return passed;
}


pair<Lepton*, Lepton*> WZEvent::GetZLeptons()
{
  Lepton* zL1 = fLeptons.at(fZLeptonsIndex.first);
  Lepton* zL2 = fLeptons.at(fZLeptonsIndex.second);
  pair<Lepton*, Lepton*> zL = make_pair(zL1, zL2);
  return zL;
}


void WZEvent::Dump(ostream& out, int verbosity)
{
  out << run << ":" << lumis << ":" << event;

// Preselection (verbosity = 5) format for event listing:
// l1_pt:l1_eta:l1_phi:l1_iso:l2_pt:l2_eta:l2_phi:l2_iso:l3_pt:l3_eta:l3_phi:l3_iso:
// deltaR_l1_l2:deltaR_l1_l3:deltaR_l2_l3:0:MET:MET_phi:3l_mass

  if (verbosity == 5)
  {
    const unsigned int nTight = fTightLeptonsIndex.size();
    if (nTight == N_TIGHTLEPTONS) {
      Lepton* tight1 = fLeptons.at(fTightLeptonsIndex.at(0));
      Lepton* tight2 = fLeptons.at(fTightLeptonsIndex.at(1));
      Lepton* tight3 = fLeptons.at(fTightLeptonsIndex.at(2));
      vector<Lepton*> tightLeptons = { tight1, tight2, tight3 };

      sort(tightLeptons.begin(), tightLeptons.end(), HigherPt);

      TLorentzVector total3L;
      for (vector<Lepton*>::const_iterator lIt = tightLeptons.begin();
           lIt != tightLeptons.end(); ++lIt) {
        out << fixed << setprecision(4)
            << ":" << (*lIt)->Pt() << ":" << (*lIt)->Eta() << ":" << (*lIt)->Phi()
            << ":" << (*lIt)->GetRelIso();
        total3L += *(*lIt);
      }

      for (vector<Lepton*>::const_iterator lIt1 = tightLeptons.begin();
           lIt1 != tightLeptons.end(); ++lIt1) {
        for (vector<Lepton*>::const_iterator lIt2 = tightLeptons.begin() + 1;
             lIt2 != tightLeptons.end(); ++lIt2)
          out << fixed << setprecision(4) << ":" << (*lIt1)->DeltaR(*(*lIt2));
      }

      const double mass3L = total3L.M();
      out << fixed << setprecision(4)
          << ":" << 0 << ":" << pfMET << ":" << pfMETPhi << ":" << mass3L;
    }
  }

// Event listing format for Z Selection (verbosity = 7) and W/Full/Final Selection (verbosity >= 10):
//  Zl1_pt:Zl1_eta:Zl1_phi:Zl1_iso:Zl2_pt:Zl2_eta:Zl2_phi:Zl2_iso:Wl_pt:Wl_eta:Wl_phi:Wl_iso:
//  deltaR_Zl1_Zl2:deltaR_Zl1_Wl:deltaR_Zl2_Wl:Z_mass:MET:MET_phi:3l_mass

  if (verbosity == 7)
  {
    if (fZLeptonsIndex.first < (int)fLeptons.size() && !(fZLeptonsIndex.first < 0) &&
        fZLeptonsIndex.second < (int)fLeptons.size() && !(fZLeptonsIndex.second < 0)) {
      out << fixed << setprecision(4)
          << ":" << GetZLeptons().first->Pt()
          << ":" << GetZLeptons().first->Eta()
          << ":" << GetZLeptons().first->Phi()
          << ":" << GetZLeptons().first->GetRelIso()
          << ":" << GetZLeptons().second->Pt()
          << ":" << GetZLeptons().second->Eta()
          << ":" << GetZLeptons().second->Phi()
          << ":" << GetZLeptons().second->GetRelIso();

      for (vector<unsigned int>::const_iterator iIt = fTightLeptonsIndex.begin();
           iIt != fTightLeptonsIndex.end(); ++iIt) {
        if ((int)*iIt != fZLeptonsIndex.first && (int)*iIt != fZLeptonsIndex.second) {
          out << ":" << fLeptons.at(*iIt)->Pt()
              << ":" << fLeptons.at(*iIt)->Eta()
              << ":" << fLeptons.at(*iIt)->Phi()
              << ":" << fLeptons.at(*iIt)->GetRelIso()
              << ":" << GetZLeptons().first->DeltaR(*(GetZLeptons().second))
              << ":" << GetZLeptons().first->DeltaR(*(fLeptons.at(*iIt)))
              << ":" << GetZLeptons().second->DeltaR(*(fLeptons.at(*iIt)))
              << ":" << (*(GetZLeptons().first) + *(GetZLeptons().second)).M()
              << ":" << pfMET << ":" << pfMETPhi
              << ":" << (*(GetZLeptons().first) + *(GetZLeptons().second) + *(fLeptons.at(*iIt))).M();
        }
      }
    }
  }

  if (!(verbosity < 10))
  {
    vector<int> indexWZLeptons =
      { fZLeptonsIndex.first, fZLeptonsIndex.second, fWLeptonIndex };
    for (vector<int>::iterator iIt = indexWZLeptons.begin(); iIt != indexWZLeptons.end(); ++iIt) {
      if (*iIt < (int)fLeptons.size() && !(*iIt < 0)) {
        out << fixed << setprecision(4)
            << ":" << fLeptons.at(*iIt)->Pt()
            << ":" << fLeptons.at(*iIt)->Eta()
            << ":" << fLeptons.at(*iIt)->Phi()
            << ":" << fLeptons.at(*iIt)->GetRelIso();
      }
    }
    out << fixed << setprecision(4)
        << ":" << GetZLeptons().first->DeltaR(*(GetZLeptons().second))
        << ":" << GetZLeptons().first->DeltaR(*(GetWLepton()))
        << ":" << GetZLeptons().second->DeltaR(*(GetWLepton()))
        << ":" << (*(GetZLeptons().first) + *(GetZLeptons().second)).M()
        << ":" << pfMET << ":" << pfMETPhi
        << ":" << (*(GetZLeptons().first) + *(GetZLeptons().second) + *(GetWLepton())).M();
  }

  out << endl;
}

