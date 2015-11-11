#include "WZEvent.h"
#include "Constants.h"


using namespace std;


// Initialize static members
WZEvent* Particle::fWZTree = 0;


void
Particle::SetWZTree(WZEvent* wztree)
{
  fWZTree = wztree;
}


Particle::Particle(unsigned int index, double pt, double eta, double phi) 
{
  SetPtEtaPhiM(pt, eta, phi, 0.);
  fIndex = index;
}


Particle::Particle(unsigned int index, double pt, double eta, double phi, double en)
{
  SetPtEtaPhiE(pt, eta, phi, en);
  fIndex = index;
}


Jet::Jet(unsigned int index, double pt, double eta, double phi, double en)
  : Particle(index, pt, eta, phi, en)
{

}


bool
Jet::IsLoose()
{
  bool loose = false;
  if (fWZTree->jetPFLooseId->at(fIndex))  loose = true;
  return loose;
}


bool
Jet::PassesPtCut()
{
  bool passPt = false;
  if (Pt() > JET_PTCUT)  passPt = true;
  return passPt;
}


bool
Jet::PassesEtaCut()
{
  bool passEta = false;
  if (abs(Eta()) > JET_ETACUT)  passEta = true;
  return passEta;
}


Lepton::Lepton(unsigned int index, double pt, double eta, double phi, double charge)
  : Particle(index, pt, eta, phi) 
{
  fCharge = charge;
  fScaleFactor = 0.; // Waiting for real implementation
}


double
EffArea(double absEleSCEta)
{
  double effA = 0;
// Slide 4 in https://indico.cern.ch/event/370507/contribution/1/attachments/1140657/1633761/Rami_eleCB_ID_25ns.pdf
// and Slide 12 in https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf
// and Line 407 for abs(eleSCEta) in https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_747/ElectronNtupler/plugins/SimpleElectronNtupler.cc
  if (absEleSCEta >= 0.0 && absEleSCEta < 1.0)         effA = 0.1752;
  else if (absEleSCEta >= 1.0 && absEleSCEta < 1.479)  effA = 0.1862;
  else if (absEleSCEta >= 1.479 && absEleSCEta < 2.0)  effA = 0.1411;
  else if (absEleSCEta >= 2.0 && absEleSCEta < 2.2)    effA = 0.1534;
  else if (absEleSCEta >= 2.2 && absEleSCEta < 2.3)    effA = 0.1903;
  else if (absEleSCEta >= 2.3 && absEleSCEta < 2.4)    effA = 0.2243;
  else if (absEleSCEta >= 2.4 && absEleSCEta < 2.5)    effA = 0.2687;
  else  effA = 0;

  return effA;
}


Electron::Electron(unsigned int index, double pt, double eta, double phi, double charge)
  : Lepton(index, pt, eta, phi, charge)
{
  fPdgId = 11;
  const double absEleSCEta = abs(fWZTree->eleSCEta->at(index));
// Slide 2 in https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf
// relIsoWithEA = 1/pt * (pfIso.sumChargedHadronPt +
//                        max(0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * EffArea)
  const double relIsoEffA = (fWZTree->elePFChIso->at(index) + max(fWZTree->elePFPhoIso->at(index) +
                             fWZTree->elePFNeuIso->at(index) - fWZTree->rho * EffArea(absEleSCEta), 0.0))
                            / fWZTree->elePt->at(index);
  fRelIso = relIsoEffA;
}


bool
Electron::IsLoose()
{
  bool loose = false;

  if (fWZTree == 0) {
    cout << "WZTree pointer is ZERO!!!! \n";
    return loose;
  }

  if (!(fWZTree->nEle))  return loose;

  const double absEleSCEta = abs(fWZTree->eleSCEta->at(fIndex));
  if (absEleSCEta <= ETASCBARREL) {
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_BARREL_LOOSE &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_BARREL_LOOSE &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_BARREL_LOOSE &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_BARREL_LOOSE &&
        fRelIso < RELISO_BARREL_LOOSE &&
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_BARREL_LOOSE &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_BARREL_LOOSE &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_BARREL_LOOSE &&
        fWZTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_BARREL &&
        fWZTree->eleConvVeto->at(fIndex) == true)
      loose = true;
  } else if (absEleSCEta > ETASCBARREL && absEleSCEta < ETASCENDCAP) {
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_ENDCAP_LOOSE &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_ENDCAP_LOOSE &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_ENDCAP_LOOSE &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_ENDCAP_LOOSE &&
        fRelIso < RELISO_ENDCAP_LOOSE &&
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_ENDCAP_LOOSE &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_ENDCAP_LOOSE &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_ENDCAP_LOOSE &&
        fWZTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_ENDCAP &&
        fWZTree->eleConvVeto->at(fIndex) == true)
      loose = true;
  }

  bool looseVID = false;
  if (fWZTree->eleIDbit->at(fIndex)>>ELELOOSE_BIT&1)  looseVID = true;
  if (loose != looseVID) {
    cout << "Error: VID different from Cut Based for Ele LOOSE !!!" << endl;
    cout << "VID : " << looseVID << ", Cut Based : " << loose << endl;
  }

  return loose;
}


bool
Electron::IsTight()
{
  bool tight = false;

  if (fWZTree == 0) {
    cout << "WZTree pointer is ZERO!!!! \n";
    return tight;
  }

  if (!(fWZTree->nEle))  return tight;

  const double absEleSCEta = abs(fWZTree->eleSCEta->at(fIndex));
  if (absEleSCEta <= ETASCBARREL) {
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_BARREL_MEDIUM &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_BARREL_MEDIUM &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_BARREL_MEDIUM &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_BARREL_MEDIUM &&
        fRelIso < RELISO_BARREL_MEDIUM &&
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_BARREL_MEDIUM &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_BARREL_MEDIUM &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_BARREL_MEDIUM &&
        fWZTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_BARREL &&
        fWZTree->eleConvVeto->at(fIndex) == true)
      tight = true;
  } else if (absEleSCEta > ETASCBARREL && absEleSCEta < ETASCENDCAP) {
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_ENDCAP_MEDIUM &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_ENDCAP_MEDIUM &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_ENDCAP_MEDIUM &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_ENDCAP_MEDIUM &&
        fRelIso < RELISO_ENDCAP_MEDIUM &&
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_ENDCAP_MEDIUM &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_ENDCAP_MEDIUM &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_ENDCAP_MEDIUM &&
        fWZTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_ENDCAP &&
        fWZTree->eleConvVeto->at(fIndex) == true)
      tight = true;
  }

  bool mediumVID = false;
  if (fWZTree->eleIDbit->at(fIndex)>>ELEMEDIUM_BIT&1)  mediumVID = true;
  if (tight != mediumVID) {
    cout << "Error: VID different from Cut Based for Ele MEDIUM !!!" << endl;
    cout << "VID : " << mediumVID << ", Cut Based: " << tight << endl;
  }

  return tight;
}


bool
Electron::PassesPtCut()
{
  bool passPt = false;
  if (Pt() > ELE_PTCUT)  passPt = true;
  return passPt;
}


bool
Electron::PassesEtaCut()
{
  bool passEta = false;
  if (abs(Eta()) < ELE_ETACUT)  passEta = true;
  return passEta;
}


Muon::Muon(unsigned int index, double pt, double eta, double phi, double charge)
  : Lepton(index, pt, eta, phi, charge)
{
  fPdgId = 13;
// From UW Twiki https://twiki.cern.ch/twiki/bin/viewauth/CMS/WZ13TeV
// relIsoDeltaBeta = 1/pt * (mu.chargedHadronIso() +
//                           max(mu.photonIso() + mu.neutralHadronIso() - 0.5 * mu.puChargedHadronIso,0.0))
  const double relIsoDeltaB = (fWZTree->muPFChIso->at(index) + max(fWZTree->muPFPhoIso->at(index) +
                               fWZTree->muPFNeuIso->at(index) - 0.5 * fWZTree->muPFPUIso->at(index), 0.0))
                              / fWZTree->muPt->at(index);
  fRelIso = relIsoDeltaB;
}


bool
Muon::IsLoose()
{
  bool loose = false;

  if (fWZTree == 0) {
    cout << "Tree pointer is ZERO!!!! \n";
    return loose;
  }

  if (!(fWZTree->nMu))  return loose;

  bool looseCut = false;
  if ((fWZTree->muType->at(fIndex)>>GLOBALMUON_BIT & 1 || fWZTree->muType->at(fIndex)>>TRACKERMUON_BIT & 1) &&
      fWZTree->muType->at(fIndex)>>PFMUON_BIT & 1)
    looseCut = true;

  bool looseVID = false;
  if (fWZTree->muIsLooseID->at(fIndex))  looseVID = true;
  if (looseCut != looseVID) {
    cout << "Error: VID different from Cut Based for Mu LOOSE !!!" << endl;
    cout << "VID : " << looseVID << ", Cut Based : " << looseCut << endl;
  }

  if (looseCut && fRelIso < MU_RELISO_LOOSE)  loose = true;

  return loose;
}


bool
Muon::IsTight()
{
  bool tight = false;

  if (fWZTree == 0) {
    cout << "Tree pointer is ZERO!!!! \n";
    return tight;
  }

  if (!(fWZTree->nMu))  return tight;

  bool tightCut = false;
  if (fWZTree->muChi2NDF->at(fIndex) < 10. &&
      fWZTree->muMuonHits->at(fIndex) > 0 &&
      fWZTree->muStations->at(fIndex) > 1 &&
      abs(fWZTree->muD0->at(fIndex)) < 0.2 &&
      abs(fWZTree->muDz->at(fIndex)) < 0.5 &&
      fWZTree->muPixelHits->at(fIndex) > 0 &&
      fWZTree->muTrkLayers->at(fIndex) > 5 &&
      fWZTree->muType->at(fIndex)>>GLOBALMUON_BIT & 1 &&
      fWZTree->muType->at(fIndex)>>PFMUON_BIT & 1)
    tightCut = true;

  bool tightVID = false;
  if (fWZTree->muIsTightID->at(fIndex))  tightVID = true;
  if (tightCut != tightVID) {
    cout << "Error: VID different from Cut Based for Mu TIGHT !!!" << endl;
    cout << "VID : " << tightVID << ", Cut Based : " << tightCut << endl;
  }

  if (tightCut && fRelIso < MU_RELISO_TIGHT)  tight = true;

  return tight;
}


bool
Muon::PassesPtCut()
{
  bool passPt = false;
  if (Pt() > MU_PTCUT)  passPt = true;
  return passPt;
}


bool
Muon::PassesEtaCut()
{
  bool passEta = false;
  if (abs(Eta()) < MU_ETACUT)  passEta = true;
  return passEta;
}

