#include "Event.h"
#include "Constants.h"


using namespace std;


// Initialize static members
Event* Particle::fTree = 0;


void
Particle::SetTree(Event* tree)
{
  fTree = tree;
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
  if (fTree->jetPFLooseId->at(fIndex))  loose = true;
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
  const double absEleSCEta = abs(fTree->eleSCEta->at(index));
// Slide 2 in https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf
// relIsoWithEA = 1/pt * (pfIso.sumChargedHadronPt +
//                        max(0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * EffArea)
  const double relIsoEffA = (fTree->elePFChIso->at(index) + max(fTree->elePFPhoIso->at(index) +
                             fTree->elePFNeuIso->at(index) - fTree->rho * EffArea(absEleSCEta), 0.0))
                            / fTree->elePt->at(index);
  fRelIso = relIsoEffA;
}


bool
Electron::IsLoose()
{
  bool loose = false;

  if (fTree == 0) {
    cout << "Tree pointer is ZERO!!!! \n";
    return loose;
  }

  if (!(fTree->nEle))  return loose;

  const double absEleSCEta = abs(fTree->eleSCEta->at(fIndex));
  if (absEleSCEta <= ETASCBARREL) {
    if (fTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_BARREL_LOOSE &&
        abs(fTree->eledEtaAtVtx->at(fIndex)) < DETAIN_BARREL_LOOSE &&
        abs(fTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_BARREL_LOOSE &&
        fTree->eleHoverE->at(fIndex) < HOVERE_BARREL_LOOSE &&
        fRelIso < RELISO_BARREL_LOOSE &&
        abs(fTree->eleEoverPInv->at(fIndex)) < OOEMOOP_BARREL_LOOSE &&
        abs(fTree->eleD0->at(fIndex)) < D0_BARREL_LOOSE &&
        abs(fTree->eleDz->at(fIndex)) < DZ_BARREL_LOOSE &&
        fTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_BARREL &&
        fTree->eleConvVeto->at(fIndex) == true)
      loose = true;
  } else if (absEleSCEta > ETASCBARREL && absEleSCEta < ETASCENDCAP) {
    if (fTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_ENDCAP_LOOSE &&
        abs(fTree->eledEtaAtVtx->at(fIndex)) < DETAIN_ENDCAP_LOOSE &&
        abs(fTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_ENDCAP_LOOSE &&
        fTree->eleHoverE->at(fIndex) < HOVERE_ENDCAP_LOOSE &&
        fRelIso < RELISO_ENDCAP_LOOSE &&
        abs(fTree->eleEoverPInv->at(fIndex)) < OOEMOOP_ENDCAP_LOOSE &&
        abs(fTree->eleD0->at(fIndex)) < D0_ENDCAP_LOOSE &&
        abs(fTree->eleDz->at(fIndex)) < DZ_ENDCAP_LOOSE &&
        fTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_ENDCAP &&
        fTree->eleConvVeto->at(fIndex) == true)
      loose = true;
  }

  bool looseVID = false;
  if (fTree->eleIDbit->at(fIndex)>>ELELOOSE_BIT & 1)  looseVID = true;
  if (loose != looseVID) {
    cout << "Error: VID different from Cut Based for Ele LOOSE !!!" << endl;
    cout << "VID : " << looseVID << ", Cut Based : " << loose << endl;
  }

  return loose;
}


bool
Electron::IsMedium()
{
  bool medium = false;

  if (fTree == 0) {
    cout << "Tree pointer is ZERO!!!! \n";
    return medium;
  }

  if (!(fTree->nEle))  return medium;

  const double absEleSCEta = abs(fTree->eleSCEta->at(fIndex));
  if (absEleSCEta <= ETASCBARREL) {
    if (fTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_BARREL_MEDIUM &&
        abs(fTree->eledEtaAtVtx->at(fIndex)) < DETAIN_BARREL_MEDIUM &&
        abs(fTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_BARREL_MEDIUM &&
        fTree->eleHoverE->at(fIndex) < HOVERE_BARREL_MEDIUM &&
        fRelIso < RELISO_BARREL_MEDIUM &&
        abs(fTree->eleEoverPInv->at(fIndex)) < OOEMOOP_BARREL_MEDIUM &&
        abs(fTree->eleD0->at(fIndex)) < D0_BARREL_MEDIUM &&
        abs(fTree->eleDz->at(fIndex)) < DZ_BARREL_MEDIUM &&
        fTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_BARREL &&
        fTree->eleConvVeto->at(fIndex) == true)
      medium = true;
  } else if (absEleSCEta > ETASCBARREL && absEleSCEta < ETASCENDCAP) {
    if (fTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_ENDCAP_MEDIUM &&
        abs(fTree->eledEtaAtVtx->at(fIndex)) < DETAIN_ENDCAP_MEDIUM &&
        abs(fTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_ENDCAP_MEDIUM &&
        fTree->eleHoverE->at(fIndex) < HOVERE_ENDCAP_MEDIUM &&
        fRelIso < RELISO_ENDCAP_MEDIUM &&
        abs(fTree->eleEoverPInv->at(fIndex)) < OOEMOOP_ENDCAP_MEDIUM &&
        abs(fTree->eleD0->at(fIndex)) < D0_ENDCAP_MEDIUM &&
        abs(fTree->eleDz->at(fIndex)) < DZ_ENDCAP_MEDIUM &&
        fTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_ENDCAP &&
        fTree->eleConvVeto->at(fIndex) == true)
      medium = true;
  }

  bool mediumVID = false;
  if (fTree->eleIDbit->at(fIndex)>>ELEMEDIUM_BIT & 1)  mediumVID = true;
  if (medium != mediumVID) {
    cout << "Error: VID different from Cut Based for Ele MEDIUM !!!" << endl;
    cout << "VID : " << mediumVID << "; Cut Based : " << medium << endl;
  }

  return medium;
}


bool
Electron::IsTight()
{
  bool tight = false;

  if (fTree == 0) {
    cout << "Tree pointer is ZERO!!!! \n";
    return tight;
  }

  if (!(fTree->nEle))  return tight;

  const double absEleSCEta = abs(fTree->eleSCEta->at(fIndex));
  if (absEleSCEta <= ETASCBARREL) {
    if (fTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_BARREL_TIGHT &&
        abs(fTree->eledEtaAtVtx->at(fIndex)) < DETAIN_BARREL_TIGHT &&
        abs(fTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_BARREL_TIGHT &&
        fTree->eleHoverE->at(fIndex) < HOVERE_BARREL_TIGHT &&
        fRelIso < RELISO_BARREL_TIGHT &&
        abs(fTree->eleEoverPInv->at(fIndex)) < OOEMOOP_BARREL_TIGHT &&
        abs(fTree->eleD0->at(fIndex)) < D0_BARREL_TIGHT &&
        abs(fTree->eleDz->at(fIndex)) < DZ_BARREL_TIGHT &&
        fTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_BARREL &&
        fTree->eleConvVeto->at(fIndex) == true)
      tight = true;
  } else if (absEleSCEta > ETASCBARREL && absEleSCEta < ETASCENDCAP) {
    if (fTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_ENDCAP_TIGHT &&
        abs(fTree->eledEtaAtVtx->at(fIndex)) < DETAIN_ENDCAP_TIGHT &&
        abs(fTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_ENDCAP_TIGHT &&
        fTree->eleHoverE->at(fIndex) < HOVERE_ENDCAP_TIGHT &&
        fRelIso < RELISO_ENDCAP_TIGHT &&
        abs(fTree->eleEoverPInv->at(fIndex)) < OOEMOOP_ENDCAP_TIGHT &&
        abs(fTree->eleD0->at(fIndex)) < D0_ENDCAP_TIGHT &&
        abs(fTree->eleDz->at(fIndex)) < DZ_ENDCAP_TIGHT &&
        fTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_ENDCAP &&
        fTree->eleConvVeto->at(fIndex) == true)
      tight = true;
  }

  bool tightVID = false;
  if (fTree->eleIDbit->at(fIndex)>>ELETIGHT_BIT & 1)  tightVID = true;
  if (tight != tightVID) {
    cout << "Error: VID different from Cut Based for Ele TIGHT !!!" << endl;
    cout << "VID : " << tightVID << "; Cut Based : " << tight << endl;
  }

  return tight;
}


bool
Electron::IsFOMLoose()
{
  bool fomLoose = false;

  if (fTree == 0) {
    cout << "Tree pointer is ZERO!!!! \n";
    return fomLoose;
  }

  if (!(fTree->nEle))  return fomLoose;

  const double absEleSCEta = abs(fTree->eleSCEta->at(fIndex));
  if (absEleSCEta <= ETASCBARREL) {
    if (fTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_BARREL_MEDIUM &&
        abs(fTree->eledEtaAtVtx->at(fIndex)) < DETAIN_BARREL_MEDIUM &&
        abs(fTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_BARREL_MEDIUM &&
        fTree->eleHoverE->at(fIndex) < HOVERE_BARREL_MEDIUM &&
        abs(fTree->eleEoverPInv->at(fIndex)) < OOEMOOP_BARREL_MEDIUM &&
        abs(fTree->eleD0->at(fIndex)) < D0_BARREL_MEDIUM &&
        abs(fTree->eleDz->at(fIndex)) < DZ_BARREL_MEDIUM &&
        fTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_BARREL &&
        fTree->eleConvVeto->at(fIndex) == true &&
        fRelIso < RELISO_BARREL_LOOSE)
      fomLoose = true;
  } else if (absEleSCEta > ETASCBARREL && absEleSCEta < ETASCENDCAP) {
    if (fTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_ENDCAP_MEDIUM &&
        abs(fTree->eledEtaAtVtx->at(fIndex)) < DETAIN_ENDCAP_MEDIUM &&
        abs(fTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_ENDCAP_MEDIUM &&
        fTree->eleHoverE->at(fIndex) < HOVERE_ENDCAP_MEDIUM &&
        abs(fTree->eleEoverPInv->at(fIndex)) < OOEMOOP_ENDCAP_MEDIUM &&
        abs(fTree->eleD0->at(fIndex)) < D0_ENDCAP_MEDIUM &&
        abs(fTree->eleDz->at(fIndex)) < DZ_ENDCAP_MEDIUM &&
        fTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_ENDCAP &&
        fTree->eleConvVeto->at(fIndex) == true &&
        fRelIso < RELISO_ENDCAP_LOOSE)
      fomLoose = true;
  }

  return fomLoose;
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
  const double relIsoDeltaB = (fTree->muPFChIso->at(index) + max(fTree->muPFPhoIso->at(index) +
                               fTree->muPFNeuIso->at(index) - 0.5 * fTree->muPFPUIso->at(index), 0.0))
                              / fTree->muPt->at(index);
  fRelIso = relIsoDeltaB;
}


bool
Muon::IsLoose()
{
  bool loose = false;

  if (fTree == 0) {
    cout << "Tree pointer is ZERO!!!! \n";
    return loose;
  }

  if (!(fTree->nMu))  return loose;

  bool looseCut = false;
  if ((fTree->muType->at(fIndex)>>GLOBALMUON_BIT & 1 || fTree->muType->at(fIndex)>>TRACKERMUON_BIT & 1) &&
      fTree->muType->at(fIndex)>>PFMUON_BIT & 1)
    looseCut = true;

/*
  bool looseVID = false;
  if (fTree->muIsLooseID->at(fIndex))  looseVID = true;
  if (looseCut != looseVID) {
    cout << "Error: VID different from Cut Based for Mu LOOSE !!!" << endl;
    cout << "VID : " << looseVID << ", Cut Based : " << looseCut << endl;
  }
*/

  if (looseCut && fRelIso < MU_RELISO_LOOSE)  loose = true;

  return loose;
}


bool
Muon::IsMedium()
{
  bool mediumCut = false;

  if (fTree == 0) {
    cout << "Tree pointer is ZERO!!!! \n";
    return mediumCut;
  }

  if (!(fTree->nEle))  return mediumCut;

  bool good = false;
  if (fTree->muType->at(fIndex)>>GLOBALMUON_BIT & 1 &&
      fTree->muChi2NDF->at(fIndex) < 3. &&
      fTree->muchi2LocalPosition->at(fIndex) < 12. &&
      fTree->mutrkKink->at(fIndex) < 20. &&
      fTree->musegmentCompatibility->at(fIndex) > 0.303)
    good = true;

  if (IsLoose() && fTree->muInnervalidFraction->at(fIndex) > 0.8 &&
      (good || fTree->musegmentCompatibility->at(fIndex) > 0.451))
    mediumCut = true;
/*
  bool mediumVID = false;
  if (fTree->muIsMediumID->at(fIndex))  mediumVID = true;
  if (mediumCut != mediumVID) {
    cout << "Error: VID different from Cut Based for Mu MEDIUM !!!" << endl;
    cout << "VID : " << mediumVID << ", Cut Based : " << mediumCut << endl;
  }
*/
  return mediumCut;
}


bool
Muon::IsTight()
{
  bool tight = false;

  if (fTree == 0) {
    cout << "Tree pointer is ZERO!!!! \n";
    return tight;
  }

  if (!(fTree->nMu))  return tight;

  bool tightCut = false;
  if (fTree->muChi2NDF->at(fIndex) < 10. &&
      fTree->muMuonHits->at(fIndex) > 0 &&
      fTree->muStations->at(fIndex) > 1 &&
      abs(fTree->muD0->at(fIndex)) < 0.2 &&
      abs(fTree->muDz->at(fIndex)) < 0.5 &&
      fTree->muPixelHits->at(fIndex) > 0 &&
      fTree->muTrkLayers->at(fIndex) > 5 &&
      fTree->muType->at(fIndex)>>GLOBALMUON_BIT & 1 &&
      fTree->muType->at(fIndex)>>PFMUON_BIT & 1)
    tightCut = true;

/*
  bool tightVID = false;
  if (fTree->muIsTightID->at(fIndex))  tightVID = true;
  if (tightCut != tightVID) {
    cout << "Error: VID different from Cut Based for Mu TIGHT !!!" << endl;
    cout << "VID : " << tightVID << ", Cut Based : " << tightCut << endl;
  }
*/

  if (tightCut && fRelIso < MU_RELISO_TIGHT)  tight = true;

  return tight;
}


bool
Muon::IsFOMLoose()
{
  bool fomLoose = false;

  if (fTree == 0) {
    cout << "Tree pointer is ZERO!!!! \n";
    return fomLoose;
  }

  if (!(fTree->nMu))  return fomLoose;

  if (fTree->muChi2NDF->at(fIndex) < 10. &&
      fTree->muMuonHits->at(fIndex) > 0 &&
      fTree->muStations->at(fIndex) > 1 &&
      abs(fTree->muD0->at(fIndex)) < 0.2 &&
      abs(fTree->muDz->at(fIndex)) < 0.5 &&
      fTree->muPixelHits->at(fIndex) > 0 &&
      fTree->muTrkLayers->at(fIndex) > 5 &&
      fTree->muType->at(fIndex)>>GLOBALMUON_BIT & 1 &&
      fTree->muType->at(fIndex)>>PFMUON_BIT & 1 &&
      fRelIso < MU_RELISO_LOOSE)
    fomLoose = true;

  return fomLoose;
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

