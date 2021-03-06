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


double
Electron::ComputeLeptonScaleFactor()
{
  double sf = 0;

  if (IsTight()) {
    if (0 <= abs(Eta()) && abs(Eta()) < 1.5) {
      if (10 < Pt() && Pt() <= 20)          sf = 1.0217;
      else if (20 < Pt() && Pt() <= 30)     sf = 0.9768;
      else if (30 < Pt() && Pt() <= 40)     sf = 0.9868;
      else if (40 < Pt() && Pt() <= 50)     sf = 0.9819;
      else if (50 < Pt() && Pt() <= 100)    sf = 0.9782;
      else if (100 < Pt() && Pt() <= 1000)  sf = 0.9844;
    } else if (1.5 <= abs(Eta()) && abs(Eta()) < 2.5) {
      if (10 < Pt() && Pt() <= 20)          sf = 0.9588;
      else if (20 < Pt() && Pt() <= 30)     sf = 0.9669;
      else if (30 < Pt() && Pt() <= 40)     sf = 0.9786;
      else if (40 < Pt() && Pt() <= 50)     sf = 0.9998;
      else if (50 < Pt() && Pt() <= 100)    sf = 1.0019;
      else if (100 < Pt() && Pt() <= 1000)  sf = 1.0297;
    }
  } else if (IsLoose()) {
    cout << "WARNING: Lepton scale factor for LOOSE electron is computed!" << endl;
    if (0 <= abs(Eta()) && abs(Eta()) < 1.5) {
      if (10 < Pt() && Pt() <= 20)          sf = 1.0270;
      else if (20 < Pt() && Pt() <= 30)     sf = 0.9791;
      else if (30 < Pt() && Pt() <= 40)     sf = 0.9918;
      else if (40 < Pt() && Pt() <= 50)     sf = 0.9987;
      else if (50 < Pt() && Pt() <= 100)    sf = 0.9851;
      else if (100 < Pt() && Pt() <= 1000)  sf = 0.9894;
    } else if (1.5 <= abs(Eta()) && abs(Eta()) < 2.5) {
      if (10 < Pt() && Pt() <= 20)          sf = 0.9886;
      else if (20 < Pt() && Pt() <= 30)     sf = 0.9806;
      else if (30 < Pt() && Pt() <= 40)     sf = 0.9983;
      else if (40 < Pt() && Pt() <= 50)     sf = 1.0095;
      else if (50 < Pt() && Pt() <= 100)    sf = 1.0053;
      else if (100 < Pt() && Pt() <= 1000)  sf = 1.1020;
    }
  }

  return sf;
}


double
Electron::ComputePromptRate()
{
  double rate = 0.;

  const double absEleSCEta = abs(fWZTree->eleSCEta->at(fIndex));
  if (absEleSCEta <= ETASCBARREL) {
    if (10 < Pt() && Pt() <= 20)        rate = 1.;
    else if (20 < Pt() && Pt() <= 30)   rate = 1.;
    else if (30 < Pt() && Pt() <= 200)  rate = 1.;
  } else if (absEleSCEta > ETASCBARREL && absEleSCEta < ETASCENDCAP) {
    if (10 < Pt() && Pt() <= 20)        rate = 1.;
    else if (20 < Pt() && Pt() <= 30)   rate = 1.;
    else if (30 < Pt() && Pt() <= 200)  rate = 1.;
  }

  return  rate;
}


double
Electron::ComputeFakeRate()
{
  double fake = 0.;
  const double absEleSCEta = abs(fWZTree->eleSCEta->at(fIndex));
  if (absEleSCEta <= ETASCBARREL) {
    if (10 < Pt() && Pt() <= 20)        fake = 0.247;
    else if (20 < Pt() && Pt() <= 30)   fake = 0.242;
    else if (30 < Pt() && Pt() <= 200)  fake = 0.197;
  } else if (absEleSCEta > ETASCBARREL && absEleSCEta < ETASCENDCAP) {
    if (10 < Pt() && Pt() <= 20)        fake = 0.274;
    else if (20 < Pt() && Pt() <= 30)   fake = 0.291;
    else if (30 < Pt() && Pt() <= 200)  fake = 0.313;
  }

  double rate = fake / (1. + fake);

  return  rate;
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

  bool looseWP = false;
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
      looseWP = true;
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
      looseWP = true;
  }

  bool looseVID = false;
  if (fWZTree->eleIDbit->at(fIndex)>>ELELOOSE_BIT&1)  looseVID = true;
  if (looseWP != looseVID) {
    cout << "WARNING: VID different from Cut Based for Ele LOOSE !!!" << endl;
    cout << "VID : " << looseVID << ", Cut Based : " << looseWP << endl;
  }

  bool looseCB = false;
  if (absEleSCEta <= ETASCBARREL) {
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_BARREL_LOOSE &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_BARREL_LOOSE &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_BARREL_LOOSE &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_BARREL_LOOSE &&
        // fRelIso < RELISO_BARREL_LOOSE &&  ->  no RelIso requirement for Loose !!!
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_BARREL_LOOSE &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_BARREL_LOOSE &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_BARREL_LOOSE &&
        fWZTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_BARREL &&
        fWZTree->eleConvVeto->at(fIndex) == true)
      looseCB = true;
  } else if (absEleSCEta > ETASCBARREL && absEleSCEta < ETASCENDCAP) {
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_ENDCAP_LOOSE &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_ENDCAP_LOOSE &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_ENDCAP_LOOSE &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_ENDCAP_LOOSE &&
        // fRelIso < RELISO_ENDCAP_LOOSE &&  ->  no RelIso requirement for Loose !!!
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_ENDCAP_LOOSE &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_ENDCAP_LOOSE &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_ENDCAP_LOOSE &&
        fWZTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_ENDCAP &&
        fWZTree->eleConvVeto->at(fIndex) == true)
      looseCB = true;
  }

  const double isoTk = fWZTree->eleDr03TkSumPt->at(fIndex) / Pt();
  const double isoEcal = fWZTree->eleDr03EcalRecHitSumEt->at(fIndex) / Pt();
  const double isoHcal = fWZTree->eleDr03HcalTowerSumEt->at(fIndex) / Pt();

  if (looseCB && PassesPtCut() && PassesEtaCut() &&
      isoTk < 0.2 && isoEcal < 0.5 && isoHcal < 0.3)
    loose = true;

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

  bool mediumCB = false;
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
      mediumCB = true;
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
      mediumCB = true;
  }

  bool mediumVID = false;
  if (fWZTree->eleIDbit->at(fIndex)>>ELEMEDIUM_BIT&1)  mediumVID = true;
  if (mediumCB != mediumVID) {
    cout << "WARNING: VID different from Cut Based for Ele MEDIUM !!!" << endl;
    cout << "VID : " << mediumVID << ", Cut Based: " << mediumCB << endl;
  }

  if (mediumCB && PassesPtCut() && PassesEtaCut())  tight = true;

//  if (tight && !IsLoose())  cout << "WARNING: Electron passes Tight but fails Loose !!!" << endl;

  return tight;
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


double
Muon::ComputeLeptonScaleFactor()
{
  double sf = 0;

  if (IsTight()) {
    if (0 <= abs(Eta()) && abs(Eta()) < 1.5) {
      if (10 < Pt() && Pt() <= 20)          sf = 0.9843;
      else if (20 < Pt() && Pt() <= 30)     sf = 0.9891;
      else if (30 < Pt() && Pt() <= 40)     sf = 0.9884;
      else if (40 < Pt() && Pt() <= 50)     sf = 0.9884;
      else if (50 < Pt() && Pt() <= 100)    sf = 0.9861;
      else if (100 < Pt() && Pt() <= 1000)  sf = 0.9811;
    } else if (1.5 <= abs(Eta()) && abs(Eta()) < 2.5) {
      if (10 < Pt() && Pt() <= 20)          sf = 0.9887;
      else if (20 < Pt() && Pt() <= 30)     sf = 0.9853;
      else if (30 < Pt() && Pt() <= 40)     sf = 0.9883;
      else if (40 < Pt() && Pt() <= 50)     sf = 0.9893;
      else if (50 < Pt() && Pt() <= 100)    sf = 0.9892;
      else if (100 < Pt() && Pt() <= 1000)  sf = 0.9770;
    }
  } else if (IsLoose()) {
    cout << "WARNING: Lepton scale factor for LOOSE muon is computed!" << endl;
    if (0 <= abs(Eta()) && abs(Eta()) < 1.5) {
      if (10 < Pt() && Pt() <= 20)          sf = 1.0025;
      else if (20 < Pt() && Pt() <= 30)     sf = 0.9990;
      else if (30 < Pt() && Pt() <= 40)     sf = 0.9991;
      else if (40 < Pt() && Pt() <= 50)     sf = 0.9998;
      else if (50 < Pt() && Pt() <= 100)    sf = 0.9986;
      else if (100 < Pt() && Pt() <= 1000)  sf = 1.0004;
    } else if (1.5 <= abs(Eta()) && abs(Eta()) < 2.5) {
      if (10 < Pt() && Pt() <= 20)          sf = 0.9943;
      else if (20 < Pt() && Pt() <= 30)     sf = 0.9986;
      else if (30 < Pt() && Pt() <= 40)     sf = 0.9983;
      else if (40 < Pt() && Pt() <= 50)     sf = 0.9983;
      else if (50 < Pt() && Pt() <= 100)    sf = 0.9990;
      else if (100 < Pt() && Pt() <= 1000)  sf = 1.0059;
    }
  }

  return sf;
}


double
Muon::ComputePromptRate()
{
  double rate = 0.;

  if (abs(Eta()) <= MU_ETA_BARREL) {
    if (10 < Pt() && Pt() <= 20)        rate = 1.;
    else if (20 < Pt() && Pt() <= 30)   rate = 1.;
    else if (30 < Pt() && Pt() <= 200)  rate = 1.;
  } else if (abs(Eta()) > MU_ETA_BARREL && abs(Eta()) < MU_ETA_ENDCAP) {
    if (10 < Pt() && Pt() <= 20)        rate = 1.;
    else if (20 < Pt() && Pt() <= 30)   rate = 1.;
    else if (30 < Pt() && Pt() <= 200)  rate = 1.;
  }

  return rate;
}


double
Muon::ComputeFakeRate()
{
  double fake = 0.;

  if (abs(Eta()) <= MU_ETA_BARREL) {
    if (10 < Pt() && Pt() <= 20)        fake = 0.067;
    else if (20 < Pt() && Pt() <= 30)   fake = 0.085;
    else if (30 < Pt() && Pt() <= 200)  fake = 0.163;
  } else if (abs(Eta()) > MU_ETA_BARREL && abs(Eta()) < MU_ETA_ENDCAP) {
    if (10 < Pt() && Pt() <= 20)        fake = 0.103;
    else if (20 < Pt() && Pt() <= 30)   fake = 0.148;
    else if (30 < Pt() && Pt() <= 200)  fake = 0.229;
  }

  double rate = fake / (1. + fake);

  return rate;
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

  bool tightCB = false;
  if (fWZTree->muChi2NDF->at(fIndex) < 10. &&
      fWZTree->muMuonHits->at(fIndex) > 0 &&
      fWZTree->muStations->at(fIndex) > 1 &&
      abs(fWZTree->muD0->at(fIndex)) < 0.2 &&
      abs(fWZTree->muDz->at(fIndex)) < 0.5 &&
      fWZTree->muPixelHits->at(fIndex) > 0 &&
      fWZTree->muTrkLayers->at(fIndex) > 5 &&
      fWZTree->muType->at(fIndex)>>GLOBALMUON_BIT & 1 &&
      fWZTree->muType->at(fIndex)>>PFMUON_BIT & 1)
    tightCB = true;

  bool tightVID = false;
  if (fWZTree->muIsTightID->at(fIndex))  tightVID = true;
  if (tightCB != tightVID) {
    cout << "WARNING: VID different from Cut Based for Muon TIGHT !!!" << endl;
    cout << "VID : " << tightVID << ", Cut Based : " << tightCB << endl;
  }

  const double isoTk = fWZTree->muIsoTrk->at(fIndex) / Pt();

  if (tightCB && PassesPtCut() && PassesEtaCut() && isoTk < 0.4)  loose = true;

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

  bool tightCB = false;
  if (fWZTree->muChi2NDF->at(fIndex) < 10. &&
      fWZTree->muMuonHits->at(fIndex) > 0 &&
      fWZTree->muStations->at(fIndex) > 1 &&
      abs(fWZTree->muD0->at(fIndex)) < 0.2 &&
      abs(fWZTree->muDz->at(fIndex)) < 0.5 &&
      fWZTree->muPixelHits->at(fIndex) > 0 &&
      fWZTree->muTrkLayers->at(fIndex) > 5 &&
      fWZTree->muType->at(fIndex)>>GLOBALMUON_BIT & 1 &&
      fWZTree->muType->at(fIndex)>>PFMUON_BIT & 1)
    tightCB = true;

  bool tightVID = false;
  if (fWZTree->muIsTightID->at(fIndex))  tightVID = true;
  if (tightCB != tightVID) {
    cout << "WARNING: VID different from Cut Based for Muon TIGHT !!!" << endl;
    cout << "VID : " << tightVID << ", Cut Based : " << tightCB << endl;
  }

  if (tightCB && PassesPtCut() && PassesEtaCut() && fRelIso < MU_RELISO_TIGHT)  tight = true;

//  if (tight && !IsLoose())  cout << "WARNING: Muon passes Tight but fails Loose !!!" << endl;

  return tight;
}

