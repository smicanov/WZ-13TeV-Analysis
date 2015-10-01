#ifndef Particles_h
#define Particles_h

#include "TLorentzVector.h"


class Event;


class Particle : public TLorentzVector
{

public:

  Particle(unsigned int index, double pt, double eta, double phi);
  Particle(unsigned int index, double pt, double eta, double phi, double en);

  int GetIndex() { return fIndex; }

  virtual bool IsLoose() = 0;
  virtual bool PassesPtCut() = 0;
  virtual bool PassesEtaCut() = 0;

  static void SetTree(Event* tree);


protected:

  int fIndex;

  static Event* fTree;

};


class Jet : public Particle
{

public:

  Jet(unsigned int index, double pt, double eta, double phi, double en, int partonId);
  bool IsLoose();
  bool PassesPtCut();
  bool PassesEtaCut();


protected:

  int fPartonId;

};


class Lepton : public Particle
{

public:

  Lepton(unsigned int index, double pt, double eta, double phi, double ch);

  virtual bool IsVeto() = 0;
  virtual bool IsMedium() = 0;
  virtual bool IsTight() = 0;
  virtual bool IsFOMLoose() = 0;
  virtual bool IsFOMTight() = 0;

  int GetPdgId() { return fPdgId; }
  double GetCharge() { return fCharge; }
  double GetRelIso() { return fRelIso; }


protected:

  int fPdgId;
  double fCharge;
  double fRelIso;

};


class Electron : public Lepton
{

public:

  Electron(unsigned int ind, double pt, double eta, double phi, double ch);

  bool IsVeto();
  bool IsLoose();
  bool IsMedium();
  bool IsTight();
  bool IsFOMLoose();
  bool IsFOMTight();

  bool PassesPtCut();
  bool PassesEtaCut();

};


class Muon : public Lepton
{

public:

  Muon(unsigned int ind, double pt, double eta, double phi, double ch);

  bool IsVeto();
  bool IsLoose();
  bool IsMedium();
  bool IsTight();
  bool IsFOMLoose();
  bool IsFOMTight();

  bool PassesPtCut();
  bool PassesEtaCut();

};


inline
bool HigherPt(const Particle* particle1, const Particle* particle2)
{
  return particle1->Pt() > particle2->Pt();
}


inline
bool HigherPtJets(const Jet j1, const Jet j2)
{
  return j1.Pt() > j2.Pt();
}


#endif

