#ifndef Particles_h
#define Particles_h

#include "TLorentzVector.h"
#include <iostream>


class WZEvent;


class Particle : public TLorentzVector
{

public:

  Particle(unsigned int index, double pt, double eta, double phi);
  Particle(unsigned int index, double pt, double eta, double phi, double en);

  int GetIndex() { return fIndex; }

  virtual bool IsLoose() = 0;
  virtual bool PassesPtCut() = 0;
  virtual bool PassesEtaCut() = 0;

  static void SetWZTree(WZEvent* wzTree);


protected:

  unsigned int fIndex;

  static WZEvent* fWZTree;

};


class Jet : public Particle
{

public:

  Jet(unsigned int index, double pt, double eta, double phi, double en);
  bool IsLoose();
  bool PassesPtCut();
  bool PassesEtaCut();

};


class Lepton : public Particle
{

public:

  Lepton(unsigned int index, double pt, double eta, double phi, double ch);

  virtual double ComputeLeptonScaleFactor() = 0;
  virtual double ComputeFakeRate() = 0;
  virtual double ComputePromptRate() = 0;

  virtual bool IsTight() = 0;

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

  bool IsLoose();
  bool IsTight();

  double ComputeLeptonScaleFactor();
  double ComputeFakeRate();
  double ComputePromptRate();

  bool PassesPtCut();
  bool PassesEtaCut();

};


class Muon : public Lepton
{

public:

  Muon(unsigned int ind, double pt, double eta, double phi, double ch);

  bool IsLoose();
  bool IsTight();

  double ComputeLeptonScaleFactor();
  double ComputeFakeRate();
  double ComputePromptRate();

  bool PassesPtCut();
  bool PassesEtaCut();

};


inline
bool HigherPt(const Particle* p1, const Particle* p2)
{
  return p1->Pt() > p2->Pt();
}


#endif

