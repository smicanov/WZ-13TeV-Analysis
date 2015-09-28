#ifndef Generic_h
#define Generic_h

#include "WZEvent.h"
#include "TH1D.h"
#include "TH2D.h"

#include <fstream>


class Generic
{

public:

  Generic(WZEvent* e, TFile* fout = 0);

  virtual void Init();
  virtual void Analysis();
  virtual void Finish();

  TH1D* bookTH1D(TString key, TString title, unsigned int nbins, double min, double max);
  TH2D* bookTH2D(TString key, TString title,
                 unsigned int nbinsx, double xmin, double xmax,
                 unsigned int nbinsy, double ymin, double ymax);

  void WriteRootFile();

protected:

  WZEvent* fWZEvent;
  
  std::vector<TObject*> fHistos;
  TFile* fOutputRootFile;

};

#endif
