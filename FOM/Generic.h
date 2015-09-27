#ifndef Generic_h
#define Generic_h

#include "Event.h"
#include "TH1D.h"
#include "TH2D.h"

#include <fstream>


class Generic
{

public:

  Generic(Event* event, TFile* outputFile = 0);

  virtual void Init();
  virtual void Analysis();

  TH1D* bookTH1D(TString key, TString title, unsigned int nbins, double min, double max);
  TH2D* bookTH2D(TString key, TString title,
                 unsigned int nbinsx, double xmin, double xmax,
                 unsigned int nbinsy, double ymin, double ymax);

  void WriteRootFile();

protected:

  Event* fEvent;
  
  vector<TObject*> fHistos;
  TFile* fOutputRootFile;

};

#endif

