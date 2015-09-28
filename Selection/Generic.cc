#include "Generic.h"

#include <ios>


using namespace std;


Generic::Generic(WZEvent* wzEvent, TFile* outputFile)
{
  fWZEvent = wzEvent;
  fOutputRootFile = outputFile;
}


TH1D* Generic::bookTH1D(TString key, TString title, unsigned int nbins, double min, double max)
{
  TH1D* h = new TH1D(key, title, nbins, min, max);
  h->SetFillColor(kCyan);
  fHistos.push_back(h);
  return h;
}


TH2D* Generic::bookTH2D(TString key, TString title,
                                 unsigned int nbinsx, double xmin, double xmax,
                                 unsigned int nbinsy, double ymin, double ymax)
{
  TH2D* h = new TH2D(key, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
  h->SetOption("COLZ");
  fHistos.push_back(h);
  return h;
}


void
Generic::Init()
{

}


void
Generic::Analysis()
{

}


void
Generic::Finish()
{

}


void
Generic::WriteRootFile()
{
  cout << "Writing to ROOT.... \n";

  if (fOutputRootFile) {
    fOutputRootFile->cd();
    cout << "Objects to write : " << fHistos.size() << endl;
    for (unsigned int i = 0; i < fHistos.size(); i++) {
      if (typeid(*(fHistos.at(i))) == typeid(TH2D))  fHistos.at(i)->Draw();
      else if (typeid(*(fHistos.at(i))) == typeid(TH1D))  fHistos.at(i)->Draw();
//      cout << "Writing object : " << fHistos.at(i)->GetName() << endl;
      fHistos.at(i)->Write();
    }
  }
}

