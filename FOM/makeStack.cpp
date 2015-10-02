#include "XSections.h"
#include "MyStyle.h"

#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TLatex.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TMath.h"

#include <iomanip>
#include <cmath>
#include <string>
#include <boost/lexical_cast.hpp>


using namespace std;


int
MakeStack(TFile* f_WZ, TFile* f_ZZ4L, TFile* f_ZZ2L2Q, TFile* f_WW, TFile* f_DYM50,
          TFile* f_TT, TFile* f_DYM10, TFile* f_WJets,
          TString histKey, unsigned int channel, TString xAxisTitle = "", double binWidth = 1.)
{
  TH1D* h_WZ = (TH1D*)(f_WZ->Get(histKey))->Clone(histKey + "_WZ");
  TH1D* h_ZZ4L = (TH1D*)(f_ZZ4L->Get(histKey))->Clone(histKey + "_ZZ4L");
  TH1D* h_ZZ2L2Q = (TH1D*)(f_ZZ2L2Q->Get(histKey))->Clone(histKey + "_ZZ2L2Q");
  TH1D* h_WW = (TH1D*)(f_WW->Get(histKey))->Clone(histKey + "_WW");
  TH1D* h_DYM50 = (TH1D*)(f_DYM50->Get(histKey))->Clone(histKey + "_DYM50");
  TH1D* h_DYM10 = (TH1D*)(f_DYM10->Get(histKey))->Clone(histKey + "_DYM10");
  TH1D* h_TT = (TH1D*)(f_TT->Get(histKey))->Clone(histKey + "_TT");
  TH1D* h_WJets = (TH1D*)(f_WJets->Get(histKey))->Clone(histKey + "_WJets");

  const double scale_WZ = EXPCMS_LUMINOSITY / (NG_WZTo3LNu / XS_WZTo3LNu);
  const double scale_ZZ4L = EXPCMS_LUMINOSITY /  (NG_ZZTo4L / XS_ZZTo4L);
  const double scale_ZZ2L2Q = EXPCMS_LUMINOSITY /  (NG_ZZTo2L2Q / XS_ZZTo2L2Q);
  const double scale_WW = EXPCMS_LUMINOSITY /  (NG_WWTo2L2Nu / XS_WWTo2L2Nu);
  const double scale_DYM50 = EXPCMS_LUMINOSITY /  (NG_DYJetsToLL_M50 / XS_DYJetsToLL_M50);
  const double scale_DYM10 = EXPCMS_LUMINOSITY /  (NG_DYJetsToLL_M10to50 / XS_DYJetsToLL_M10to50);
  const double scale_TT = EXPCMS_LUMINOSITY /  (NG_TTJets / XS_TTJets);
  const double scale_WJets = EXPCMS_LUMINOSITY /  (NG_WJetsToLNu / XS_WJetsToLNu);

  h_WZ->Scale(scale_WZ);
  h_ZZ4L->Scale(scale_ZZ4L);
  h_ZZ2L2Q->Scale(scale_ZZ2L2Q);
  h_WW->Scale(scale_WW);
  h_DYM50->Scale(scale_DYM50);
  h_DYM10->Scale(scale_DYM10);
  h_TT->Scale(scale_TT);
  h_WJets->Scale(scale_WJets);

  THStack* stack = new THStack(histKey, histKey);
  stack->Add(h_WW);
  stack->Add(h_ZZ2L2Q);
  stack->Add(h_TT);
  stack->Add(h_ZZ4L);
  stack->Add(h_DYM50);
  stack->Add(h_DYM10);
  stack->Add(h_WZ);
  stack->Add(h_WJets);

  h_WZ->SetFillColor(kOrange-2);
  h_ZZ4L->SetFillColor(kRed+1);
  h_ZZ2L2Q->SetFillColor(kRed+2);
  h_WW->SetFillColor(kGreen);
  h_DYM50->SetFillColor(kGray+1);
  h_DYM10->SetFillColor(kGray+3);
  h_TT->SetFillColor(kAzure);
  h_WJets->SetFillColor(kMagenta-7);

  stack->SetMinimum(0.0);
  stack->SetMaximum(stack->GetMaximum() * 1.33);

  stack->Draw();
  stack->GetXaxis()->SetTitle(xAxisTitle);
  if (binWidth >= 0.1)  stack->GetYaxis()->SetTitle(Form("Events / %1.1f GeV", binWidth));
  else if (binWidth > 0. && binWidth < 0.1)
    stack->GetYaxis()->SetTitle(Form("Events / %1.3f", binWidth));
  else  stack->GetYaxis()->SetTitle("Events");

  stack->GetXaxis()->SetLabelFont(132);
  stack->GetYaxis()->SetLabelFont(132);
  stack->GetXaxis()->SetLabelOffset(0.007);
  stack->GetYaxis()->SetLabelOffset(0.007);
  stack->GetXaxis()->SetLabelSize(0.03);
  stack->GetYaxis()->SetLabelSize(0.03);
  stack->GetXaxis()->SetTitleFont(132);
  stack->GetYaxis()->SetTitleFont(132);
  stack->GetXaxis()->SetTitleSize(0.045);
  stack->GetYaxis()->SetTitleSize(0.04);
  stack->GetXaxis()->SetTitleOffset(1.);
  stack->GetYaxis()->SetTitleOffset(1.5);
  stack->GetXaxis()->SetNdivisions(510);
  stack->GetYaxis()->SetNdivisions(510);

  TLegend* leg = new TLegend(0.7, 0.65, 0.9, 0.85);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetShadowColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);

  leg->AddEntry(h_WJets, "W + jets", "f");
  leg->AddEntry(h_WZ, "WZ", "f");
  leg->AddEntry(h_DYM50, "Drell Yan (Z + jets)", "f");
  leg->AddEntry(h_ZZ4L, "ZZ #rightarrow 4l", "f");
  leg->AddEntry(h_TT, "TTJets", "f");
  leg->AddEntry(h_ZZ2L2Q, "ZZ #rightarrow 2l2q", "f");
  leg->AddEntry(h_WW, "WW", "f");
  leg->AddEntry(h_DYM10, "Z + jets (low mass)", "f");

  leg->Draw();

  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextAlign(12);
  latexLabel.SetTextSize(0.04);
  latexLabel.DrawLatex(0.19, 0.9, Form("#font[132]{#intL dt = %1.1f fb^{-1}}", 1.));

  latexLabel.SetTextSize(0.05);
  if (channel == 0) latexLabel.DrawLatex(0.44, 0.95, "#font[132]{Channel undefined}");
  if (channel == 1) latexLabel.DrawLatex(0.4, 0.95, "#font[132]{Channel W ele - Fake ele}");
  if (channel == 2) latexLabel.DrawLatex(0.4, 0.95, "#font[132]{Channel W ele - Fake #mu}");
  if (channel == 3) latexLabel.DrawLatex(0.4, 0.95, "#font[132]{Channel W #mu - Fake ele}");
  if (channel == 4) latexLabel.DrawLatex(0.4, 0.95, "#font[132]{Channel W #mu - Fake #mu}");
  if (channel == 5) latexLabel.DrawLatex(0.46, 0.95, "#font[132]{All Channels}");

  return 1;
}


int
MakeStackDeltaR(TFile* f_WZ, TFile* f_ZZ4L, TFile* f_ZZ2L2Q, TFile* f_WW, TFile* f_DYM50,
                TFile* f_TT, TFile* f_DYM10, TFile* f_WJets,
                TString histKey, TString xAxisTitle, double max)
{
  TH1D* h_WZ = (TH1D*)(f_WZ->Get(histKey))->Clone(histKey + "_WZ");
  TH1D* h_ZZ4L = (TH1D*)(f_ZZ4L->Get(histKey))->Clone(histKey + "_ZZ4L");
  TH1D* h_ZZ2L2Q = (TH1D*)(f_ZZ2L2Q->Get(histKey))->Clone(histKey + "_ZZ2L2Q");
  TH1D* h_WW = (TH1D*)(f_WW->Get(histKey))->Clone(histKey + "_WW");
  TH1D* h_DYM50 = (TH1D*)(f_DYM50->Get(histKey))->Clone(histKey + "_DYM50");
  TH1D* h_DYM10 = (TH1D*)(f_DYM10->Get(histKey))->Clone(histKey + "_DYM10");
  TH1D* h_TT = (TH1D*)(f_TT->Get(histKey))->Clone(histKey + "_TT");
  TH1D* h_WJets = (TH1D*)(f_WJets->Get(histKey))->Clone(histKey + "_WJets");

  const double scale_WZ = EXPCMS_LUMINOSITY / (NG_WZTo3LNu / XS_WZTo3LNu);
  const double scale_ZZ4L = EXPCMS_LUMINOSITY /  (NG_ZZTo4L / XS_ZZTo4L);
  const double scale_ZZ2L2Q = EXPCMS_LUMINOSITY /  (NG_ZZTo2L2Q / XS_ZZTo2L2Q);
  const double scale_WW = EXPCMS_LUMINOSITY /  (NG_WWTo2L2Nu / XS_WWTo2L2Nu);
  const double scale_DYM50 = EXPCMS_LUMINOSITY /  (NG_DYJetsToLL_M50 / XS_DYJetsToLL_M50);
  const double scale_DYM10 = EXPCMS_LUMINOSITY /  (NG_DYJetsToLL_M10to50 / XS_DYJetsToLL_M10to50);
  const double scale_TT = EXPCMS_LUMINOSITY /  (NG_TTJets / XS_TTJets);
  const double scale_WJets = EXPCMS_LUMINOSITY /  (NG_WJetsToLNu / XS_WJetsToLNu);

  h_WZ->Scale(scale_WZ);
  h_ZZ4L->Scale(scale_ZZ4L);
  h_ZZ2L2Q->Scale(scale_ZZ2L2Q);
  h_WW->Scale(scale_WW);
  h_DYM50->Scale(scale_DYM50);
  h_DYM10->Scale(scale_DYM10);
  h_TT->Scale(scale_TT);
  h_WJets->Scale(scale_WJets);

  THStack* stack = new THStack(histKey, histKey);
  stack->Add(h_WW);
  stack->Add(h_ZZ2L2Q);
  stack->Add(h_TT);
  stack->Add(h_ZZ4L);
  stack->Add(h_DYM50);
  stack->Add(h_DYM10);
  stack->Add(h_WZ);
  stack->Add(h_WJets);

  h_WZ->SetFillColor(kOrange-2);
  h_ZZ4L->SetFillColor(kRed+1);
  h_ZZ2L2Q->SetFillColor(kRed+2);
  h_WW->SetFillColor(kGreen);
  h_DYM50->SetFillColor(kGray+1);
  h_DYM10->SetFillColor(kGray+3);
  h_TT->SetFillColor(kAzure);
  h_WJets->SetFillColor(kMagenta-7);

  stack->SetMinimum(0.0);
  stack->SetMaximum(max);

  stack->Draw();
  stack->GetXaxis()->SetTitle(xAxisTitle);
  stack->GetYaxis()->SetTitle("Jets / 0.02");

  stack->GetXaxis()->SetLabelFont(132);
  stack->GetYaxis()->SetLabelFont(132);
  stack->GetXaxis()->SetLabelOffset(0.007);
  stack->GetYaxis()->SetLabelOffset(0.007);
  stack->GetXaxis()->SetLabelSize(0.03);
  stack->GetYaxis()->SetLabelSize(0.03);
  stack->GetXaxis()->SetTitleFont(132);
  stack->GetYaxis()->SetTitleFont(132);
  stack->GetXaxis()->SetTitleSize(0.045);
  stack->GetYaxis()->SetTitleSize(0.04);
  stack->GetXaxis()->SetTitleOffset(1.);
  stack->GetYaxis()->SetTitleOffset(1.5);
  stack->GetXaxis()->SetNdivisions(510);
  stack->GetYaxis()->SetNdivisions(510);

  TLegend* leg = new TLegend(0.7, 0.65, 0.9, 0.85);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetShadowColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);

  leg->AddEntry(h_WJets, "W + jets", "f");
  leg->AddEntry(h_WZ, "WZ", "f");
  leg->AddEntry(h_DYM50, "Drell Yan (Z + jets)", "f");
  leg->AddEntry(h_ZZ4L, "ZZ #rightarrow 4l", "f");
  leg->AddEntry(h_TT, "TTJets", "f");
  leg->AddEntry(h_ZZ2L2Q, "ZZ #rightarrow 2l2q", "f");
  leg->AddEntry(h_WW, "WW", "f");
  leg->AddEntry(h_DYM10, "Z + jets (low mass)", "f");

  leg->Draw();

  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextAlign(12);
  latexLabel.SetTextSize(0.04);
  latexLabel.DrawLatex(0.19, 0.9, Form("#font[132]{#intL dt = %1.1f fb^{-1}}", 1.));

  return 1;
}


pair<double, double>
SignalAboveCutBin(TFile* f_WZ, TFile* f_ZZ4L, TFile* f_ZZ2L2Q, TFile* f_WW, TFile* f_DYM50,
                  TFile* f_TT, TFile* f_DYM10, TFile* f_WJets,
                  TString histKey, TFile* f_Signal, double scaleSignal, unsigned int cutBin)
{
  TH1D* h_WZ = (TH1D*)(f_WZ->Get(histKey))->Clone(histKey + "_WZ");
  TH1D* h_ZZ4L = (TH1D*)(f_ZZ4L->Get(histKey))->Clone(histKey + "_ZZ4L");
  TH1D* h_ZZ2L2Q = (TH1D*)(f_ZZ2L2Q->Get(histKey))->Clone(histKey + "_ZZ2L2Q");
  TH1D* h_WW = (TH1D*)(f_WW->Get(histKey))->Clone(histKey + "_WW");
  TH1D* h_DYM50 = (TH1D*)(f_DYM50->Get(histKey))->Clone(histKey + "_DYM50");
  TH1D* h_DYM10 = (TH1D*)(f_DYM10->Get(histKey))->Clone(histKey + "_DYM10");
  TH1D* h_TT = (TH1D*)(f_TT->Get(histKey))->Clone(histKey + "_TT");
  TH1D* h_WJets = (TH1D*)(f_WJets->Get(histKey))->Clone(histKey + "_WJets");

  const double scale_WZ = EXPCMS_LUMINOSITY / (NG_WZTo3LNu / XS_WZTo3LNu);
  const double scale_ZZ4L = EXPCMS_LUMINOSITY /  (NG_ZZTo4L / XS_ZZTo4L);
  const double scale_ZZ2L2Q = EXPCMS_LUMINOSITY /  (NG_ZZTo2L2Q / XS_ZZTo2L2Q);
  const double scale_WW = EXPCMS_LUMINOSITY /  (NG_WWTo2L2Nu / XS_WWTo2L2Nu);
  const double scale_DYM50 = EXPCMS_LUMINOSITY /  (NG_DYJetsToLL_M50 / XS_DYJetsToLL_M50);
  const double scale_DYM10 = EXPCMS_LUMINOSITY /  (NG_DYJetsToLL_M10to50 / XS_DYJetsToLL_M10to50);
  const double scale_TT = EXPCMS_LUMINOSITY /  (NG_TTJets / XS_TTJets);
  const double scale_WJets = EXPCMS_LUMINOSITY /  (NG_WJetsToLNu / XS_WJetsToLNu);

  h_WZ->Scale(scale_WZ);
  h_ZZ4L->Scale(scale_ZZ4L);
  h_ZZ2L2Q->Scale(scale_ZZ2L2Q);
  h_WW->Scale(scale_WW);
  h_DYM50->Scale(scale_DYM50);
  h_DYM10->Scale(scale_DYM10);
  h_TT->Scale(scale_TT);
  h_WJets->Scale(scale_WJets);

  TH1D* h_All = (TH1D*) h_WJets->Clone("h_All");
  h_All->Add(h_WZ);
  h_All->Add(h_WW);
  h_All->Add(h_ZZ2L2Q);
  h_All->Add(h_TT);
  h_All->Add(h_ZZ4L);
  h_All->Add(h_DYM50);
  h_All->Add(h_DYM10);

  TH1D* h_Signal = (TH1D*) (f_Signal->Get(histKey))->Clone("h_Signal");
  h_Signal->Scale(scaleSignal);

  double signal = h_Signal->Integral(cutBin, 300);
  double total = h_All->Integral(cutBin, 300);
  double signalFraction = 0;
  double significance = 0;
  if (total) {
    signalFraction = signal / total;
    significance = signal / sqrt(total);
  } else {
    signalFraction = 0;
    significance = 0;
  }

  pair<double, double> result = make_pair(signalFraction, significance);
  return result;
}


pair<double, double>
TotalSignalAboveCutBin(TFile* f_WZ, TFile* f_ZZ4L, TFile* f_ZZ2L2Q, TFile* f_WW, TFile* f_DYM50,
                       TFile* f_TT, TFile* f_DYM10, TFile* f_WJets,
                       TString histKey, unsigned int cutBin)
{
  TH1D* h_WZ = (TH1D*)(f_WZ->Get(histKey))->Clone(histKey + "_WZ");
  TH1D* h_ZZ4L = (TH1D*)(f_ZZ4L->Get(histKey))->Clone(histKey + "_ZZ4L");
  TH1D* h_ZZ2L2Q = (TH1D*)(f_ZZ2L2Q->Get(histKey))->Clone(histKey + "_ZZ2L2Q");
  TH1D* h_WW = (TH1D*)(f_WW->Get(histKey))->Clone(histKey + "_WW");
  TH1D* h_DYM50 = (TH1D*)(f_DYM50->Get(histKey))->Clone(histKey + "_DYM50");
  TH1D* h_DYM10 = (TH1D*)(f_DYM10->Get(histKey))->Clone(histKey + "_DYM10");
  TH1D* h_TT = (TH1D*)(f_TT->Get(histKey))->Clone(histKey + "_TT");
  TH1D* h_WJets = (TH1D*)(f_WJets->Get(histKey))->Clone(histKey + "_WJets");

  const double scale_WZ = EXPCMS_LUMINOSITY / (NG_WZTo3LNu / XS_WZTo3LNu);
  const double scale_ZZ4L = EXPCMS_LUMINOSITY /  (NG_ZZTo4L / XS_ZZTo4L);
  const double scale_ZZ2L2Q = EXPCMS_LUMINOSITY /  (NG_ZZTo2L2Q / XS_ZZTo2L2Q);
  const double scale_WW = EXPCMS_LUMINOSITY /  (NG_WWTo2L2Nu / XS_WWTo2L2Nu);
  const double scale_DYM50 = EXPCMS_LUMINOSITY /  (NG_DYJetsToLL_M50 / XS_DYJetsToLL_M50);
  const double scale_DYM10 = EXPCMS_LUMINOSITY /  (NG_DYJetsToLL_M10to50 / XS_DYJetsToLL_M10to50);
  const double scale_TT = EXPCMS_LUMINOSITY /  (NG_TTJets / XS_TTJets);
  const double scale_WJets = EXPCMS_LUMINOSITY /  (NG_WJetsToLNu / XS_WJetsToLNu);

  h_WZ->Scale(scale_WZ);
  h_ZZ4L->Scale(scale_ZZ4L);
  h_ZZ2L2Q->Scale(scale_ZZ2L2Q);
  h_WW->Scale(scale_WW);
  h_DYM50->Scale(scale_DYM50);
  h_DYM10->Scale(scale_DYM10);
  h_TT->Scale(scale_TT);
  h_WJets->Scale(scale_WJets);

  TH1D* h_All = (TH1D*) h_WJets->Clone("h_All");
  h_All->Add(h_WZ);
  h_All->Add(h_WW);
  h_All->Add(h_ZZ2L2Q);
  h_All->Add(h_TT);
  h_All->Add(h_ZZ4L);
  h_All->Add(h_DYM50);
  h_All->Add(h_DYM10);

  double signal = h_All->Integral(cutBin, 300);
  double total = h_All->Integral(0, 300);
  double signalFraction = 0.;
  (total != 0.0) ? signalFraction = signal / total : signalFraction = 0.;

  double signal_WJets = h_WJets->Integral(cutBin, 300);
  double totalSignal_WJets = h_WJets->Integral(0, 300);
  double sF_WJets = 0.;
  (totalSignal_WJets != 0.0) ? sF_WJets = signal_WJets / totalSignal_WJets : sF_WJets = 0.;

  pair<double, double> sF = make_pair(signalFraction, sF_WJets);
  return sF;
}


int
MakeSignificanceGraphs(TFile* f_WZ, TFile* f_ZZ4L, TFile* f_ZZ2L2Q, TFile* f_WW, TFile* f_DYM50,
                       TFile* f_TT, TFile* f_DYM10, TFile* f_WJets,
                       TString histKey, unsigned int channel, TCanvas* canvas, TString xAxisTitle,
                       double start, double end, double binWidth)
{
  const double scale_DYM50 = EXPCMS_LUMINOSITY /  (NG_DYJetsToLL_M50 / XS_DYJetsToLL_M50);
  const double scale_DYM10 = EXPCMS_LUMINOSITY /  (NG_DYJetsToLL_M50 / XS_DYJetsToLL_M50);
  const double scale_TT = EXPCMS_LUMINOSITY /  (NG_TTJets / XS_TTJets);
  const double scale_WJets = EXPCMS_LUMINOSITY /  (NG_WJetsToLNu / XS_WJetsToLNu);

  vector<double> xValues;
  vector<double> totalSignalFractions;
  vector<double> signalFractions_WJets, signalFractions_DYM50, signalFractions_DYM10, signalFractions_TT;
  vector<double> significances_WJets, significances_DYM50, significances_DYM10, significances_TT;
  const unsigned int endBin = (unsigned int)((end - start) / binWidth);
  for (unsigned int i = 0; i < endBin+1; i++) {
    const double xValue = start + i * binWidth;
    xValues.push_back(xValue);
    const double signalFraction_WJets =
      SignalAboveCutBin(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT, f_DYM10, f_WJets,
                        histKey, f_WJets, scale_WJets, i).first;
    signalFractions_WJets.push_back(signalFraction_WJets);
    const double significance_WJets =
      SignalAboveCutBin(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT, f_DYM10, f_WJets,
                        histKey, f_WJets, scale_WJets, i).second;
    significances_WJets.push_back(significance_WJets);
    const double signalFraction_DYM50 =
      SignalAboveCutBin(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT, f_DYM10, f_WJets,
                        histKey, f_DYM50, scale_DYM50, i).first;
    signalFractions_DYM50.push_back(signalFraction_DYM50);
    const double significance_DYM50 =
      SignalAboveCutBin(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT, f_DYM10, f_WJets,
                        histKey, f_DYM50, scale_DYM50, i).second;
    significances_DYM50.push_back(significance_DYM50);
    const double signalFraction_DYM10 =
      SignalAboveCutBin(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT, f_DYM10, f_WJets,
                        histKey, f_DYM10, scale_DYM10, i).first;
    signalFractions_DYM10.push_back(signalFraction_DYM10);
    const double significance_DYM10 =
      SignalAboveCutBin(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT, f_DYM10, f_WJets,
                        histKey, f_DYM10, scale_DYM10, i).second;
    significances_DYM10.push_back(significance_DYM10);
    const double signalFraction_TT =
      SignalAboveCutBin(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT, f_DYM10, f_WJets,
                        histKey, f_TT, scale_TT, i).first;
    signalFractions_TT.push_back(signalFraction_TT);
    const double significance_TT =
      SignalAboveCutBin(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT, f_DYM10, f_WJets,
                        histKey, f_TT, scale_TT, i).second;
    significances_TT.push_back(significance_TT);
    const double signalFraction_Total =
      TotalSignalAboveCutBin(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT, f_DYM10, f_WJets,
                             histKey, i).second;
    totalSignalFractions.push_back(signalFraction_Total);
  }

  canvas->Divide(1, 2);
  
  canvas->cd(1);
  TGraph* g_sF_WJets = new TGraph(xValues.size(), &xValues.at(0), &signalFractions_WJets.at(0));
  TGraph* g_sF_DYM50 = new TGraph(xValues.size(), &xValues.at(0), &signalFractions_DYM50.at(0));
  TGraph* g_sF_DYM10 = new TGraph(xValues.size(), &xValues.at(0), &signalFractions_DYM10.at(0));
  TGraph* g_sF_TT = new TGraph(xValues.size(), &xValues.at(0), &signalFractions_TT.at(0));
  TGraph* g_sF_Total = new TGraph(xValues.size(), &xValues.at(0), &totalSignalFractions.at(0));

  g_sF_WJets->SetMarkerColor(kMagenta-7);
  g_sF_WJets->SetMarkerSize(0.6);
  g_sF_WJets->SetMaximum(1.44);
  g_sF_WJets->SetMinimum(0.);
  g_sF_WJets->Draw("AP");
  g_sF_WJets->GetYaxis()->SetTitle("S/(S+B)");
  g_sF_WJets->GetXaxis()->SetTitle(xAxisTitle);

  g_sF_WJets->GetXaxis()->SetLabelFont(132);
  g_sF_WJets->GetYaxis()->SetLabelFont(132);
  g_sF_WJets->GetXaxis()->SetLabelOffset(0.007);
  g_sF_WJets->GetYaxis()->SetLabelOffset(0.007);
  g_sF_WJets->GetXaxis()->SetLabelSize(0.03);
  g_sF_WJets->GetYaxis()->SetLabelSize(0.03);
  g_sF_WJets->GetXaxis()->SetTitleFont(132);
  g_sF_WJets->GetYaxis()->SetTitleFont(132);
  g_sF_WJets->GetXaxis()->SetTitleSize(0.045);
  g_sF_WJets->GetYaxis()->SetTitleSize(0.04);
  g_sF_WJets->GetXaxis()->SetTitleOffset(1.);
  g_sF_WJets->GetYaxis()->SetTitleOffset(1.5);
  g_sF_WJets->GetXaxis()->SetNdivisions(510);
  g_sF_WJets->GetYaxis()->SetNdivisions(510);

  g_sF_DYM50->SetMarkerColor(kGray+1);
  g_sF_DYM50->SetMarkerSize(0.6);
  g_sF_DYM50->SetMaximum(1.44);
  g_sF_DYM50->SetMinimum(0.);
  g_sF_DYM50->Draw("P");

  g_sF_DYM10->SetMarkerColor(kGray+3);
  g_sF_DYM10->SetMarkerSize(0.6);
  g_sF_DYM10->SetMaximum(1.44);
  g_sF_DYM10->SetMinimum(0.);
  g_sF_DYM10->Draw("P");

  g_sF_TT->SetMarkerColor(kAzure);
  g_sF_TT->SetMarkerSize(0.6);
  g_sF_TT->SetMaximum(1.44);
  g_sF_TT->SetMinimum(0.);
  g_sF_TT->Draw("P");

  g_sF_Total->SetMarkerColor(kBlack);
  g_sF_Total->SetMarkerSize(0.4);
  g_sF_Total->SetMaximum(1.44);
  g_sF_Total->SetMinimum(0.);
  g_sF_Total->Draw("*");

  TLegend* leg_sF = new TLegend(0.18, 0.75, 0.38, 0.95);
  leg_sF->SetFillColor(0);
  leg_sF->SetFillStyle(0);
  leg_sF->SetShadowColor(0);
  leg_sF->SetBorderSize(0);
  leg_sF->SetTextSize(0.03);

  leg_sF->AddEntry(g_sF_WJets, "W + jets", "p");
  leg_sF->AddEntry(g_sF_DYM50, "Drell Yan (Z + jets)", "p");
  leg_sF->AddEntry(g_sF_DYM10, "Z + jets (low mass)", "p");
  leg_sF->AddEntry(g_sF_TT, "TTJets", "p");
  leg_sF->AddEntry(g_sF_Total, "Total W+jets Signal", "p");

  leg_sF->Draw();

  TLatex latexLabel_sF;
  latexLabel_sF.SetNDC();
  latexLabel_sF.SetTextAlign(12);
  latexLabel_sF.SetTextSize(0.05);
  latexLabel_sF.DrawLatex(0.8, 0.9, "Signal Fraction");

  if (channel == 0) latexLabel_sF.DrawLatex(0.46, 0.95, "#font[132]{Channel undefined}");
  if (channel == 1) latexLabel_sF.DrawLatex(0.46, 0.95, "#font[132]{Channel W ele - Fake ele}");
  if (channel == 2) latexLabel_sF.DrawLatex(0.46, 0.95, "#font[132]{Channel W ele - Fake #mu}");
  if (channel == 3) latexLabel_sF.DrawLatex(0.46, 0.95, "#font[132]{Channel W #mu - Fake ele}");
  if (channel == 4) latexLabel_sF.DrawLatex(0.46, 0.95, "#font[132]{Channel W #mu - Fake #mu}");
  if (channel == 5) latexLabel_sF.DrawLatex(0.46, 0.95, "#font[132]{All Channels}");

  canvas->cd(2);
  TGraph* g_significance_WJets = new TGraph(xValues.size(), &xValues.at(0), &significances_WJets.at(0));
  TGraph* g_significance_DYM50 = new TGraph(xValues.size(), &xValues.at(0), &significances_DYM50.at(0));
  TGraph* g_significance_DYM10 = new TGraph(xValues.size(), &xValues.at(0), &significances_DYM10.at(0));
  TGraph* g_significance_TT = new TGraph(xValues.size(), &xValues.at(0), &significances_TT.at(0));

  vector<double> max_significances =
    { TMath::MaxElement(g_significance_WJets->GetN(), g_significance_WJets->GetY()),
      TMath::MaxElement(g_significance_DYM50->GetN(), g_significance_DYM50->GetY()),
      TMath::MaxElement(g_significance_DYM10->GetN(), g_significance_DYM10->GetY()),
      TMath::MaxElement(g_significance_TT->GetN(), g_significance_TT->GetY()) };
  double max_significance = *(max_element(max_significances.begin(), max_significances.end()));

  g_significance_WJets->SetMarkerColor(kMagenta-7);
  g_significance_WJets->SetMarkerSize(0.6);
  g_significance_WJets->SetMaximum(max_significance * 1.44);
  g_significance_WJets->SetMinimum(0.);
  g_significance_WJets->Draw("AP");
  g_significance_WJets->GetYaxis()->SetTitle("S/#sqrt{S+B}");
  g_significance_WJets->GetXaxis()->SetTitle(xAxisTitle);

  g_significance_WJets->GetXaxis()->SetLabelFont(132);
  g_significance_WJets->GetYaxis()->SetLabelFont(132);
  g_significance_WJets->GetXaxis()->SetLabelOffset(0.007);
  g_significance_WJets->GetYaxis()->SetLabelOffset(0.007);
  g_significance_WJets->GetXaxis()->SetLabelSize(0.03);
  g_significance_WJets->GetYaxis()->SetLabelSize(0.03);
  g_significance_WJets->GetXaxis()->SetTitleFont(132);
  g_significance_WJets->GetYaxis()->SetTitleFont(132);
  g_significance_WJets->GetXaxis()->SetTitleSize(0.045);
  g_significance_WJets->GetYaxis()->SetTitleSize(0.04);
  g_significance_WJets->GetXaxis()->SetTitleOffset(1.);
  g_significance_WJets->GetYaxis()->SetTitleOffset(1.5);
  g_significance_WJets->GetXaxis()->SetNdivisions(510);
  g_significance_WJets->GetYaxis()->SetNdivisions(510);

  g_significance_DYM50->SetMarkerColor(kGray+1);
  g_significance_DYM50->SetMarkerSize(0.6);
  g_significance_DYM50->SetMaximum(max_significance * 1.4);
  g_significance_DYM50->SetMinimum(0.);
  g_significance_DYM50->Draw("P");

  g_significance_DYM10->SetMarkerColor(kGray+3);
  g_significance_DYM10->SetMarkerSize(0.6);
  g_significance_DYM10->SetMaximum(max_significance * 1.44);
  g_significance_DYM10->SetMinimum(0.);
  g_significance_DYM10->Draw("P");

  g_significance_TT->SetMarkerColor(kAzure);
  g_significance_TT->SetMarkerSize(0.6);
  g_significance_TT->SetMaximum(max_significance * 1.44);
  g_significance_TT->SetMinimum(0.);
  g_significance_TT->Draw("P");

  TLegend* leg_significance = new TLegend(0.18, 0.75, 0.38, 0.95);
  leg_significance->SetFillColor(0);
  leg_significance->SetFillStyle(0);
  leg_significance->SetShadowColor(0);
  leg_significance->SetBorderSize(0);
  leg_significance->SetTextSize(0.03);

  leg_significance->AddEntry(g_sF_WJets, "W + jets", "p");
  leg_significance->AddEntry(g_sF_DYM50, "Drell Yan (Z + jets)", "p");
  leg_significance->AddEntry(g_sF_DYM10, "Z + jets (low mass)", "p");
  leg_significance->AddEntry(g_sF_TT, "TTJets", "p");

  leg_significance->Draw();

  TLatex latexLabel_significance;
  latexLabel_significance.SetNDC();
  latexLabel_significance.SetTextAlign(12);
  latexLabel_significance.SetTextSize(0.05);
  latexLabel_significance.DrawLatex(0.82, 0.9, "Significance");

  if (channel == 0) latexLabel_significance.DrawLatex(0.46, 0.95, "#font[132]{Channel undefined}");
  if (channel == 1) latexLabel_significance.DrawLatex(0.46, 0.95, "#font[132]{Channel W ele - Fake ele}");
  if (channel == 2) latexLabel_significance.DrawLatex(0.46, 0.95, "#font[132]{Channel W ele - Fake #mu}");
  if (channel == 3) latexLabel_significance.DrawLatex(0.46, 0.95, "#font[132]{Channel W #mu - Fake ele}");
  if (channel == 4) latexLabel_significance.DrawLatex(0.46, 0.95, "#font[132]{Channel W #mu - Fake #mu}");
  if (channel == 5) latexLabel_significance.DrawLatex(0.46, 0.95, "#font[132]{All Channels}");

  return 1;
}


int
main()
{
  const MyStyle rootStyle(600);

  TFile* f_WZ   = TFile::Open("/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/fom/GoodJetPt35GeV_LeadJetPt40GeV/NoRelIso/SSSelection/WZ.root");
  TFile* f_ZZ4L   = TFile::Open("/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/fom/GoodJetPt35GeV_LeadJetPt40GeV/NoRelIso/SSSelection/ZZ4L.root");
  TFile* f_ZZ2L2Q   = TFile::Open("/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/fom/GoodJetPt35GeV_LeadJetPt40GeV/NoRelIso/SSSelection/ZZ2L2Q.root");
  TFile* f_WW   = TFile::Open("/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/fom/GoodJetPt35GeV_LeadJetPt40GeV/NoRelIso/SSSelection/WW.root");
  TFile* f_DYM50   = TFile::Open("/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/fom/GoodJetPt35GeV_LeadJetPt40GeV/NoRelIso/SSSelection/DYM50.root");
  TFile* f_TT   = TFile::Open("/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/fom/GoodJetPt35GeV_LeadJetPt40GeV/NoRelIso/SSSelection/TT.root");
  TFile* f_DYM10   = TFile::Open("/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/fom/GoodJetPt35GeV_LeadJetPt40GeV/NoRelIso/SSSelection/DYM10to50.root");
  TFile* f_WJets   = TFile::Open("/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/fom/GoodJetPt35GeV_LeadJetPt40GeV/NoRelIso/SSSelection/WJets.root");

  unsigned int n = 6;

  vector<string> histoName =
    { "WPt", "WEta", "WPhi", "WRelIso", "FakePt", "FakeEta", "FakePhi", "FakeRelIso",
      "DeltaPhiWFake", "DeltaRWFake", "MET", "METPhi",
      "DeltaPhiWMET", "DeltaRWMET", "DeltaPhiFakeMET", "DeltaRFakeMET",
      "WMt", "FakeMt",
      "DeltaPhiWBosonWl", "DeltaRWBosonWl", "DeltaPhiWBosonFakel", "DeltaRWBosonFakel",
      "DeltaPhiFakeBosonWl", "DeltaRFakeBosonWl", "DeltaPhiFakeBosonFakel", "DeltaRFakeBosonFakel",
      "2LMass",
      "GoodJets", "DRminGoodJetWl", "DRminGoodJetFakel",
      "GoodJetsCut", "DRminGoodJetCutWl", "DRminGoodJetCutFakel",
      "GoodJetsLeadCut", "DRminGoodJetWlLeadCut", "DRminGoodJetFakelLeadCut",
      "LeadJetPt", "LeadJetEta", "LeadJetPhi", "LeadJetEt" };

  const string path = "/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/fom/GoodJetPt35GeV_LeadJetPt40GeV/NoRelIso/SSSelection/plots/stack/";
  vector<string> xAxisHisto =
    { "W lepton p_{t} [GeV]", "W lepton #eta", "W lepton #Phi", "W lepton relIso",
      "Fake lepton p_{t} [GeV]", "Fake lepton #eta", "Fake lepton #Phi", "Fake lepton relIso",
      "|#Delta#Phi (W l, Fake l)|", "#DeltaR (W l, Fake l)",
      "missing E_{t} [GeV]", "missing E_{t} #Phi" ,
      "|#Delta#Phi (W l, miss E_{t})|", "#DeltaR (W l, miss E_{t})", "|#Delta#Phi (Fake l, miss E_{t})|", "#DeltaR (Fake l, miss E_{t})",
       "W M_{t} [GeV]", "Fake M_{t} [GeV]",
       "|#Delta#Phi (W boson, W l)|", "#DeltaR (W boson, W l)", "|#Delta#Phi (W boson, Fake l)|", "#DeltaR (W boson, Fake l)",
       "|#Delta#Phi (Fake boson, W l)|", "#DeltaR (Fake boson, W l)", "|#Delta#Phi (Fake boson, Fake l)|", "#DeltaR (Fake boson, Fake l)",
       "mass_{2l} [GeV]",
       "Number of good Jets (p_{t}^{jet} > 10 GeV)",
       "#DeltaR_{min} (W l, good jet (p_{t}^{jet} > 10 GeV))", "#DeltaR_{min} (Fake l, good jet (p_{t}^{jet} > 10 GeV))",
       "Number of good Jets (p_{t}^{good jet} > 35 GeV)",
       "#DeltaR_{min} (W l, good jet (p_{t}^{good jet} > 35 GeV))", "#DeltaR_{min} (Fake l, good jet (p_{t}^{good jet} > 35 GeV))",
       "Number of good Jets (p_{t}^{lead jet} > 40 GeV)",
       "#DeltaR_{min} (W l, good jet (p_{t}^{lead jet} > 40 GeV))", "#DeltaR_{min} (Fake l, good jet (p_{t}^{lead jet} > 40 GeV))",
       "Lead jet p_{t} [GeV]", "Lead jet #eta", "Lead jet #phi", "Lead jet E_{t} [GeV]" };

  vector<double> binWidthHisto = { 2.0, 0.1, 0.087267, 0.002, 1.0, 0.1, 0.087267, 0.05, 0.04363, 0.05,
                                   2.0, 0.087267, 0.04363, 0.05, 0.04363, 0.05, 2.0, 2.0,
                                   0.04363, 0.04, 0.04363, 0.05, 0.04363, 0.05, 0.04363, 0.04,
                                   2.0,
                                   0., 0.01, 0.01 , 0., 0.05, 0.05 , 0., 0.01, 0.01,
                                   2.0, 0.1, 0.087267, 2.0 };

  for (unsigned int histo = 0; histo < histoName.size(); histo++) {
    TCanvas* canvas[n];
    for (unsigned int i = 0; i < n; i++) {
      canvas[i] = new TCanvas((histoName.at(histo) + "_" + boost::lexical_cast<string>(i)).c_str(),
                              histoName.at(histo).c_str());
      canvas[i]->cd();
      MakeStack(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT, f_DYM10, f_WJets,
                ("h" + histoName.at(histo) + "_" + boost::lexical_cast<string>(i)).c_str(),
                i, xAxisHisto.at(histo).c_str(), binWidthHisto.at(histo));
      canvas[i]->SaveAs((path + histoName.at(histo) + "_" +
                        boost::lexical_cast<string>(i) + ".pdf").c_str());
      delete canvas[i];
    }
  }
/*
  for (unsigned int jets = 0; jets < jetsName.size(); jets++) {
    TCanvas* canvas[n];
    for (unsigned int i = 0; i < n; i++) {
      canvas[i] = new TCanvas((jetsName.at(jets) + "_" + boost::lexical_cast<string>(i)).c_str(),
                              jetsName.at(jets).c_str());
      canvas[i]->cd();
      MakeStack(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT, f_DYM10, f_WJets,
                ("h" + jetsName.at(jets) + "_" + boost::lexical_cast<string>(i)).c_str(),
                i, xAxisJets.c_str(), 0.);
      canvas[i]->SaveAs((path + jetsName.at(jets) + "_" +
                        boost::lexical_cast<string>(i) + ".pdf").c_str());
      delete canvas[i];
    }
  }
*/
/*
  vector<string> deltaRJets = { "DeltaRL", "DeltaRMu", "DeltaREle" };
  vector<string> xAxisDeltaR =
    { "#deltaR_{min}(jet, lepton)", "#deltaR_{min}(jet, muon)", "#deltaR_{min}(jet, electron)" };
  for (unsigned int i = 0; i < deltaRJets.size(); i++) {
    TCanvas* canvas = new TCanvas(deltaRJets.at(i).c_str(), deltaRJets.at(i).c_str());
    canvas->cd();
    double max = 10000.;
    i ? max = 6000. : max = 12000.;
    MakeStackDeltaR(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT, f_DYM10, f_WJets
                    ("h" + deltaRJets.at(i)).c_str(), xAxisDeltaR.at(i).c_str(), max);
    canvas->SaveAs((path + deltaRJets.at(i) + ".pdf").c_str());
    delete canvas;
  }
*/

  const string pathGraph = "/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/fom/GoodJetPt35GeV_LeadJetPt40GeV/NoRelIso/SSSelection/plots/cuts/";
  vector<string> graphName = { "LeadJetPt", "LeadJetEta", "LeadJetPhi", "LeadJetEt" };
  vector<string> xAxisGraph = { "Lead jet p_{t} [GeV]", "Lead jet #eta", "Lead jet #phi", "Lead jet E_{t} [GeV]" };
  vector<double> startGraph = { 0., -2.5, -3.1416, 0. };
  vector<double> endGraph = { 250., 2.5, 3.1416, 250. };
  vector<double> binWidthGraph = { 2., 0.1, 0.087267, 2. };
  for (unsigned int graph = 0; graph < graphName.size(); graph++) {
    TCanvas* canvas[n];
    for (unsigned int i = 0; i < n; i++) {
      canvas[i] = new TCanvas((graphName.at(graph) + "_" + boost::lexical_cast<string>(i)).c_str(),
                              graphName.at(graph).c_str());
      canvas[i]->cd();
      MakeSignificanceGraphs(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT, f_DYM10, f_WJets,
                             ("h" + graphName.at(graph) + "_" + boost::lexical_cast<string>(i)).c_str(),
                             i, canvas[i], xAxisGraph.at(graph).c_str(),
                             startGraph.at(graph), endGraph.at(graph), binWidthGraph.at(graph));
      canvas[i]->SaveAs((pathGraph + graphName.at(graph) + "_" +
                        boost::lexical_cast<string>(i) + ".pdf").c_str());
      delete canvas[i];
    }
  }

  return 1;
}

