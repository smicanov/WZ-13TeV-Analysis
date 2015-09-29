#include "WJetsSelection.h"
#include "MyStyle.h"

#include <string>


using namespace std;


int main(int argc, char **argv)
{
  char* outputFileName(0);
  char* inputFileName(0);
  char* fileList(0);
  bool useInputList = false;
  bool gotInput  = false;
  bool gotOutput = false;
  char c;

  while ((c = getopt (argc, argv, "i:o:l:")) != -1)
    switch (c) {
      case 'o':
        gotOutput = true;
        outputFileName = new char[strlen(optarg)+1];
        strcpy(outputFileName, optarg);
        break;
      case 'i':
        gotInput = true;
        inputFileName = new char[strlen(optarg)+1];
        strcpy(inputFileName, optarg);
        break;
      case 'l':
        useInputList = true;
        fileList = new char[strlen(optarg)+1];
        strcpy(fileList, optarg);
        break;
      default:
        cout << "usage: [-k|-g|-l] [-v] [-b <binWidth>]   -i <input> -o <output> \n";
        abort ();
    }

// OUTPUT ROOT FILE
  TFile* outputFile;
  if (gotOutput)  outputFile = new TFile(outputFileName, "RECREATE");
  else            outputFile = new TFile("/users/msasa/work/cms/wz/ggAna/code/WZ-13TeV-Analysis/output/fom/default.root", "RECREATE");

// INPUT TREES
  vector<TString> inputName;
  TChain wjets("ggNtuplizer/EventTree");

  if (useInputList) {
    ifstream list(fileList);
    TString name;
    while (list>>name)  inputName.push_back(name);
  } else if (gotInput)  inputName.push_back(inputFileName);
  else {
    cout << "Got no input ROOT file: quit \n";
    return 1;
  }

  for (unsigned int input = 0; input < inputName.size(); input++) {
    cout << "Adding: " << input << endl;
    wjets.Add(inputName[input]);
    cout << "Added \n";
  }

  TTree* wjets_tTree = (TTree*)& wjets;
  Event* eventWJets = new Event(wjets_tTree);
  Int_t events = wjets_tTree->GetEntries();
  cout << endl << "Total number of events: " << events << endl << endl;
  const MyStyle rootStyle(800, 1);

  WJetsSelection* selection = new WJetsSelection(eventWJets, outputFile);
  selection->Init();

// Event loop
  unsigned int nEvents = 0;

  for  (Int_t k = 0; k < events /*&& k < 200000*/; k++) {
    cerr <<  "  " << int(100 * 100 * (k+1.) / events + 0.5) / 100. << " %            \r" << flush;
    nEvents++;
    wjets_tTree->GetEntry(k);
    eventWJets->Read();

    selection->Analysis();
  }

  cerr << "  100%" << endl << endl;

  cout << "Events : " << nEvents << endl;

  selection->Finish();
  selection->WriteRootFile();

  outputFile->Close();
}

