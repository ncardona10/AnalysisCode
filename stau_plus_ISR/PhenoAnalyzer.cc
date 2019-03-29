/*
@file PhenoAnalyzer.cc
@author Andres Florez
Code used to perform phenomenological analysis of Heavy Neutrinos in the tau channel
*/

#include "PhenoAnalyzer.h"
#include <string>

int main(int argc, char *argv[])
{

  //TApplication app("App",&argc, argv);
  // gSystem->Load("libDelphes.so");
  TChain chain("Delphes");
  chain.Add(argv[1]);
  TFile *HistoOutputFile = new TFile(argv[2], "RECREATE");
  int nDir = 1;
  TDirectory *theDirectory[nDir];
  theDirectory[0] = HistoOutputFile->mkdir("No_cuts");
  PhenoAnalysis BSM_analysis(chain, HistoOutputFile, theDirectory, nDir);
}

using namespace std;
PhenoAnalysis::PhenoAnalysis(TChain &chain, TFile *theFile, TDirectory *cdDir[], int nDir)
{
  ifstream inFile;
  inFile.open("config.in", ios::in);

  if (!inFile)
  {
    cerr << "ERROR: Can't open input file: " << endl;
    exit(1);
  }

  string inputType = "";

  //This set of lines are used to open and read the "config.in" file.
  ///////////////////////////////////////////////////////////////////////
  TEnv *params = new TEnv("config_file");
  params->ReadFile("config.in", kEnvChange);

  double b_jet_pt_min = params->GetValue("b_jet_pt_min", 30.0);

  createHistoMaps(nDir);

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");

  MissingET *METpointer; //Energia transversal perdida

  for (Int_t entry = 0; entry < numberOfEntries; ++entry)
  {

    treeReader->ReadEntry(entry);
    int pass_cuts[nDir];
    TLorentzVector jetLeadingVec(0., 0., 0., 0.); //tc stands for tau channel

    //Search for taus and bjets
    for (int j = 0; j < branchJet->GetEntriesFast(); j++)
    {
      Jet *jet = (Jet *)branchJet->At(j);
    }
    //Check if taus overlap with muons

    for (int muo = 0; muo < branchMuon->GetEntriesFast(); muo++)
    {

    } // for (int muo = 0; muo < branchMuon->GetEntriesFast(); muo++)

    //Check if taus overlap with electrons

    for (int el = 0; el < branchElectron->GetEntriesFast(); el++)
    {
    }

    //////// Apply cuts /////////
    pass_cuts[0] = 1;

    //Fill histograms (MM: No sé como llenar los histogramas de masa...)
    for (Int_t i = 0; i < nDir; i++)
    {
      _hmap_Nevents[i]->Fill(0.0);
    } // end entry loop for tau channel

    theFile->cd();
    for (int d = 0; d < nDir; d++)
    {
      cdDir[d]->cd();
      _hmap_Nevents[d]->Write();
    }
    theFile->Close();
  }
}
PhenoAnalysis::~PhenoAnalysis()
{
  // do anything here that needs to be done at desctruction time
}

double PhenoAnalysis::calculateE(double eta, double pt, double mass)
{

  double theta = TMath::ATan(TMath::Exp(-eta));
  double sin_theta = TMath::Sin(2 * theta);
  double p = pt / sin_theta;
  double e = sqrt(pow(p, 2) + pow(mass, 2));

  return e;
}

double PhenoAnalysis::calculate_deltaR(TLorentzVector vector, Track *track)
{

  double eta1 = vector.Eta();
  double phi1 = vector.Phi();
  double eta2 = track->Eta;
  double phi2 = track->Phi;
  double deltaR = sqrt(pow(eta1 - eta2, 2) + pow(phi1 - phi2, 2));
  return deltaR;
}

double PhenoAnalysis::normalizedDphi(double phi)
{
  const double PI = 3.141592653589793238463;
  double twoPI = 2.0 * PI;
  if (phi < -PI)
  {
    phi += twoPI;
  }
  if (phi > PI)
  {
    phi = twoPI - phi;
  }
  else
    phi = TMath::Abs(phi);
  return phi;
}
void PhenoAnalysis::createHistoMaps(int directories)
{
  for (int i = 0; i < directories; i++)
  {
    std::string si = "N_events_" + std::to_string(i);
    char char_array[si.length() + 1];
    strcpy(char_array, si.c_str());
    printf(char_array);
    printf("\n");
    _hmap_Nevents[i] = new TH1F(char_array, char_array, 3, 0.0, 3.0);
  }
  printf("finished createHistoMaps\n`");
}
