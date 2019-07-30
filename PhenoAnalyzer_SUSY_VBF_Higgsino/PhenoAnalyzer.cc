/*
@file PhenoAnalyzer.cc
@author Andres Florez
@author Carlos Miguel Patino
@author Nathalia Cardona 
@date April 2, 2017

Code used to perform phenomenological analysis of Heavy Neutrinos in the tau channel
*/

#include "PhenoAnalyzer.h"
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <iomanip>
#include <bits/stdc++.h>

using namespace std;

void nLeptonAnalysis(ExRootTreeReader *treeReader, TFile *theFile, TDirectory *cdDir, int PTUpperCut, map<string, TClonesArray *> branchDict)
{

  cout << "n = " << PTUpperCut << endl;

  Long64_t numberOfEntries = treeReader->GetEntries();

  // create histogram
  string si = "Number of Leptons, 8 < P_{T} < " + to_string(PTUpperCut);
  char char_array[si.length() + 1];
  strcpy(char_array, si.c_str());

  string histoName = "# of leptons PT < " + to_string(PTUpperCut);
  char charArrayHistoName[histoName.length() + 1];
  strcpy(charArrayHistoName, histoName.c_str()); 

  TH1 *nLeptonHistogram = new TH1F(charArrayHistoName, char_array, 15, 0.0, 15.0);

  

  // cout << "Number of entries: " << numberOfEntries << endl;

  for (Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    float leptonCount = 0;
    // print percentage of completion
    cout << "\r" << (100.0 * entry) / numberOfEntries << "%";
    treeReader->ReadEntry(entry);

    // electrons
    for (int leaf = 0; leaf < branchDict["Electron"]->GetEntries(); leaf++)
    {
      Electron *lepton = (Electron *)branchDict["Electron"]->At(leaf);
      if (lepton->PT > 8.0 && lepton->PT < PTUpperCut)
      {
        leptonCount +=1.0;
      }
    }

    // muons
    for (int leaf = 0; leaf < branchDict["Muon"]->GetEntries(); leaf++)
    {
      Muon *lepton = (Muon *)branchDict["Muon"]->At(leaf);
      if (lepton->PT > 8.0 && lepton->PT < PTUpperCut)
      {
        leptonCount +=1.0;
      }
    }

    // taus
    // for (int leaf = 0; leaf < branchDict["Tau"]->GetEntries(); leaf++)
    // {
    //   Tau *lepton = (Tau *)branchDict["Tau"]->At(leaf);
    //   if (lepton->PT > 8.0 && lepton->PT < PTUpperCut)
    //   {
    //     leptonCount +=1.0;
    //   }
    // }

    // write histograms
    nLeptonHistogram->Fill(leptonCount);
    
  }
  
  nLeptonHistogram->Write();

  cout << endl;
}

int main(int argc, char *argv[])
{
  cout << "Starting phenoanalyzer..." << endl;

  // standardize print to 2 dp
  cout << fixed;
  cout << setprecision(2);

  // Importing Delphes data
  TChain chain("Delphes");
  chain.Add(argv[1]);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);

  vector<int> ns = {15, 20, 30, 40, 50};

  // output file manager
  TFile *HistoOutputFile = new TFile(argv[2], "RECREATE");
  // directory to store the histograms
  TDirectory *theDirectory = HistoOutputFile->mkdir("nLeptons");

  cout << "processing.." << endl;

  // get tree info
  

  vector<string> branches = {
      "Electron",
      "Muon",
      // "Tau"
      };

  map<string, TClonesArray *> branchDict;
  // create a dictionary with the branches
  for (int i = 0; (unsigned)i < branches.size(); i++)
  {
    TClonesArray *branch = treeReader->UseBranch(branches[i].c_str());
    branchDict[branches[i]] = branch;
  }

  /*
    File Structure:
    - Test.root
      - nLeptons
        - 8 < pt < 50 
        - 8 < pt < 40
            .
            .
            .
  */
  HistoOutputFile->cd();
  theDirectory->cd();

  for (int i = 0; (unsigned) i < ns.size(); i++)
  {
    nLeptonAnalysis(treeReader, HistoOutputFile, theDirectory, ns[i], branchDict);
  }

  // close output file
  HistoOutputFile->Close();

  cout << "DONE." << endl;
}