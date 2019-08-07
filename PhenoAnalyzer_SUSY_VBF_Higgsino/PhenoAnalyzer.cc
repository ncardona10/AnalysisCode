/*                             __
                           .--()°'.'
 Author: Nathalia Cardona '|, . ,'
                           !_-(_\
*/

#include "PhenoAnalyzer.h"
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <iomanip>
#include <bits/stdc++.h>
#include "Analysis/LeptonCounter.h"
#include "Analysis/VBF_Cuts.h"

using namespace std;

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
  TDirectory *nLeptonsDirectory = HistoOutputFile->mkdir("nLeptons");
  TDirectory *VBF_CutsDirectory = HistoOutputFile->mkdir("VBF_Cuts");

  cout << "processing.." << endl;

  // get tree info
  vector<string> branches = {
      "Electron",
      "Muon",
      "Jet",
      "MissingET"};

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
  // open output file
  HistoOutputFile->cd();

  nLeptonsDirectory->cd();

  drawLeptonCount(treeReader, ns, branchDict, noFilter);

  ptEtaPhi(treeReader, branchDict, noFilter);


  VBF_CutsDirectory->cd();

  drawLeptonCount(treeReader, ns, branchDict, cut1);

  ptEtaPhi(treeReader, branchDict, cut1);

  // close output file
  HistoOutputFile->Close();

  cout << "DONE." << endl;
}