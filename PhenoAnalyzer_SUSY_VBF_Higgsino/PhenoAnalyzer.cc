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
#include "Analysis/Cuts.h"

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
  TDirectory *CutsDirectory = HistoOutputFile->mkdir("Extra_Cuts");
  TDirectory *single_e = HistoOutputFile->mkdir("single_e");
  TDirectory *single_mu = HistoOutputFile->mkdir("single_mu");
  TDirectory *single_tau = HistoOutputFile->mkdir("single_tau");

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
  cout<<"nLeptons"<<endl;
  drawLeptonCount(treeReader, ns, branchDict, noFilter);
  ptEtaPhiMjjMt(treeReader, branchDict, noFilter);

  cout<<"nLeptons done."<<endl; 

  VBF_CutsDirectory->cd();
  cout<<"VBF_Cuts"<<endl;
  drawLeptonCount(treeReader, ns, branchDict, vbfCut);
  ptEtaPhiMjjMt(treeReader, branchDict, vbfCut);

  CutsDirectory->cd();
  cout<<"Extra_cuts"<<endl;
  drawLeptonCount(treeReader, ns, branchDict, cuts);
  ptEtaPhiMjjMt(treeReader, branchDict, cuts);

  single_e->cd();
  cout<<"single_e"<<endl;
  drawLeptonCount(treeReader, ns, branchDict, cut_e);
  ptEtaPhiMjjMt(treeReader, branchDict, cut_e);

  single_mu->cd();
  cout<<"single_mu"<<endl;
  drawLeptonCount(treeReader, ns, branchDict, cut_mu);
  ptEtaPhiMjjMt(treeReader, branchDict, cut_mu);

  single_tau->cd();
  cout<<"single_tau"<<endl;
  drawLeptonCount(treeReader, ns, branchDict, cut_tau);
  ptEtaPhiMjjMt(treeReader, branchDict, cut_tau);

  // close output file
  HistoOutputFile->Close();

  cout << "DONE." << endl;
}