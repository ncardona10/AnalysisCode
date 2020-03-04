/*                             __
                           .--()°"."
 Author: Nathalia Cardona "|, . ,"
                           !_-(_\
*/

#include "PhenoAnalyzer.h"
#include "Analysis/LeptonCounter.h"
#include "Analysis/Cuts.h"
#include "Analysis/Physics.h"
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <iomanip>
#include <bits/stdc++.h>
#include <utility>

using namespace std;

void writeCsv(int count, string path, string cut)
{
  ofstream outfile;
  outfile.open("/home/n.cardonac/AnalysisCode/PhenoAnalyzer_SUSY_VBF_Higgsino/counts.csv", ios_base::app); // append instead of overwrite
  outfile << path << "," << cut << "," << count << "\n";
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
  TDirectory *nLeptonsDirectory = HistoOutputFile->mkdir("nLeptons");

  TDirectory *JetsDirectory = HistoOutputFile->mkdir("Jets");
  TDirectory *METDirectory = HistoOutputFile->mkdir("MET");
  TDirectory *BJetsDirectory = HistoOutputFile->mkdir("BJets");
  TDirectory *VBFDirectory = HistoOutputFile->mkdir("VBF");

  TDirectory *single_e_met_bjets_vbf = HistoOutputFile->mkdir("single_e_met_bjets_vbf");
  TDirectory *single_mu_met_bjets_vbf = HistoOutputFile->mkdir("single_mu_met_bjets_vbf");

  // TDirectory *di_e_met_bjets_vbf = HistoOutputFile->mkdir("di_e_met_bjets_vbf");
  // TDirectory *di_mu_met_bjets_vbf = HistoOutputFile->mkdir("di_mu_met_bjets_vbf");

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

  // boolean mask to avoid over computation
  vector<int> cutsArr;

  // vector of pairs for jet pair that max mjj
  vector<pair<Jet *, Jet *>> maxMjjJetPairs;

  cout<<"creating pairs"<<endl;

  for (int i = 0; (unsigned)i < treeReader->GetEntries(); i++)
  {
    cutsArr.push_back(1);
    maxMjjJetPairs.push_back(getMaxMjjJetPair(treeReader, branchDict, i));
  }

  cout<< "done creating pairs"<<endl;

  int nEvents;

  // write number of events to csv
  nEvents = (int)treeReader->GetEntries();
  writeCsv(nEvents, string(argv[1]), "C0");

  // open output file
  HistoOutputFile->cd();

  // ---------------------------------No cuts--------------------------------------------

  nLeptonsDirectory->cd();
  cout << "No cuts" << endl;
  //drawLeptonCount(treeReader, ns, branchDict, cutsArr,noFilter);
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, cutsArr,maxMjjJetPairs, noFilter);
  cout << "No cuts done." << endl;

  writeCsv(nEvents, string(argv[1]), "nLeptons");

  // ---------------------------------Jets--------------------------------------------

  JetsDirectory->cd();
  cout << "Jets" << endl;
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, cutsArr, maxMjjJetPairs, jetConditions);
  cout << "Jets done." << endl;

  writeCsv(nEvents, string(argv[1]), "Jets");

  // ---------------------------------BJets--------------------------------------------

  BJetsDirectory->cd();
  cout << "BJets" << endl;
  //drawLeptonCount(treeReader, ns, branchDict, cutsArr, bjets);
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, cutsArr, maxMjjJetPairs, bjets);
  cout << "BJets done." << endl;

  writeCsv(nEvents, string(argv[1]), "BJets");

    // ---------------------------------MET--------------------------------------------

  METDirectory->cd();
  cout << "MET" << endl;
  //drawLeptonCount(treeReader, ns, branchDict, cutsArr, met);
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, cutsArr, maxMjjJetPairs, met);
  cout << "MET done." << endl;

  writeCsv(nEvents, string(argv[1]), "MET");

  // ------------------------------------VBF----------------------------------------
  VBFDirectory->cd();
  cout << "VBF" << endl;
  //drawLeptonCount(treeReader, ns, branchDict, cutsArr, vbfCut);
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, cutsArr,maxMjjJetPairs,  vbfCut);
  cout << "VBF done." << endl;

  writeCsv(nEvents, string(argv[1]), "VBF");

  // ---------------------------------SINGLE E--------------------------------------------

  single_e_met_bjets_vbf->cd();
  cout << "single_e_met_bjets_vbf" << endl;
  //drawLeptonCount(treeReader, ns, branchDict, cutsArr, mono_e);
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, cutsArr, maxMjjJetPairs, mono_e);
  cout << "single_e_met_bjets_vbf done." << endl;

  writeCsv(nEvents, string(argv[1]), "single_e_met_bjets_vbf");

  // ---------------------------------SINGLE MU--------------------------------------------

  single_mu_met_bjets_vbf->cd();
  cout << "single_mu_met_bjets_vbf" << endl;
  //drawLeptonCount(treeReader, ns, branchDict, cutsArr,mono_mu);
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, cutsArr, maxMjjJetPairs, mono_mu);
  cout << "single_mu_met_bjets_vbf done." << endl;

  writeCsv(nEvents, string(argv[1]), "single_mu_met_bjets_vbf");

  // ---------------------------------DI LEPTON---------------------------------------------

  // di_e_met_bjets_vbf->cd();
  // cout << "di_e_met_bjets_vbf" << endl;
  // //drawLeptonCount(treeReader, ns, branchDict, cutsArr, di_e);
  // nEvents = ptEtaPhiMjjMt(treeReader, branchDict, cutsArr, di_e);
  // cout << "di_e_met_bjets_vbf done." << endl;

  // writeCsv(nEvents, string(argv[1]),"di_e_met_bjets_vbf");

  // di_mu_met_bjets_vbf->cd();
  // cout << "di_mu_met_bjets_vbf" << endl;
  // //drawLeptonCount(treeReader, ns, branchDict, cutsArr,di_mu);
  // nEvents = ptEtaPhiMjjMt(treeReader, branchDict, cutsArr,di_mu);
  // cout << "di_mu_met_bjets_vbf done." << endl;

  // writeCsv(nEvents, string(argv[1]),"di_mu_met_bjets_vbf");

  // di_tau_met_bjets_vbf->cd();
  // cout << "di_tau_met_bjets_vbf" << endl;
  // //drawLeptonCount(treeReader, ns, branchDict, cutsArr, di_tau);
  // nEvents = ptEtaPhiMjjMt(treeReader, branchDict, cutsArr, di_tau);
  // cout << "di_tau_met_bjets_vbf done." << endl;

  // writeCsv(nEvents, string(argv[1]),"di_tau_met_bjets_vbf");

  // --------------------------------------------------------------------------------------------

  // close output file
  cout << "closing output file" << endl;
  HistoOutputFile->Close();

  //write to file as log
  cout << "Writing to log file" << endl;
  ofstream myfile;
  myfile.open("finishedProcesses.dat", ios::app);
  myfile << argv[1] << "\n";
  myfile.close();

  cout << "DONE." << endl;
}
