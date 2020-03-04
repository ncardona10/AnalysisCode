/*                             __
                           .--()°"."
 Author: Nathalia Cardona "|, . ,"
                           !_-(_\
*/

#include "PhenoAnalyzer.h"
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <iomanip>
#include <bits/stdc++.h>
#include "Analysis/LeptonCounter.h"
#include "Analysis/Cuts.h"

using namespace std;


void writeCsv(int count, string path, string cut){
  ofstream outfile;
  outfile.open("/home/n.cardonac/AnalysisCode/PhenoAnalyzer_InvisibleVBFDM/counts.csv", ios_base::app); // append instead of overwrite
  outfile << path<<","<<cut<<","<<count<<"\n"; 
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

  // TDirectory *single_e = HistoOutputFile->mkdir("single_e");
  // TDirectory *single_mu = HistoOutputFile->mkdir("single_mu");
  // TDirectory *single_tau = HistoOutputFile->mkdir("single_tau");

  TDirectory *METDirectory = HistoOutputFile->mkdir("MET");
  TDirectory *BJetsDirectory = HistoOutputFile->mkdir("BJets");
  TDirectory *VBFDirectory = HistoOutputFile->mkdir("VBF");
  TDirectory *ZeroLeptonsDirectory = HistoOutputFile->mkdir("ZeroLeptons");
  
  // TDirectory *single_e_met_bjets_vbf = HistoOutputFile->mkdir("single_e_met_bjets_vbf");
  // TDirectory *single_mu_met_bjets_vbf = HistoOutputFile->mkdir("single_mu_met_bjets_vbf");
  // // TDirectory *single_tau_met_bjets_vbf = HistoOutputFile->mkdir("single_tau_met_bjets_vbf");

  // TDirectory *di_e_met_bjets_vbf = HistoOutputFile->mkdir("di_e_met_bjets_vbf");
  // TDirectory *di_mu_met_bjets_vbf = HistoOutputFile->mkdir("di_mu_met_bjets_vbf");
  // TDirectory *di_tau_met_bjets_vbf = HistoOutputFile->mkdir("di_tau_met_bjets_vbf");

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
  vector<int> vbfCutsArr;
  vector<int> cutsArr;
  vector<int> vbfCutsArr_nocuts;
  vector<int> cutsArr_nocuts;

  for(int i = 0 ; (unsigned) i < treeReader->GetEntries(); i ++)
  {
    vbfCutsArr.push_back(1);
    cutsArr.push_back(1);
    vbfCutsArr_nocuts.push_back(1);
    cutsArr_nocuts.push_back(1);
  }

  int nEvents ;

  // write number of events to csv
  nEvents = (int) treeReader->GetEntries();
  writeCsv(nEvents, string(argv[1]),"C0");

  // open output file
  HistoOutputFile->cd();

  nLeptonsDirectory->cd();
  cout << "nLeptons" << endl;
  //drawLeptonCount(treeReader, ns, branchDict, vbfCutsArr, cutsArr,noFilter);
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, vbfCutsArr, cutsArr,noFilter);
  cout << "nLeptons done." << endl;

  writeCsv(nEvents, string(argv[1]),"nLeptons");

  // -----------------------------------------------------------------------------------------
  // need a different boolean array to avoid filtering problems

  // single_e->cd();
  // cout << "single_e" << endl;
  // //drawLeptonCount(treeReader, ns, branchDict, vbfCutsArr_nocuts, cutsArr_nocuts, mono_e);
  // nEvents = ptEtaPhiMjjMt(treeReader, branchDict, vbfCutsArr_nocuts, cutsArr_nocuts, mono_e);
  // cout << "single_e done." << endl;

  // writeCsv(nEvents, string(argv[1]),"single_e");

  // single_mu->cd();
  // cout << "single_mu" << endl;
  // //drawLeptonCount(treeReader, ns, branchDict, vbfCutsArr_nocuts, cutsArr_nocuts,mono_mu);
  // nEvents = ptEtaPhiMjjMt(treeReader, branchDict, vbfCutsArr_nocuts, cutsArr_nocuts,mono_mu);
  // cout << "single_mu done." << endl;

  // writeCsv(nEvents, string(argv[1]),"single_mu");


  // single_tau->cd();
  // cout << "single_tau" << endl;
  // //drawLeptonCount(treeReader, ns, branchDict, vbfCutsArr_nocuts, cutsArr_nocuts, mono_tau);
  // nEvents = ptEtaPhiMjjMt(treeReader, branchDict, vbfCutsArr_nocuts, cutsArr_nocuts, mono_tau);
  // cout << "single_tau done." << endl;

  // writeCsv(nEvents, string(argv[1]),"single_tau");

  // -----------------------------------------------------------------------------------------

  METDirectory->cd();
  cout << "MET" << endl;
  //drawLeptonCount(treeReader, ns, branchDict, vbfCutsArr, cutsArr, met);
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, vbfCutsArr, cutsArr, met);
  cout << "MET done." << endl;

  writeCsv(nEvents, string(argv[1]),"MET");

  BJetsDirectory->cd();
  cout << "BJets" << endl;
  //drawLeptonCount(treeReader, ns, branchDict, vbfCutsArr, cutsArr, bjets);
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, vbfCutsArr, cutsArr, bjets);
  cout << "BJets done." << endl;

  writeCsv(nEvents, string(argv[1]),"BJets");
  
  VBFDirectory->cd();
  cout << "VBF" << endl;
  //drawLeptonCount(treeReader, ns, branchDict, vbfCutsArr, cutsArr, vbfCut);
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, vbfCutsArr, cutsArr, vbfCut);
  cout << "VBF done." << endl;

  writeCsv(nEvents, string(argv[1]),"VBF");

  ZeroLeptonsDirectory->cd();
  cout << "ZeroLeptons" << endl;
  //drawLeptonCount(treeReader, ns, branchDict, vbfCutsArr, cutsArr, zero_leptons);
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, vbfCutsArr, cutsArr, zero_leptons);
  cout << "ZeroLeptons done." << endl;

  writeCsv(nEvents, string(argv[1]),"ZeroLeptons");

  // // -----------------------------------------------------------------------------------------

  // single_e_met_bjets_vbf->cd();
  // cout << "single_e_met_bjets_vbf" << endl;
  // //drawLeptonCount(treeReader, ns, branchDict, vbfCutsArr, cutsArr, mono_e);
  // nEvents = ptEtaPhiMjjMt(treeReader, branchDict, vbfCutsArr, cutsArr, mono_e);
  // cout << "single_e_met_bjets_vbf done." << endl;

  // writeCsv(nEvents, string(argv[1]),"single_e_met_bjets_vbf");

  // single_mu_met_bjets_vbf->cd();
  // cout << "single_mu_met_bjets_vbf" << endl;
  // //drawLeptonCount(treeReader, ns, branchDict, vbfCutsArr, cutsArr,mono_mu);
  // nEvents = ptEtaPhiMjjMt(treeReader, branchDict, vbfCutsArr, cutsArr,mono_mu);
  // cout << "single_mu_met_bjets_vbf done." << endl;
  
  // writeCsv(nEvents, string(argv[1]),"single_mu_met_bjets_vbf");

  // single_tau_met_bjets_vbf->cd();
  // cout << "single_tau_met_bjets_vbf" << endl;
  // //drawLeptonCount(treeReader, ns, branchDict, vbfCutsArr, cutsArr, mono_tau);
  // nEvents = ptEtaPhiMjjMt(treeReader, branchDict, vbfCutsArr, cutsArr, mono_tau);
  // cout << "single_tau_met_bjets_vbf done." << endl;
  
  // writeCsv(nEvents, string(argv[1]),"single_tau_met_bjets_vbf");

  // -----------------------------------------------------------------------------------------

  // di_e_met_bjets_vbf->cd();
  // cout << "di_e_met_bjets_vbf" << endl;
  // //drawLeptonCount(treeReader, ns, branchDict, vbfCutsArr, cutsArr, di_e);
  // nEvents = ptEtaPhiMjjMt(treeReader, branchDict, vbfCutsArr, cutsArr, di_e);
  // cout << "di_e_met_bjets_vbf done." << endl;

  // writeCsv(nEvents, string(argv[1]),"di_e_met_bjets_vbf");

  // di_mu_met_bjets_vbf->cd();
  // cout << "di_mu_met_bjets_vbf" << endl;
  // //drawLeptonCount(treeReader, ns, branchDict, vbfCutsArr, cutsArr,di_mu);
  // nEvents = ptEtaPhiMjjMt(treeReader, branchDict, vbfCutsArr, cutsArr,di_mu);
  // cout << "di_mu_met_bjets_vbf done." << endl;

  // writeCsv(nEvents, string(argv[1]),"di_mu_met_bjets_vbf");

  // // di_tau_met_bjets_vbf->cd();
  // cout << "di_tau_met_bjets_vbf" << endl;
  // //drawLeptonCount(treeReader, ns, branchDict, vbfCutsArr, cutsArr, di_tau);
  // nEvents = ptEtaPhiMjjMt(treeReader, branchDict, vbfCutsArr, cutsArr, di_tau);
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
