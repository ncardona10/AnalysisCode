#include "PhenoAnalyzer.h"
#include <bits/stdc++.h>
#include "MyHistograms.h"

using namespace std;

int main(int argc, char const *argv[])
{
  cout << "creating output file" << endl;
  TFile *HistoOutputFile = new TFile("/home/n.cardonac/AnalysisCode_SUSY_VBF_Higgsino/plots/results/finalHistosnoetanomjj.root", "RECREATE");

  vector<TString> inputFolders = {"/home/n.cardonac/AnalysisCode_SUSY_VBF_Higgsino/RunCode/outPuts/outputfilesnoetanomjjtest/",
                                  "/home/n.cardonac/AnalysisCode_SUSY_VBF_Higgsino/RunCode/outPuts/outputfileseta55nomjjtest/"};

  vector<TString> inputFolderHistoName = {"w/o deta","with deta"};

  vector<TString> folderNames = {
      "m_n2_100_c1_80_n1_60",
      "m_n2_100_c1_75_n1_50",
      "m_n2_400_c1_385_n1_370",
      "m_n2_200_c1_175_n1_150",
      "wz",
      "zz",
      "ww",
      "w+jets",
      "z+jets",
      "ttbar"};

  map<TString, TFile *> sfiles;

  cout << "reading input files..." << endl;

  for (int i = 0; (unsigned)i < folderNames.size(); i++)
  {
    for (int j = 0; (unsigned)j < inputFolders.size(); j++)
    {
      sfiles[inputFolders[j] + folderNames[i]] = new TFile(inputFolders[j] + folderNames[i] + "/" + folderNames[i] + "_nocut-1.root", "READ");
      cout << "reading file: " << inputFolders[j] + folderNames[i] + "/" + folderNames[i] + "_nocut-1.root" << endl;
    }
  }

  vector<TString> dirNames = {
      // "nLeptons",
      // "MET",
      // "BJets",
      "VBF",
      // "single_e_met_bjets_vbf",
      // "single_mu_met_bjets_vbf",
      // "di_e_met_bjets_vbf",
      // "di_mu_met_bjets_vbf"
  };

  //create ouput folders in output root
  cout << "creating output subdirectories" << endl;
  HistoOutputFile->cd();
  map<TString, TDirectory *> outputDirs;
  for (int i = 0; (unsigned)i < folderNames.size(); i++)
  {
    outputDirs[folderNames[i]] = HistoOutputFile->mkdir(folderNames[i]);
  }

  //histogram names
  vector<TString> histoNames = {
      // "# of leptons PT < 15",
      //                               "# of electrons PT < 15",
      //                               "# of muons PT < 15",
      //                               "# of taus PT < 15",
      //                               "# of leptons PT < 20",
      //                               "# of electrons PT < 20",
      //                               "# of muons PT < 20",
      //                               "# of taus PT < 20",
      //                               "# of leptons PT < 30",
      //                               "# of electrons PT < 30",
      //                               "# of muons PT < 30",
      //                               "# of taus PT < 30",
      //                               "# of leptons PT < 40",
      //                               "# of electrons PT < 40",
      //                               "# of muons PT < 40",
      //                               "# of taus PT < 40",
      //                               "# of leptons PT < 50",
      //                               "# of electrons PT < 50",
      //                               "# of muons PT < 50",
      //                               "# of taus PT < 50",
      //                               "ptelectron",
      //                               "ptmuon",
      //                               "pttau",
      //                               "ptjet",
      //                               "etaelectron",
      //                               "etamuon",
      //                               "etatau",
      //                               "etajet",
      //                               "phielectron",
      //                               "phimuon",
      //                               "phitau",
      //                               "phijet",
      //                               "Mtelectron",
      //                               "Mtmuon",
      //                               "Mttau",
      //                               "Mtjet",
      "Mjj",
      // "MET"
  };

  /*

for every folder( ww,wz,...):
  for every directory (nleptons,...):
    for every histogram in folderhistogram in folder:
      get histogram
      change name to folder
    plot

*/

  /*
for folderNames

  get both input file histos

  graph them in one
*/

  for (int i = 0; (unsigned) i < folderNames.size(); i++)
  {
    cout<<i+1<<"/"<<folderNames.size()<<endl;
    cout << "processing " << folderNames[i] << endl;
    TObjArray histos;

    // get the data

    for (int j = 0; (unsigned) j < inputFolders.size(); j++)
    {

      // cout << "changing to " << folderNames[folder_i] << " directory and " << dirNames[directory_i] << " subfolder" << endl;
      sfiles[inputFolders[j] + folderNames[i]]->cd("VBF");

      cout << "getting histogram info" << endl;
      TH1F *histo = (TH1F *)sfiles[inputFolders[j] + folderNames[i]]->Get("VBF/Mjj");

      cout << "setting title and name" << endl;
      histo->SetTitle(inputFolderHistoName[j]);

      // histo->SetName(folderNames[folder_i]);
      histo->SetName(inputFolders[j]);

      cout << "adding to histogram list" << endl;
      histos.AddLast(histo);
    }

    outputDirs[folderNames[i]]->cd();

    TCanvas *cl = new TCanvas("VBF Mjj", "VBF Mjj", 600, 500);

    // cl->Divide(2,2); //create subplots if needed

    string histoNameString = (string)histoNames[0];
    Draw_Normalised(histos, (TPad *)cl->cd(0), true, histoNameString);

    cl->Write();

    cout<<"done with "<<folderNames[i] << endl;
  }


  HistoOutputFile->Close();

  for (int i = 0; (unsigned)i < folderNames.size(); i++)
  {
    for (int j = 0; (unsigned)j < inputFolders.size(); j++)
    {
      cout << "closing: " << inputFolders[j] + folderNames[i] << endl;
      sfiles[inputFolders[j] + folderNames[i]]->Close();
    }
  }
  cout << "DONE" << endl;

  cout << "all good." << endl;
  return 0;
}
