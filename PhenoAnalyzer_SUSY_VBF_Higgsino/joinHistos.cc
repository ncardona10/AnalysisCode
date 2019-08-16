#include "PhenoAnalyzer.h"
#include <bits/stdc++.h>
#include "./Plots/MyHistograms.h"

using namespace std;

int main(int argc, char const *argv[])
{
  TFile *HistoOutputFile = new TFile("FinalPlots.root", "RECREATE");
  TDirectory *outputDir = HistoOutputFile->mkdir("finalPlots");

  TFile *sfile = new TFile("Test.root");

  vector<TString> dirNames = {"nLeptons",
                              "VBF_Cuts",
                              "Extra_Cuts",
                              "single_e",
                              "single_mu",
                              "single_tau"};

  sfile->cd("nLeptons");

  //get histogram names
  vector<TString> histoNames = {"# of leptons PT < 15",
                                "# of electrons PT < 15",
                                "# of muons PT < 15",
                                "# of taus PT < 15",
                                "# of leptons PT < 20",
                                "# of electrons PT < 20",
                                "# of muons PT < 20",
                                "# of taus PT < 20",
                                "# of leptons PT < 30",
                                "# of electrons PT < 30",
                                "# of muons PT < 30",
                                "# of taus PT < 30",
                                "# of leptons PT < 40",
                                "# of electrons PT < 40",
                                "# of muons PT < 40",
                                "# of taus PT < 40",
                                "# of leptons PT < 50",
                                "# of electrons PT < 50",
                                "# of muons PT < 50",
                                "# of taus PT < 50",
                                "ptelectron",
                                "ptmuon",
                                "pttau",
                                "ptjet",
                                "etaelectron",
                                "etamuon",
                                "etatau",
                                "etajet",
                                "phielectron",
                                "phimuon",
                                "phitau",
                                "phijet",
                                "Mtelectron",
                                "Mtmuon",
                                "Mttau",
                                "Mtjet",
                                "Mjj",
                                "MET"};

  for (int i = 0; (unsigned)i < histoNames.size(); i++)
  {
    cout << i * 100.0 / histoNames.size() << endl;
    ;

    TObjArray histos;

    for (int j = 0; (unsigned)j < dirNames.size(); j++)
    {

      TString dirName = dirNames[j];

      sfile->cd(dirName);

      TH1F *histo = (TH1F *)sfile->Get(dirName + "/" + histoNames[i]);

      histo->SetTitle(dirName);
      histo->SetName(dirName);

      histos.AddLast(histo);
    }

    outputDir->cd();

    TCanvas *cl = new TCanvas(histoNames[i], histoNames[i], 600, 500);

    // cl->Divide(2,2); //create subplots

    string histoNameString = (string)histoNames[i];
    Draw_Normalised(histos, (TPad *)cl->cd(0), false, histoNameString);

    cl->Write();
  }

  HistoOutputFile->Close();
  sfile->Close();
  cout << endl;
  return 0;
}
