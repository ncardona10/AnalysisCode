/*                             __
                           .--()°'.'
 Author: Nathalia Cardona '|, . ,'
                           !_-(_\
*/

#include "../ROOTFunctions.h"
#include "../DelphesFunctions.h"

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
        leptonCount += 1.0;
      }
    }

    // muons
    for (int leaf = 0; leaf < branchDict["Muon"]->GetEntries(); leaf++)
    {
      Muon *lepton = (Muon *)branchDict["Muon"]->At(leaf);
      if (lepton->PT > 8.0 && lepton->PT < PTUpperCut)
      {
        leptonCount += 1.0;
      }
    }

    // taus
    for (int leaf = 0; leaf < branchDict["Jet"]->GetEntries(); leaf++)
    {
      Jet *jet = (Jet *)branchDict["Jet"]->At(leaf);
      if (jet->TauTag == 1)
      {
        if (jet->PT > 8.0 && jet->PT < PTUpperCut)
        {
          leptonCount += 1.0;
        }
      }
    }

    // write histograms
    nLeptonHistogram->Fill(leptonCount);
  }

  nLeptonHistogram->Write();

  cout << endl;
}

