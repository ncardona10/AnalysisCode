/*                             __
                           .--()°'.'
 Author: Nathalia Cardona '|, . ,'
                           !_-(_\
 Counts the number of leptons in different PT ranges.                        
*/

#include "../ROOTFunctions.h"
#include "../DelphesFunctions.h"
#include "../Plots/MyHistograms.h"
#include "../Analysis/Physics.h"
#include <string>
#include <map>

bool elecOverlap(ExRootTreeReader *treeReader,
                 map<string, TClonesArray *> branchDict,
                 Jet *jet)
{
  bool ans = false;

  double jetE = calculateE(jet->PT, jet->Eta, jet->Mass);
  TLorentzVector jetTLV(jet->PT, jet->Eta, jet->Phi, jetE);

  int leaf = 0;

  while (!ans && leaf < branchDict["Electron"]->GetEntries())
  {

    Electron *electron = (Electron *)branchDict["Electron"]->At(leaf);

    double electE = calculateE(electron->PT, electron->Eta, 0.000510998902);
    TLorentzVector elecTLV(electron->PT, electron->Eta, electron->Phi, electE);

    double dr = dR(jetTLV, elecTLV);

    ans = overlap(dr);

    leaf++;
  }
  return ans;
}

bool muonOverlap(ExRootTreeReader *treeReader,
                 map<string, TClonesArray *> branchDict,
                 Jet *jet)
{
  bool ans = false;

  double jetE = calculateE(jet->PT, jet->Eta, jet->Mass);
  TLorentzVector jetTLV(jet->PT, jet->Eta, jet->Phi, jetE);

  int leaf = 0;

  while (!ans && leaf < branchDict["Muon"]->GetEntries())
  {

    Muon *muon = (Muon *)branchDict["Muon"]->At(leaf);

    double muonE = calculateE(muon->PT, muon->Eta, 0.1056583715);
    TLorentzVector muonTLV(muon->PT, muon->Eta, muon->Phi, muonE);

    double dr = dR(jetTLV, muonTLV);

    ans = overlap(dr);

    leaf++;
  }
  return ans;
}

void fillHisto(TH1 *histo, float value)
{
  if (value > 0)
  {
    histo->Fill(value);
  }
}

map<string, TH1 *> nLeptonAnalysis(ExRootTreeReader *treeReader,
                                   int PTUpperCut,
                                   map<string, TClonesArray *> branchDict)
{

  cout << "n = " << PTUpperCut << endl;

  Long64_t numberOfEntries = treeReader->GetEntries();

  // create histograms

  //lepton histogram
  TH1 *nLeptonHistogram = blankHistogram("Number of Leptons, 8 < P_{T} < " + to_string(PTUpperCut),
                                         "# of leptons PT < " + to_string(PTUpperCut),
                                         15, 0.0, 15.0);

  // electron  histogram
  TH1 *nElecHistogram = blankHistogram("Number of Electrons, 8 < P_{T} < " + to_string(PTUpperCut),
                                       "# of electrons PT < " + to_string(PTUpperCut),
                                       15, 0.0, 15.0);

  // muon  histogram
  TH1 *nMuonHistogram = blankHistogram("Number of Muons, 8 < P_{T} < " + to_string(PTUpperCut),
                                       "# of muons PT < " + to_string(PTUpperCut),
                                       15, 0.0, 15.0);

  // tau  histogram
  TH1 *nTauHistogram = blankHistogram("Number of Taus, 8 < P_{T} < " + to_string(PTUpperCut),
                                      "# of taus PT < " + to_string(PTUpperCut),
                                      15, 0.0, 15.0);

  map<string, TH1 *> histograms = {{"lepton", nLeptonHistogram},
                                   {"electron", nElecHistogram},
                                   {"muon", nMuonHistogram},
                                   {"tau", nTauHistogram}};

  for (Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    float leptonCount = 0;
    float electronCount = 0;
    float muonCount = 0;
    float tauCount = 0;

    // print percentage of completion
    cout << "\r" << (100.0 * entry) / numberOfEntries << "%";
    treeReader->ReadEntry(entry);

    // electrons
    for (int leaf = 0; leaf < branchDict["Electron"]->GetEntries(); leaf++)
    {
      Electron *lepton = (Electron *)branchDict["Electron"]->At(leaf);
      if (lepton->PT > 8.0 && lepton->PT < PTUpperCut)
      {
        electronCount += 1.0;
        leptonCount += 1.0;
      }
    }

    // muons
    for (int leaf = 0; leaf < branchDict["Muon"]->GetEntries(); leaf++)
    {
      Muon *lepton = (Muon *)branchDict["Muon"]->At(leaf);
      if (lepton->PT > 8.0 && lepton->PT < PTUpperCut)
      {
        muonCount += 1.0;
        leptonCount += 1.0;
      }
    }

    // taus
    for (int leaf = 0; leaf < branchDict["Jet"]->GetEntries(); leaf++)
    {
      Jet *jet = (Jet *)branchDict["Jet"]->At(leaf);
      if (jet->TauTag == 1)
      {
        if (!elecOverlap(treeReader, branchDict, jet))
        {
          if (!muonOverlap(treeReader, branchDict, jet))
          {
            // its a tau!
            if (jet->PT > 8.0 && jet->PT < PTUpperCut)
            {
              leptonCount += 1.0;
              tauCount += 1.0;
            }
          }
          else
          {
            //muon overlap
            muonCount -= 1.0;
            leptonCount -= 1.0;
          }
        }
        else
        {
          //electron overlap
          electronCount -= 1.0;
          leptonCount -= 1.0;
        }
      }
    }

    // write histograms

    fillHisto(nLeptonHistogram, leptonCount);
    fillHisto(nElecHistogram, electronCount);
    fillHisto(nMuonHistogram, muonCount);
    fillHisto(nTauHistogram, tauCount);
  }

  cout << endl;
  return histograms;
}

void drawMultiHistos(TObjArray histos, string title, string particleType)
{

  char charTitle[title.length() + 1];
  strcpy(charTitle, title.c_str());

  TCanvas *cl = new TCanvas(charTitle, charTitle, 600, 500);

  // cl->Divide(2,2); //create subplots

  Draw_Normalised(histos, (TPad *)cl->cd(0), false, "# of " + particleType + "s under different P_{T} cuts", 10);

  cl->Write();
}

void drawLeptonCount(ExRootTreeReader *treeReader, vector<int> ns, map<string, TClonesArray *> branchDict)
{

  vector<string> particleTypes = {"lepton", "electron", "muon", "tau"};
  map<string, TObjArray> histos;

  for (int i = 0; (unsigned)i < particleTypes.size(); i++)
  {
    histos[particleTypes[i]] = TObjArray();
  }

  for (int i = 0; (unsigned)i < ns.size(); i++)
  {
    map<string, TH1 *> histoOutput = nLeptonAnalysis(treeReader, ns[i], branchDict);

    for (int j = 0; (unsigned)j < particleTypes.size(); j++)
    {
      histos[particleTypes[j]].AddLast(histoOutput[particleTypes[j]]);
    }
  }

  for (int i = 0; (unsigned)i < particleTypes.size(); i++)
  {
    drawMultiHistos(histos[particleTypes[i]], "#" + particleTypes[i] + "s", particleTypes[i]);
  }
}
