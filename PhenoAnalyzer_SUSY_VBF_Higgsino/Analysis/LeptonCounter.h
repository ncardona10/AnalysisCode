/*                             __
                           .--()°'.'
 Author: Nathalia Cardona '|, . ,'
                           !_-(_\
 Counts the number of leptons in different PT ranges.                        
*/


#ifndef LEPTONCOUNTER_H
#define LEPTONCOUNTER_H

#include "../ROOTFunctions.h"
#include "../DelphesFunctions.h"
#include "../Plots/MyHistograms.h"
#include "./Physics.h"
#include "./Overlaps.h"
#include <string>
#include <map>
#include <vector>
#include <set>

void fillHisto(TH1 *histo, float value)
{
  if (value > 0)
  {
    histo->Fill(value);
  }
}

map<string, TH1 *> nLeptonAnalysis(ExRootTreeReader *treeReader,
                     int PTUpperCut,
                     map<string, TClonesArray *> branchDict,
                     vector<int> &cutsArr,
                     vector<pair<Jet*,Jet*>> &jetPairs,
                     bool (*filter)(ExRootTreeReader *,
                                    map<string, TClonesArray *>,
                                    int,
                                    vector<int> &,
                                    vector<pair<Jet*,Jet*>> &))
{

  //  cout << "nLepton Analysis with n = " << PTUpperCut << endl;

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
    //    // cout << "\r" << (100.0 * entry) / numberOfEntries << "%";
    treeReader->ReadEntry(entry);

    if (filter(treeReader, branchDict, entry, cutsArr,jetPairs))
    {

      // electrons
      for (int leaf = 0; leaf < branchDict["Electron"]->GetEntries(); leaf++)
      {
        Electron *lepton = (Electron *)branchDict["Electron"]->At(leaf);
        if (lepton->PT > 8.0 && lepton->PT < PTUpperCut && abs(lepton->Eta) < 2.4)
        {
          electronCount += 1.0;
          leptonCount += 1.0;
        }
      }

      // muons
      for (int leaf = 0; leaf < branchDict["Muon"]->GetEntries(); leaf++)
      {
        Muon *lepton = (Muon *)branchDict["Muon"]->At(leaf);
        if (lepton->PT > 5.0 && lepton->PT < PTUpperCut && abs(lepton->Eta) < 2.4)
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
          if (elecOverlap(treeReader, branchDict, jet) == -1) // if ==-1, then there is no overlap
          {
            if (muonOverlap(treeReader, branchDict, jet) == -1)
            {
              // its a tau!
              // well constructed taus have pt>20
              if (jet->PT > 20.0 && jet->PT < PTUpperCut && abs(jet->Eta) < 2.4)
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
  }

  //  // cout << endl;
  //  cout << "nLepton Analysis with n = " << PTUpperCut << " done." << endl;
  return histograms;
}

bool inSet(int val, set<int> theSet)
{
  return theSet.count(val) > 0;
}

int ptEtaPhiMjjMt(ExRootTreeReader *treeReader,
                     map<string, TClonesArray *> branchDict,
                     vector<int> &cutsArr,
                     vector<pair<Jet*,Jet*>> &jetPairs,
                     bool (*filter)(ExRootTreeReader *,
                                    map<string, TClonesArray *>,
                                    int,
                                    vector<int> &,
                                    vector<pair<Jet*,Jet*>> &))
{

  int numEvents = 0;

  vector<string> variables = {"pt", "eta", "phi", "Mt"};
  vector<string> particleTypes = {"electron",
                                  "muon",
                                  "tau",
                                  "jet"};

  // create histograms
  map<string, TH1 *> histos;
  for (int i = 0; (unsigned)i < variables.size(); i++)
  {
    for (int j = 0; (unsigned)j < particleTypes.size(); j++)
    {
      int bins = 100;
      float x_min = 0.0;
      float x_max = 15.0;

      if (variables[i].compare("pt") == 0)
      {
        x_max = 500;
        bins = 500;
      }
      else if (variables[i].compare("phi") == 0)
      {
        x_max = TMath::Pi();
        x_min = -TMath::Pi();
      }
      else if (variables[i].compare("Mt") == 0)
      {
        x_max = 1000;
      }
      else
      {
        x_min = -5;
        x_max = 5;
      }

      histos[variables[i] + particleTypes[j]] = blankHistogram(particleTypes[j] + " " + variables[i],
                                                               variables[i] + particleTypes[j],
                                                               bins, x_min, x_max); // check the histogram limits & bins
    }
  }

  histos["mass"] = blankHistogram("Mjj", "Mjj", 100, 0, 10000);
  histos["MET"] = blankHistogram("MET", "MET", 100, 0, 1000);
  histos["DEta"] = blankHistogram("DEta", "DEta", 100, 0, 10);

  Long64_t numberOfEntries = treeReader->GetEntries();

  for (Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    treeReader->ReadEntry(entry);

    if (filter(treeReader, branchDict, entry,cutsArr,jetPairs))
    {
      numEvents += 1;

      set<int> elecIndices;
      set<int> muonIndices;

      // MET stuff
      MissingET *METPointer = (MissingET *)branchDict["MissingET"]->At(0);
      double MET = METPointer->MET;

      // taus
      for (int leaf = 0; leaf < branchDict["Jet"]->GetEntries(); leaf++)
      {
        Jet *jet = (Jet *)branchDict["Jet"]->At(leaf);
        if (jet->TauTag == 1)
        {
          int elecOverlapIndex = elecOverlap(treeReader, branchDict, jet);
          int muonOverlapIndex = muonOverlap(treeReader, branchDict, jet);

          if (elecOverlapIndex == -1)
          {
            if (muonOverlapIndex == -1)
            {
              // its a tau!
              if ((abs(jet->Eta) < 2.4) && (jet->PT > 20.0))
              {
                histos["pttau"]->Fill(jet->PT);
                histos["etatau"]->Fill(jet->Eta);
                histos["phitau"]->Fill(normalizedDphi(jet->Phi));

                double mtval = mt(jet->PT, MET, jet->Phi - METPointer->Phi);
                histos["Mttau"]->Fill(mtval);
              }
            }
            else
            {
              //muon overlap
              muonIndices.insert(muonOverlapIndex);
            }
          }
          else
          {
            //electron overlap
            elecIndices.insert(elecOverlapIndex);
          }
        }
      }

      //      cout << "electrons" << endl;

      // electrons
      for (int leaf = 0; leaf < branchDict["Electron"]->GetEntries(); leaf++)
      {
        Electron *electron = (Electron *)branchDict["Electron"]->At(leaf);
        if (!inSet(leaf, elecIndices))
        {
          if ((abs(electron->Eta) < 2.4) && (electron->PT > 8.0))
          {
            histos["ptelectron"]->Fill(electron->PT);
            histos["etaelectron"]->Fill(electron->Eta);
            histos["phielectron"]->Fill(normalizedDphi(electron->Phi));

            double mtval = mt(electron->PT, MET, electron->Phi - METPointer->Phi);
            histos["Mtelectron"]->Fill(mtval);
          }
        }
      }

      //      cout << "muons" << endl;
      // muons
      for (int leaf = 0; leaf < branchDict["Muon"]->GetEntries(); leaf++)
      {
        Muon *muon = (Muon *)branchDict["Muon"]->At(leaf);
        if (!inSet(leaf, muonIndices))
        {
          if ((abs(muon->Eta) < 2.4) && (muon->PT > 5.0))
          {
            histos["ptmuon"]->Fill(muon->PT);
            histos["etamuon"]->Fill(muon->Eta);
            histos["phimuon"]->Fill(normalizedDphi(muon->Phi));

            double mtval = mt(muon->PT, MET, muon->Phi - METPointer->Phi);
            histos["Mtmuon"]->Fill(mtval);
          }
        }
      }

      //      cout<<"jets"<<endl;
      //jets
      for (int leaf = 0; leaf < branchDict["Jet"]->GetEntries(); leaf++)
      {
        //        cout<<"jet: "<< leaf <<endl;
        Jet *jet = (Jet *)branchDict["Jet"]->At(leaf);
        if (!inSet(leaf, muonIndices) && !inSet(leaf, elecIndices) && !jet->TauTag)
        {

          //          cout<<"calculating overlaps"<<endl;
          int elecOverlapIndex = elecOverlap(treeReader, branchDict, jet);
          int muonOverlapIndex = muonOverlap(treeReader, branchDict, jet);
          //          cout<<"overlaps calculated"<<endl;

          if (elecOverlapIndex == -1 && muonOverlapIndex == -1)
          {
            // reconstruction cuts
            if (jet->PT > 30.0 && abs(jet->Eta) < 5.0)
            {
              //            cout<<"filling pt eta and phi histos"<<endl;
              histos["ptjet"]->Fill(jet->PT);
              histos["etajet"]->Fill(jet->Eta);
              histos["phijet"]->Fill(normalizedDphi(jet->Phi));
              //            cout<<"filling pt eta and phi histos DONE"<<endl;

              //            cout<<"filling mt histo"<<endl;
              //Doesnt make sense but needed to keep code symmetry
              double mtval = mt(jet->PT, MET, jet->Phi - METPointer->Phi);
              histos["Mtjet"]->Fill(mtval);
              //            cout<<"filling mt histo DONE"<<endl;
            }
          }
        }
      }

      // Mjj and deta 
      // min2JetsNotTau guarantees not null pair
      if (min2JetsNotTau(treeReader, branchDict, entry))
      {
        double mass = mjj(jetPairs[entry].first, jetPairs[entry].second);
        histos["mass"]->Fill(mass);

        float deta = deltaEta(jetPairs[entry].first, jetPairs[entry].second);
        histos["DEta"]->Fill(deta);
        
      }

      //MET
      //      cout<<"MET"<<endl;
      histos["MET"]->Fill(MET);
    }
  }

  //  cout << "entry loop done" << endl;

  for (int i = 0; (unsigned)i < variables.size(); i++)
  {
    for (int j = 0; (unsigned)j < particleTypes.size(); j++)
    {
      histos[variables[i] + particleTypes[j]]->Write();
    }
  }

  histos["mass"]->Write();
  histos["MET"]->Write();
  histos["DEta"]->Write();

  return numEvents;
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

void drawLeptonCount(ExRootTreeReader *treeReader,
                     vector<int> ns,
                     map<string, TClonesArray *> branchDict,
                     vector<int> &cutsArr,
                     vector<pair<Jet*,Jet*>> &jetPairs,
                     bool (*filter)(ExRootTreeReader *,
                                    map<string, TClonesArray *>,
                                    int,
                                    vector<int> &,
                                    vector<pair<Jet*,Jet*>> &))
{
  //  cout << "draw lepton count" << endl;
  vector<string> particleTypes = {"lepton", "electron", "muon", "tau"};

  for (int i = 0; (unsigned)i < ns.size(); i++)
  {
    map<string, TH1 *> histoOutput = nLeptonAnalysis(treeReader, ns[i], branchDict,cutsArr, jetPairs, filter);

    for (int j = 0; (unsigned)j < particleTypes.size(); j++)
    {
      // histos[particleTypes[j]].AddLast(histoOutput[particleTypes[j]]);
      histoOutput[particleTypes[j]]->Write();
    }
  }
}

#endif