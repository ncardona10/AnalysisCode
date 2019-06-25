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

bool compareTLVPTDescending(TLorentzVector tlv1, TLorentzVector tlv2);

int main(int argc, char *argv[])
{
  cout<<"Starting phenoanalyzer..."<<endl;

  // standardize print to 2 dp
  cout<<fixed;
  cout<<setprecision(2);

  // Importing Delphes data
  TChain chain("Delphes");
  chain.Add(argv[1]);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);

  vector<string> namesDirectories = { "No_cuts",
                                      "jets_kinematics",
                                      "noBjets",
                                      "METcut",
                                      "VBF_deltaEta",
                                      "Dijet_mass",
                                      "singleLepton_tau",
                                      "singleLepton_muon",
                                      "singleLepton_elec",
                                      "DiLepton_muTau",
                                      "DiLepton_muMu",
                                      "DiLepton_tauTau",
                                      "DiLepton_eTau",
                                      "DiLepton_ee",
                                      "TriLepton_MuMuMu",
                                      "TriLepton_eee"};

  int nDir = namesDirectories.size();


  // output file manager
  TFile *HistoOutputFile = new TFile(argv[2], "RECREATE");  
  TDirectory *theDirectory[nDir];

  // iterate through directory names and create them
  for (int i = 0; (unsigned) i < namesDirectories.size(); i++)
  {
    theDirectory[i] = HistoOutputFile->mkdir(namesDirectories[i].c_str());
  }
  
  cout<< "processing.."<<endl;
  PhenoAnalysis BSM_analysis(treeReader, HistoOutputFile, theDirectory, nDir);
  cout<< "DONE." << endl;
}


PhenoAnalysis::PhenoAnalysis(ExRootTreeReader *treeReader, TFile *theFile, TDirectory *cdDir[], int nDir)
{
  ifstream inFile;
  inFile.open("config.in", ios::in);

  if (!inFile)
  {
    cerr << "ERROR: Can't open input file: " << endl;
    exit(1);
  }

  string inputType = "";


  // configuration file
  cout<< "reading config file..."<<endl;
  
  TEnv *params = new TEnv("config_file");
  params->ReadFile("config.in", kEnvChange);
  // directories names vector
  vector<string> vectKeys = { "b_jet_pt_min",
                              "DR_jet_lep_max",
                              "jet_min_pt",
                              "jet_max_eta",
                              "VBF_jetPt_min",
                              "tau_pt_cut",
                              "tau_pt_cut_max",
                              "tau_eta_cut",
                              "deltaEta_diJet_cut",
                              "diJetmass_cut",
                              "MET_cut",
                              "muon_pt_cut",
                              "muon_pt_cu_max",
                              "muon_eta_cut",
                              "elec_pt_cut",
                              "elec_pt_cut_max",
                              "elec_eta_cut",
                              "muTau_mass_input",
                              "muMu_mass_input",
                              "tauTau_mass_input",
                              "elecTau_mass_input",
                              "elecElec_mass_input"};
  // default values for directories
  vector<double> defaultConfigValues = {30.0,
                                        0.3,
                                        30.0,
                                        5.,
                                        30.0,
                                        20.,
                                        40.,
                                        2.3,
                                        3.8,
                                        500.0,
                                        250.0,
                                        8.,
                                        40.,
                                        2.5,
                                        8.0,
                                        40.0,
                                        8.,
                                        10.,
                                        10.,
                                        10.,
                                        10.,
                                        10.};
  
  map<string, double> configDict;

  for (int i = 0; (unsigned) i < configDict.size(); i++)
  {
    configDict[vectKeys[i]] = params->GetValue(vectKeys[i].c_str(), defaultConfigValues[i]);
  }

  cout<< "Done reading config files." << endl;  

  createHistoMaps(nDir);
  
  Long64_t numberOfEntries = treeReader->GetEntries();

    vector<string> branches = { "Jet",
                                "Electron",
                                "Muon",
                                "MissingET",
                                "Track",
                                "Particle"};


  map<string, TClonesArray*> branchDict;
  // create a dictionary with the branches                           
  for (int i = 0; (unsigned) i < branches.size(); i++)
  {
    TClonesArray *branch = treeReader->UseBranch(branches[i].c_str());
    branchDict[ branches[i]] = branch;
  }


  MissingET *METpointer;

  vector<string> TLVectorNames = { "muTau_TLV",
                                   "tauTau_TLV",
                                   "muMu_TLV",
                                   "elecTau_TLV",
                                   "elecElec_TLV",
                                   "muMuMu_TLV",
                                   "elElEl_TLV"};

  map<string, map<unsigned int, TLorentzVector> > singleParticleTLVDict;
  map<string, map<unsigned int, TLorentzVector> > pairParticleTLVDict;

  // Create TLorentzVectors for single and paired particles
  for (int i = 0; (unsigned) i < TLVectorNames.size(); i++)
  {
    map<unsigned int, TLorentzVector> single_particle_TLV;
    map<unsigned int, TLorentzVector> pairs_particle_TLV;

    singleParticleTLVDict[TLVectorNames[i]] = single_particle_TLV;
    pairParticleTLVDict[TLVectorNames[i]] = pairs_particle_TLV;
  }

  cout<< "Number of entries: " << numberOfEntries << endl;

  for (Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // print percentage of completion
    cout<< "\r" << (100.0*entry)/numberOfEntries<< "%";
    treeReader->ReadEntry(entry);
    int pass_cuts[nDir];
    TLorentzVector jetLeadingVec(0., 0., 0., 0.);
    TLorentzVector jetSleadingVec(0., 0., 0., 0.);

    TLorentzVector Jet_1(0., 0., 0., 0.);
    TLorentzVector Jet_2(0., 0., 0., 0.);

    TLorentzVector jet_i(0., 0., 0., 0.);
    TLorentzVector elec_i(0., 0., 0., 0.);
    TLorentzVector tmp_tlv(0., 0., 0., 0.); //Vector used to swap vectors


    //Tau Vectors
    vector<TLorentzVector> hadronTausTLV;
  
    //Muon Vectors
    vector<TLorentzVector> hadronMuonsTLV;

    //Electron Vectors
    vector<TLorentzVector> hadronElectronsTLV;

    vector<TLorentzVector> hadronTausTLV_gen;

    vector<TLorentzVector> jetsList;

    vector<bool> fill_tau;

    // create and save tlvs in vectors
    for (int i = 0; i < 4; i++)
    {
      // create tlvs
      TLorentzVector hadronTauTLV(0.,0.,0.,0.);
      TLorentzVector hadronMuonTLV(0.,0.,0.,0.);
      TLorentzVector hadronElectronTLV(0.,0.,0.,0.);
      TLorentzVector hadronTauTLV_gen(0.,0.,0.,0.);
      
      // fill the vectors
      hadronTausTLV.push_back(hadronTauTLV);
      hadronMuonsTLV.push_back(hadronMuonTLV);
      hadronElectronsTLV.push_back(hadronElectronTLV);
      hadronTausTLV_gen.push_back(hadronTauTLV_gen);

      fill_tau.push_back(false);
    }
    
    
    // bool fill_tau1 = false;
    // bool fill_tau2 = false;
    // bool fill_tau3 = false;
    // bool fill_tau4 = false;

    METpointer = (MissingET *)branchDict["MissingET"]->At(0);
    double MET = METpointer->MET;
    double Met_phi = METpointer->Phi;
    double tau_transmass = 0.;
    double muon_transmass = 0.;
    double elec_transmass = 0.;
    int ntau_counter = 0;
    int nmuon_counter = 0;
    int nelec_counter = 0;
    double DiJetMass_final = 500.;
    int nBJets = 0;

    //////////////////Tau Channel///////////
    //Search for taus and bjets
    int fillcounter = 0;
    
    for (int j = 0; j < branchDict["Jet"]->GetEntriesFast(); j++)
    {

      Jet *jet = (Jet *)branchDict["Jet"]->At(j);

      if ((jet->BTag == 1) && (abs(jet->Eta) < 2.5) && (jet->PT > configDict["b_jet_pt_min"]))
      {
        nBJets++;
      }

      //Tau search
      if ((jet->TauTag == 1) && (abs(jet->Eta) < 2.5) && (jet->PT > 20.0))
      {
        ntau_counter++;

        double tau_energy = calculateE(jet->Eta, jet->PT, jet->Mass);

        // question: are they the 4 taus with the biggest PT? do we want them to be the biggest?
        if(fillcounter<4){
          hadronTausTLV[fillcounter].SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, tau_energy);
          fillcounter++;
        }
      }
    }

    
    bool tau1_elec_overlap = false;
    bool tau2_elec_overlap = false;
    bool tau3_elec_overlap = false;
    bool tau4_elec_overlap = false;

    bool fill_muon1 = false;
    bool fill_muon2 = false;
    bool fill_muon3 = false;
    bool fill_muon4 = false;

    //Check if taus overlap with muons

    int muo = 0;
    int numOverlaps  = 0;
    int fillMuonCounter = 0;
    TLorentzVector muon_i(0., 0., 0., 0.);

    while(muo < branchDict["Muon"]->GetEntriesFast() && numOverlaps<4)
    {

      Muon *muon = (Muon *)branchDict["Muon"]->At(muo);
      double muon_energy = calculateE(muon->Eta, muon->PT, 0.1056583715);
      muon_i.SetPtEtaPhiE(muon->PT, muon->Eta, muon->Phi, muon_energy);

      for(int i = 0 ; i< 4; i++){
        double DR_tau_muon = hadronTausTLV[i].DeltaR(muon_i);
        if(DR_tau_muon<configDict["DR_jet_lep_max"]){
          hadronTausTLV[i].SetPtEtaPhiE(0., 0., 0., 0.);
          numOverlaps++;
          ntau_counter--;
        }
      }

        // fill muons
      if ((muon->PT > 8.0) && (abs(muon->Eta) < 2.5))
      {
        for(int i = 0 ; i<4; i++){
          // cout<<i << ", " << fillMuonCounter << endl;


          if(fillMuonCounter==0){
            hadronMuonsTLV[fillMuonCounter] = muon_i;
            nmuon_counter++;
            fillMuonCounter++;
          }
          else if(fillMuonCounter==1){
            if (muon_i.DeltaR(hadronMuonsTLV[0]) > 0.3){
              hadronMuonsTLV[fillMuonCounter] = muon_i;
              nmuon_counter++;
              fillMuonCounter++;
            }
          }
          else if(fillMuonCounter==2){
            if ((muon_i.DeltaR(hadronMuonsTLV[1]) > 0.3) && (muon_i.DeltaR(hadronMuonsTLV[0]) > 0.3)){
              hadronMuonsTLV[fillMuonCounter] = muon_i;
              nmuon_counter++;
              fillMuonCounter++;
            }
          }
          else if(fillMuonCounter==3){
            if ((muon_i.DeltaR(hadronMuonsTLV[1]) > 0.3) && (muon_i.DeltaR(hadronMuonsTLV[0]) > 0.3) && (muon_i.DeltaR(hadronMuonsTLV[2]) > 0.3)){
              hadronMuonsTLV[fillMuonCounter] = muon_i;
              nmuon_counter++;
              fillMuonCounter++;
            }
          }

          

        }

       

        // if ((fill_muon1 == false) && (fill_muon2 == false) && (fill_muon3 == false) && (fill_muon4 == false))
        // {
        //   hadronMuonsTLV[0] = muon_i;
        //   fill_muon1 = true;
        //   nmuon_counter++;
        //   continue;
        // }
        // if ((fill_muon1 == true) && (fill_muon2 == false) && (fill_muon3 == false) && (fill_muon4 == false))
        // {
        //   if (muon_i.DeltaR(hadronMuonsTLV[0]) > 0.3)
        //   {
        //     hadronMuonsTLV[1] = muon_i;
        //     fill_muon2 = true;
        //     nmuon_counter++;
        //     continue;
        //   }
        // }
        // if ((fill_muon1 == true) && (fill_muon2 == true) && (fill_muon3 == false) && (fill_muon4 == false))
        // {
        //   if ((muon_i.DeltaR(hadronMuonsTLV[1]) > 0.3) && (muon_i.DeltaR(hadronMuonsTLV[0]) > 0.3))
        //   {
        //     hadronMuonsTLV[2] = muon_i;
        //     fill_muon3 = true;
        //     nmuon_counter++;
        //     continue;
        //   }
        // }
        // if ((fill_muon1 == true) && (fill_muon2 == true) && (fill_muon3 == true) && (fill_muon4 == false))
        // {
        //   if ((muon_i.DeltaR(hadronMuonsTLV[1]) > 0.3) && (muon_i.DeltaR(hadronMuonsTLV[0]) > 0.3) && (muon_i.DeltaR(hadronMuonsTLV[2]) > 0.3))
        //   {
        //     hadronMuonsTLV[3] = muon_i;
        //     fill_muon4 = true;
        //     nmuon_counter++;
        //     continue;
        //   }
        // }
      }


      muo++;
    }


    bool fill_elec1 = false;
    bool fill_elec2 = false;
    bool fill_elec3 = false;
    bool fill_elec4 = false;

    //Check if taus overlap with electrons

    for (int el = 0; el < branchDict["Electron"]->GetEntriesFast(); el++)
    {

      if (tau1_elec_overlap && tau2_elec_overlap && tau3_elec_overlap && tau4_elec_overlap)
      {
        break;
      }

      Electron *elec = (Electron *)branchDict["Electron"]->At(el);
      double elec_energy = calculateE(elec->Eta, elec->PT, 0.000510998902);
      elec_i.SetPtEtaPhiE(elec->PT, elec->Eta, elec->Phi, elec_energy);

      double DR_tau1_elec = hadronTausTLV[0].DeltaR(elec_i);
      double DR_tau2_elec = hadronTausTLV[1].DeltaR(elec_i);
      double DR_tau3_elec = hadronTausTLV[2].DeltaR(elec_i);
      double DR_tau4_elec = hadronTausTLV[3].DeltaR(elec_i);

      if (DR_tau4_elec < configDict["DR_jet_lep_max"])
      {
        tau4_elec_overlap = true;
        hadronTausTLV[3].SetPtEtaPhiE(0., 0., 0., 0.);
        ntau_counter--;
      }
      if (DR_tau3_elec < configDict["DR_jet_lep_max"])
      {
        tau3_elec_overlap = true;
        hadronTausTLV[2].SetPtEtaPhiE(0., 0., 0., 0.);
        ntau_counter--;
      }
      if (DR_tau2_elec < configDict["DR_jet_lep_max"])
      {
        tau2_elec_overlap = true;
        hadronTausTLV[1].SetPtEtaPhiE(0., 0., 0., 0.);
        ntau_counter--;
      }
      if (DR_tau1_elec < configDict["DR_jet_lep_max"])
      {
        tau1_elec_overlap = true;
        hadronTausTLV[0].SetPtEtaPhiE(0., 0., 0., 0.);
        ntau_counter--;
      }

      if ((elec->PT > 8.0) && (abs(elec->Eta) < 2.5))
      {
        if ((fill_elec1 == false) && (fill_elec2 == false) && (fill_elec3 == false) && (fill_elec4 == false))
        {
          hadronElectronsTLV[0] = elec_i;
          fill_elec1 = true;
          nelec_counter++;
          continue;
        }
        if ((fill_elec1 == true) && (fill_elec2 == false) && (fill_elec3 == false) && (fill_elec4 == false))
        {
          if (elec_i.DeltaR(hadronElectronsTLV[0]) > 0.3)
          {
            hadronElectronsTLV[1] = elec_i;
            fill_elec2 = true;
            nelec_counter++;
            continue;
          }
        }
        if ((fill_elec1 == true) && (fill_elec2 == true) && (fill_elec3 == false) && (fill_elec4 == false))
        {
          if ((elec_i.DeltaR(hadronElectronsTLV[1]) > 0.3) && (elec_i.DeltaR(hadronElectronsTLV[0]) > 0.3))
          {
            hadronElectronsTLV[2] = elec_i;
            fill_elec3 = true;
            nelec_counter++;
            continue;
          }
        }
        if ((fill_elec1 == true) && (fill_elec2 == true) && (fill_elec3 == true) && (fill_elec4 == false))
        {
          if ((elec_i.DeltaR(hadronElectronsTLV[2]) > 0.3) && (elec_i.DeltaR(hadronElectronsTLV[1]) > 0.3) && (elec_i.DeltaR(hadronElectronsTLV[0]) > 0.3))
          {
            hadronElectronsTLV[3] = elec_i;
            fill_elec4 = true;
            nelec_counter++;
            continue;
          }
        }
      }
    }

    //Check if jets overlap with taus
    int n_jets = 0;

    for (int l = 0; l < branchDict["Jet"]->GetEntriesFast(); l++)
    {

      Jet *jet = (Jet *)branchDict["Jet"]->At(l);

      if ((jet->PT > configDict["jet_min_pt"]) && (jet->TauTag == 0) && (jet->BTag == 0))
      {

        double jet_i_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
        jet_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_i_energy);

        double DR_tau1_jet = hadronTausTLV[0].DeltaR(jet_i);
        double DR_tau2_jet = hadronTausTLV[1].DeltaR(jet_i);
        double DR_tau3_jet = hadronTausTLV[2].DeltaR(jet_i);
        double DR_tau4_jet = hadronTausTLV[3].DeltaR(jet_i);

        double DR_muon1_jet = hadronMuonsTLV[0].DeltaR(jet_i);
        double DR_muon2_jet = hadronMuonsTLV[1].DeltaR(jet_i);
        double DR_muon3_jet = hadronMuonsTLV[2].DeltaR(jet_i);
        double DR_muon4_jet = hadronMuonsTLV[3].DeltaR(jet_i);

        double DR_elec1_jet = hadronElectronsTLV[0].DeltaR(jet_i);
        double DR_elec2_jet = hadronElectronsTLV[1].DeltaR(jet_i);
        double DR_elec3_jet = hadronElectronsTLV[2].DeltaR(jet_i);
        double DR_elec4_jet = hadronElectronsTLV[3].DeltaR(jet_i);

        if (DR_tau1_jet < configDict["DR_jet_lep_max"])
        {
          hadronTausTLV[0].SetPtEtaPhiE(0., 0., 0., 0.);
          ntau_counter--;
        }
        if (DR_tau2_jet < configDict["DR_jet_lep_max"])
        {
          hadronTausTLV[1].SetPtEtaPhiE(0., 0., 0., 0.);
          ntau_counter--;
        }
        if (DR_tau3_jet < configDict["DR_jet_lep_max"])
        {
          hadronTausTLV[2].SetPtEtaPhiE(0., 0., 0., 0.);
          ntau_counter--;
        }
        if (DR_tau4_jet < configDict["DR_jet_lep_max"])
        {
          hadronTausTLV[3].SetPtEtaPhiE(0., 0., 0., 0.);
          ntau_counter--;
        }

        // remove overlaps of jets with muons
        if (DR_muon1_jet < configDict["DR_jet_lep_max"])
        {
          hadronMuonsTLV[0].SetPtEtaPhiE(0., 0., 0., 0.);
          nmuon_counter--;
        }
        if (DR_muon2_jet < configDict["DR_jet_lep_max"])
        {
          hadronMuonsTLV[1].SetPtEtaPhiE(0., 0., 0., 0.);
          nmuon_counter--;
        }
        if (DR_muon3_jet < configDict["DR_jet_lep_max"])
        {
          hadronMuonsTLV[2].SetPtEtaPhiE(0., 0., 0., 0.);
          nmuon_counter--;
        }
        if (DR_muon4_jet < configDict["DR_jet_lep_max"])
        {
          hadronMuonsTLV[3].SetPtEtaPhiE(0., 0., 0., 0.);
          nmuon_counter--;
        }

        // remove overlaps of jets with electrons
        if (DR_elec1_jet < configDict["DR_jet_lep_max"])
        {
          hadronElectronsTLV[0].SetPtEtaPhiE(0., 0., 0., 0.);
          nelec_counter--;
        }
        if (DR_elec2_jet < configDict["DR_jet_lep_max"])
        {
          hadronElectronsTLV[1].SetPtEtaPhiE(0., 0., 0., 0.);
          nelec_counter--;
        }
        if (DR_elec3_jet < configDict["DR_jet_lep_max"])
        {
          hadronElectronsTLV[2].SetPtEtaPhiE(0., 0., 0., 0.);
          nelec_counter--;
        }
        if (DR_elec4_jet < configDict["DR_jet_lep_max"])
        {
          hadronElectronsTLV[3].SetPtEtaPhiE(0., 0., 0., 0.);
          nelec_counter--;
        }

        if ((n_jets <= 6) && (jet->TauTag == 0) && (jet->BTag == 0) && (abs(jet->Eta) < 5.0))
        {
          jetsList.push_back(jet_i);
          n_jets++;
        }
      }
    }

    //Order taus by pt - descending order
    sort(hadronTausTLV.begin(), hadronTausTLV.end(), compareTLVPTDescending); 

    //Order muons by pt- descending order
    sort(hadronMuonsTLV.begin(), hadronMuonsTLV.end(), compareTLVPTDescending); 

    //Order electrons by pt- descending order
    sort(hadronElectronsTLV.begin(), hadronElectronsTLV.end(), compareTLVPTDescending); 

    // GEN INFO
    bool fill_tau1_gen = false;
    bool fill_tau2_gen = false;
    bool fill_tau3_gen = false;

    for (int tau_i = 0; tau_i < branchDict["Particle"]->GetEntriesFast(); tau_i++)
    {
      // Look for the particle PDGID
      GenParticle *tau_gen = (GenParticle *)branchDict["Particle"]->At(tau_i);
      if (abs(tau_gen->PID) == 15)
      {
        double tau_energy = calculateE(tau_gen->Eta, tau_gen->PT, tau_gen->Mass);
        hadronTausTLV_gen[0].SetPtEtaPhiE(tau_gen->PT, tau_gen->Eta, tau_gen->Phi, tau_energy);
        if (fill_tau1_gen == false)
        {
          hadronTausTLV_gen[1] = hadronTausTLV_gen[0];
          fill_tau1_gen = true;
          continue;
        }
        if ((fill_tau1_gen == true) && (fill_tau2_gen == false))
        {
          double delta_R_gen1 = hadronTausTLV_gen[1].DeltaR(hadronTausTLV_gen[0]);
          if (delta_R_gen1 > 0.3)
          {
            hadronTausTLV_gen[2] = hadronTausTLV_gen[0];
            fill_tau2_gen = true;
          }
          continue;
        }
        if ((fill_tau2_gen == true) && (fill_tau3_gen == false))
        {
          double delta_R_gen2 = hadronTausTLV_gen[2].DeltaR(hadronTausTLV_gen[0]);
          if (delta_R_gen2 > 0.3)
          {
            hadronTausTLV_gen[3] = hadronTausTLV_gen[0];
            fill_tau3_gen = true;
          }
          continue;
        }
      }
    }

    double tau1_track_DR_min = 999.;
    double tau1_track_DR;
    //Search tau track
    for (Int_t i = 0; i < branchDict["Track"]->GetEntriesFast(); i++)
    {
      Track *track = (Track *)branchDict["Track"]->At(i);
      tau1_track_DR = calculate_deltaR(hadronTausTLV[0], track);
      if (tau1_track_DR < tau1_track_DR_min)
      {
        tau1_track_DR_min = tau1_track_DR;
      }
    }

    int dijet_index1 = 0;
    int dijet_index2 = 0;

    //Search DiJetMass
    for (UInt_t k = 0; k < jetsList.size(); k++)
    {

      if (jetsList.size() < 2)
      {
        break;
      }
      Jet_1 = jetsList[k];

      if ((Jet_1.Pt() < configDict["VBF_jetPt_min"]) || (abs(Jet_1.Eta()) > 5.0))
      {
        continue;
      }

      for (UInt_t sj = k + 1; sj < jetsList.size(); sj++)
      {
        if (sj != k)
        {
          Jet_2 = jetsList[sj];
          if ((Jet_2.Pt() < configDict["VBF_jetPt_min"]) || (abs(Jet_2.Eta()) > 5.0))
          {
            continue;
          }
          double DiJetMass = (Jet_1 + Jet_2).M();
          if (DiJetMass > DiJetMass_final)
          {
            DiJetMass_final = DiJetMass;
            jetLeadingVec = Jet_1;
            jetSleadingVec = Jet_2;
            dijet_index1 = k;
            dijet_index2 = sj;
          }
        }
      }
    }

    if (jetLeadingVec.Pt() < jetSleadingVec.Pt())
    {
      tmp_tlv = jetLeadingVec;
      jetLeadingVec = jetSleadingVec;
      jetSleadingVec = tmp_tlv;
    }

    if (jetsList.size() > 2)
    {

      //Remove the diJet pair from the jets list
      tmp_tlv = jetsList[jetsList.size() - 1];
      jetsList[jetsList.size() - 1] = jetsList[dijet_index1];
      jetsList[dijet_index1] = tmp_tlv;

      tmp_tlv = jetsList[jetsList.size() - 2];
      jetsList[jetsList.size() - 2] = jetsList[dijet_index2];
      jetsList[dijet_index2] = tmp_tlv;

      jetsList.pop_back();
      jetsList.pop_back();
    }

    //Check for jets pt condition
    int jet_pt_condition = 0;
    for (Int_t i = 0; i < (Int_t)jetsList.size(); i++)
    {
      if ((jetsList[i].Pt() > configDict["jet_min_pt"]) && (abs(jetsList[i].Eta()) < 5.0))
      {
        jet_pt_condition++;
      }
    }

    double delta_eta_diJet = abs(jetLeadingVec.Eta() - jetSleadingVec.Eta());
    double tauMass = (hadronTausTLV[0] + hadronTausTLV[1]).M();
    tau_transmass = TMath::Sqrt(TMath::Abs(2 * hadronTausTLV[0].Pt() * MET * (1 - TMath::Cos(normalizedDphi(hadronTausTLV[0].Phi() - Met_phi)))));
    muon_transmass = TMath::Sqrt(TMath::Abs(2 * hadronMuonsTLV[0].Pt() * MET * (1 - TMath::Cos(normalizedDphi(hadronMuonsTLV[0].Phi() - Met_phi)))));
    elec_transmass = TMath::Sqrt(TMath::Abs(2 * hadronElectronsTLV[0].Pt() * MET * (1 - TMath::Cos(normalizedDphi(hadronElectronsTLV[0].Phi() - Met_phi)))));
    double ht = 0.;
    double st = 0.;

    ht += jetLeadingVec.Pt() + jetSleadingVec.Pt();
    for (Int_t i = 0; i < (Int_t)jetsList.size(); i++)
    {
      ht += jetsList[i].Pt();
    }

    st += ht + hadronTausTLV[0].Pt() + hadronTausTLV[1].Pt() + hadronTausTLV[2].Pt() + jetLeadingVec.Pt() + jetSleadingVec.Pt() + MET;

    singleParticleTLVDict["muTau_TLV"][0] = hadronTausTLV[0];
    singleParticleTLVDict["muTau_TLV"][1] = hadronTausTLV[1];
    singleParticleTLVDict["muTau_TLV"][2] = hadronTausTLV[2];
    singleParticleTLVDict["muTau_TLV"][3] = hadronTausTLV[3];
    singleParticleTLVDict["muTau_TLV"][4] = hadronMuonsTLV[0];
    singleParticleTLVDict["muTau_TLV"][5] = hadronMuonsTLV[1];
    singleParticleTLVDict["muTau_TLV"][6] = hadronMuonsTLV[2];
    singleParticleTLVDict["muTau_TLV"][7] = hadronMuonsTLV[3];

    double muTau_mass_i = 0.;

    for (Int_t i = 0; i < 4; i++)
    {
      if ((singleParticleTLVDict["muTau_TLV"][i].Pt() > configDict["tau_pt_cut"]) && (abs(singleParticleTLVDict["muTau_TLV"][i].Eta()) < configDict["tau_eta_cut"]))
      {
        for (int j = 4; j < 8; j++)
        {
          if ((singleParticleTLVDict["muTau_TLV"][j].Pt() > configDict["muon_pt_cut"]) && (abs(singleParticleTLVDict["muTau_TLV"][j].Eta() < configDict["muon_eta_cut"])))
          {
            double muTau_mass = (singleParticleTLVDict["muTau_TLV"][i] + singleParticleTLVDict["muTau_TLV"][j]).M();
            if (muTau_mass > muTau_mass_i)
            {
              muTau_mass_i = muTau_mass;
              pairParticleTLVDict["muTau_TLV"][0] = singleParticleTLVDict["muTau_TLV"][i];
              pairParticleTLVDict["muTau_TLV"][1] = singleParticleTLVDict["muTau_TLV"][j];
            }
          }
        }
      }
    }

    singleParticleTLVDict["muMu_TLV"][0] = hadronMuonsTLV[0];
    singleParticleTLVDict["muMu_TLV"][1] = hadronMuonsTLV[1];
    singleParticleTLVDict["muMu_TLV"][2] = hadronMuonsTLV[2];
    singleParticleTLVDict["muMu_TLV"][3] = hadronMuonsTLV[3];

    double muMu_mass_i = 0.;

    for (Int_t i = 0; i < 4; i++)
    {
      if ((singleParticleTLVDict["muMu_TLV"][i].Pt() > configDict["muon_pt_cut"]) && (abs(singleParticleTLVDict["muMu_TLV"][i].Eta()) < configDict["muon_eta_cut"]))
      {
        for (int j = i + 1; j < 4; j++)
        {
          if ((singleParticleTLVDict["muMu_TLV"][j].Pt() > configDict["muon_pt_cut"]) && (abs(singleParticleTLVDict["muMu_TLV"][j].Eta()) < configDict["muon_eta_cut"]) && (i != j))
          {
            double muMu_mass = (singleParticleTLVDict["muMu_TLV"][i] + singleParticleTLVDict["muMu_TLV"][j]).M();
            if (muMu_mass > muMu_mass_i)
            {
              muMu_mass_i = muMu_mass;
              pairParticleTLVDict["muMu_TLV"][0] = singleParticleTLVDict["muMu_TLV"][i];
              pairParticleTLVDict["muMu_TLV"][1] = singleParticleTLVDict["muMu_TLV"][j];
            }
          }
        }
      }
    }

    singleParticleTLVDict["tauTau_TLV"][0] = hadronTausTLV[0];
    singleParticleTLVDict["tauTau_TLV"][1] = hadronTausTLV[1];
    singleParticleTLVDict["tauTau_TLV"][2] = hadronTausTLV[2];
    singleParticleTLVDict["tauTau_TLV"][3] = hadronTausTLV[3];

    double tauTau_mass_i = 0.;

    for (Int_t i = 0; i < 4; i++)
    {
      if ((singleParticleTLVDict["tauTau_TLV"][i].Pt() > configDict["tau_pt_cut"]) && (abs(singleParticleTLVDict["tauTau_TLV"][i].Eta()) < configDict["tau_eta_cut"]))
      {
        for (int j = i + 1; j < 4; j++)
        {
          if ((singleParticleTLVDict["tauTau_TLV"][j].Pt() > configDict["tau_pt_cut"]) && (abs(singleParticleTLVDict["tauTau_TLV"][j].Eta()) < configDict["tau_eta_cut"]) && (i != j))
          {
            double tauTau_mass = (singleParticleTLVDict["tauTau_TLV"][i] + singleParticleTLVDict["tauTau_TLV"][j]).M();
            if (tauTau_mass > tauTau_mass_i)
            {
              tauTau_mass_i = tauTau_mass;
              pairParticleTLVDict["tauTau_TLV"][0] = singleParticleTLVDict["tauTau_TLV"][i];
              pairParticleTLVDict["tauTau_TLV"][1] = singleParticleTLVDict["tauTau_TLV"][j];
            }
          }
        }
      }
    }

    // elec - tau pairs
    singleParticleTLVDict["elecTau_TLV"][0] = hadronElectronsTLV[0];
    singleParticleTLVDict["elecTau_TLV"][1] = hadronElectronsTLV[1];
    singleParticleTLVDict["elecTau_TLV"][2] = hadronElectronsTLV[2];
    singleParticleTLVDict["elecTau_TLV"][3] = hadronElectronsTLV[3];
    singleParticleTLVDict["elecTau_TLV"][4] = hadronTausTLV[0];
    singleParticleTLVDict["elecTau_TLV"][5] = hadronTausTLV[1];
    singleParticleTLVDict["elecTau_TLV"][6] = hadronTausTLV[2];
    singleParticleTLVDict["elecTau_TLV"][7] = hadronTausTLV[3];

    double elecTau_mass_i = 0.;

    for (Int_t i = 0; i < 4; i++)
    {
      if ((singleParticleTLVDict["elecTau_TLV"][i].Pt() > configDict["elec_pt_cut"]) && (abs(singleParticleTLVDict["elecTau_TLV"][i].Eta()) < configDict["elec_eta_cut"]))
      {
        for (int j = 4; j < 8; j++)
        {
          if ((singleParticleTLVDict["elecTau_TLV"][j].Pt() > configDict["tau_pt_cut"]) && (abs(singleParticleTLVDict["elecTau_TLV"][j].Eta()) < configDict["tau_eta_cut"]) && (i != j))
          {
            double elecTau_mass = (singleParticleTLVDict["elecTau_TLV"][i] + singleParticleTLVDict["elecTau_TLV"][j]).M();
            if (elecTau_mass > elecTau_mass_i)
            {
              elecTau_mass_i = elecTau_mass;
              pairParticleTLVDict["elecTau_TLV"][0] = singleParticleTLVDict["elecTau_TLV"][i];
              pairParticleTLVDict["elecTau_TLV"][1] = singleParticleTLVDict["elecTau_TLV"][j];
            }
          }
        }
      }
    }

    // elec - elec pairs
    singleParticleTLVDict["elecElec_TLV"][0] = hadronElectronsTLV[0];
    singleParticleTLVDict["elecElec_TLV"][1] = hadronElectronsTLV[1];
    singleParticleTLVDict["elecElec_TLV"][2] = hadronElectronsTLV[2];
    singleParticleTLVDict["elecElec_TLV"][3] = hadronElectronsTLV[3];

    double elecElec_mass_i = 0.;

    for (Int_t i = 0; i < 4; i++)
    {
      if ((singleParticleTLVDict["elecElec_TLV"][i].Pt() > configDict["elec_pt_cut"]) && (abs(singleParticleTLVDict["elecElec_TLV"][i].Eta()) < configDict["elec_eta_cut"]))
      {
        for (int j = i + 1; j < 4; j++)
        {
          if ((singleParticleTLVDict["elecElec_TLV"][j].Pt() > configDict["elec_pt_cut"]) && (abs(singleParticleTLVDict["elecElec_TLV"][j].Eta()) < configDict["elec_eta_cut"]) && (i != j))
          {
            double elecElec_mass = (singleParticleTLVDict["elecElec_TLV"][i] + singleParticleTLVDict["elecElec_TLV"][j]).M();
            if (elecElec_mass > elecElec_mass_i)
            {
              elecElec_mass_i = elecElec_mass;
              pairParticleTLVDict["elecElec_TLV"][0] = singleParticleTLVDict["elecElec_TLV"][i];
              pairParticleTLVDict["elecElec_TLV"][1] = singleParticleTLVDict["elecElec_TLV"][j];
            }
          }
        }
      }
    }

    // mu-mu-mu trio
    singleParticleTLVDict["muMuMu_TLV"][0] = hadronMuonsTLV[0];
    singleParticleTLVDict["muMuMu_TLV"][1] = hadronMuonsTLV[1];
    singleParticleTLVDict["muMuMu_TLV"][2] = hadronMuonsTLV[2];
    singleParticleTLVDict["muMuMu_TLV"][3] = hadronMuonsTLV[3];

    double muMuMu_mass_i = 0.;

    for (Int_t i = 0; i < 4; i++)
    {
      if ((singleParticleTLVDict["muMuMu_TLV"][i].Pt() > configDict["muon_pt_cut"]) && (abs(singleParticleTLVDict["muMuMu_TLV"][i].Eta()) < configDict["muon_eta_cut"]))
      {
        for (int j = i + 1; j < 4; j++)
        {
          if ((singleParticleTLVDict["muMuMu_TLV"][j].Pt() > configDict["muon_pt_cut"]) && (abs(singleParticleTLVDict["muMuMu_TLV"][j].Eta()) < configDict["muon_eta_cut"]) && (i != j))
          {
            for (int k = j + 1; k < 4; k++)
            {
              double muMuMu_mass = (singleParticleTLVDict["muMuMu_TLV"][i] + singleParticleTLVDict["muMuMu_TLV"][j] + singleParticleTLVDict["muMuMu_TLV"][k]).M();
              if (muMuMu_mass > muMuMu_mass_i)
              {
                muMuMu_mass_i = muMuMu_mass;
                pairParticleTLVDict["muMuMu_TLV"][0] = singleParticleTLVDict["muMuMu_TLV"][i];
                pairParticleTLVDict["muMuMu_TLV"][1] = singleParticleTLVDict["muMuMu_TLV"][j];
                pairParticleTLVDict["muMuMu_TLV"][2] = singleParticleTLVDict["muMuMu_TLV"][k];
              }
            }
          }
        }
      }
    }

    // elec-elec-elec trio
    singleParticleTLVDict["elElEL_TLV"][0] = hadronElectronsTLV[0];
    singleParticleTLVDict["elElEL_TLV"][1] = hadronElectronsTLV[1];
    singleParticleTLVDict["elElEL_TLV"][2] = hadronElectronsTLV[2];
    singleParticleTLVDict["elElEL_TLV"][3] = hadronElectronsTLV[3];

    double elElEl_mass_i = 0.;

    for (Int_t i = 0; i < 4; i++)
    {
      if ((singleParticleTLVDict["elElEl_TLV"][i].Pt() > configDict["elec_pt_cut"]) && (abs(singleParticleTLVDict["elElEl_TLV"][i].Eta()) < configDict["elec_eta_cut"]))
      {
        for (int j = i + 1; j < 4; j++)
        {
          if ((singleParticleTLVDict["elElEl_TLV"][j].Pt() > configDict["elec_pt_cut"]) && (abs(singleParticleTLVDict["elElEl_TLV"][j].Eta()) < configDict["elec_eta_cut"]) && (i != j))
          {
            for (int k = j + 1; k < 4; k++)
            {
              double elElEl_mass = (singleParticleTLVDict["elElEl_TLV"][i] + singleParticleTLVDict["elElEL_TLV"][j] + singleParticleTLVDict["elElEL_TLV"][k]).M();
              if (elElEl_mass > elElEl_mass_i)
              {
                elElEl_mass_i = elElEl_mass;
                pairParticleTLVDict["elElEl_TLV"][0] = singleParticleTLVDict["elElEL_TLV"][i];
                pairParticleTLVDict["elElEl_TLV"][1] = singleParticleTLVDict["elElEL_TLV"][j];
                pairParticleTLVDict["elElEl_TLV"][2] = singleParticleTLVDict["elElEL_TLV"][k];
              }
            }
          }
        }
      }
    }

    //////////////////////////////////////////////////////////////////////////////////

    //////// Apply cuts /////////
    // Events with no cuts
    pass_cuts[0] = 1;

    // Number of bjets
    if ((nBJets == 0))
    {
      pass_cuts[1] = 1;
    }
    // dijet kinematics
    if ((pass_cuts[1] == 1) && (jetLeadingVec.Pt() > configDict["jet_min_pt"]) && (jetSleadingVec.Pt() > configDict["jet_min_pt"]) && (abs(jetLeadingVec.Eta()) < configDict["jet_max_eta"
]) && (abs(jetSleadingVec.Eta()) < configDict["jet_max_eta"
]))
    {
      pass_cuts[2] = 1;
    }
    // MET cut
    if ((pass_cuts[2] == 1) && (MET > configDict["MET_cut"]))
    {
      pass_cuts[3] = 1;
    }
    // deltaEta cut
    if ((pass_cuts[3] == 1) && (delta_eta_diJet > configDict["deltaEta_diJet_cut"]))
    {
      pass_cuts[4] = 1;
    }
    // diJet mass cut
    if ((pass_cuts[4] == 1) && (DiJetMass_final > configDict["diJetmass_cut"]))
    {
      pass_cuts[5] = 1;
    }
    // Events with 1 tau with kinematic cuts
    if ((pass_cuts[5] == 1) && (ntau_counter == 1) && (nmuon_counter == 0) && (nelec_counter == 0) && (hadronTausTLV[0].Pt() > configDict["tau_pt_cut"]) && (hadronTausTLV[0].Pt() < configDict["tau_pt_cut_max"
]) && (abs(hadronTausTLV[0].Eta()) < configDict["tau_eta_cut"]))
    {
      pass_cuts[6] = 1;
    }
    // Events with 1 muon  with kinematic cuts
    if ((pass_cuts[5] == 1) && (nmuon_counter == 1) && (ntau_counter == 0) && (nelec_counter == 0) && (hadronMuonsTLV[0].Pt() > configDict["muon_pt_cut"]) && (hadronMuonsTLV[0].Pt() < configDict["muon_pt_cut_max"]) && (abs(hadronMuonsTLV[0].Eta()) < configDict["muon_eta_cut"]))
    {
      pass_cuts[7] = 1;
    }
    // Events with 1 electron with kinematic cuts
    if ((pass_cuts[5] == 1) && (nelec_counter == 1) && (ntau_counter == 0) && (nmuon_counter == 0) && (hadronElectronsTLV[0].Pt() > configDict["elec_pt_cut"]) && (hadronElectronsTLV[0].Pt() < configDict["elec_pt_cut_max"]) && (abs(hadronElectronsTLV[0].Eta()) < configDict["elec_eta_cut"]))
    {
      pass_cuts[8] = 1;
    }
    // // muTau kinematics
    // if ((pass_cuts[5] == 1) && (nmuon_counter == 1) && (ntau_counter == 1) && (nelec_counter == 0) && (pairParticleTLVDict["muTau_TLV"][0].Pt() > configDict["tau_pt_cut"]) && (pairParticleTLVDict["muTau_TLV"][1].Pt() > configDict["muon_pt_cut"]) && (abs(pairParticleTLVDict["muTau_TLV"][0].Eta()) < configDict["tau_eta_cut"]) && (abs(pairParticleTLVDict["muTau_TLV"][1].Eta()) < configDict["muon_eta_cut"]) && (muTau_mass_i > configDict["muTau_mass_input"]))
    // {
    //   pass_cuts[9] = 1;
    // }
    // // muMu kinematics
    // if ((pass_cuts[5] == 1) && (nmuon_counter == 2) && (ntau_counter == 0) && (nelec_counter == 0) && (pairParticleTLVDict["muMu_TLV"][0].Pt() > configDict["muon_pt_cut"]) && (pairParticleTLVDict["muMu_TLV"][1].Pt() > configDict["muon_pt_cut"]) && (abs(pairParticleTLVDict["muMu_TLV"][0].Eta()) < configDict["muon_eta_cut"]) && (abs(pairParticleTLVDict["muMu_TLV"][1].Eta()) < configDict["muon_eta_cut"]) && (muMu_mass_i > configDict["muMu_mass_input"]))
    // {
    //   pass_cuts[10] = 1;
    // }
    // // tauTau kinematics
    // if ((pass_cuts[5] == 1) && (ntau_counter == 2) && (nmuon_counter == 0) && (nelec_counter == 0) && (pairParticleTLVDict["tauTau_TLV"][0].Pt() > configDict["tau_pt_cut"]) && (pairParticleTLVDict["tauTau_TLV"][1].Pt() > configDict["tau_pt_cut"]) && (abs(pairParticleTLVDict["tauTau_TLV"][0].Eta()) < configDict["tau_eta_cut"]) && (abs(pairParticleTLVDict["tauTau_TLV"][1].Eta()) < configDict["tau_eta_cut"]) && (tauTau_mass_i > configDict["tauTau_mass_input"]))
    // {
    //   pass_cuts[11] = 1;
    // }
    // // elecTau kinematics
    // if ((pass_cuts[5] == 1) && (nelec_counter == 1) && (ntau_counter == 1) && (nmuon_counter == 0) && (pairParticleTLVDict["elecTau_TLV"][0].Pt() > configDict["elec_pt_cut"]) && (pairParticleTLVDict["elecTau_TLV"][1].Pt() > configDict["tau_pt_cut"]) && (abs(pairParticleTLVDict["elecTau_TLV"][0].Eta()) < configDict["elec_eta_cut"]) && (abs(pairParticleTLVDict["elecTau_TLV"][1].Eta()) < configDict["tau_eta_cut"]) && (elecTau_mass_i > configDict["elecTau_mass_input"]))
    // {
    //   pass_cuts[12] = 1;
    // }
    // // elecElec kinematics
    // if ((pass_cuts[5] == 1) && (nelec_counter == 2) && (ntau_counter == 0) && (nmuon_counter == 0) && (pairParticleTLVDict["elecElec_TLV"][0].Pt() > configDict["elec_pt_cut"]) && (pairParticleTLVDict["elecElec_TLV"][1].Pt() > configDict["elec_pt_cut"]) && (abs(pairParticleTLVDict["elecElec_TLV"][0].Eta()) < configDict["elec_eta_cut"]) && (abs(pairParticleTLVDict["elecElec_TLV"][1].Eta()) < configDict["elec_eta_cut"]) && (elecElec_mass_i > configDict["elecElec_mass_input"]))
    {
      pass_cuts[13] = 1;
    }
    // muMuMu kinematics
    if ((pass_cuts[5] == 1) && (nmuon_counter == 3) && (ntau_counter == 0) && (nelec_counter == 0) && (pairParticleTLVDict["muMuMu_TLV"][0].Pt() > configDict["muon_pt_cut"]) && (pairParticleTLVDict["muMuMu_TLV"][1].Pt() > configDict["muon_pt_cut"]) && (pairParticleTLVDict["muMuMu_TLV"][2].Pt() > configDict["muon_pt_cut"]) && (abs(pairParticleTLVDict["muMuMu_TLV"][0].Eta()) < configDict["muon_eta_cut"]) && (abs(pairParticleTLVDict["muMuMu_TLV"][1].Eta()) < configDict["muon_eta_cut"]) && (abs(pairParticleTLVDict["muMuMu_TLV"][2].Eta()) < configDict["muon_eta_cut"]))
    {
      pass_cuts[14] = 1;
    }
    // elElEl kinematics
    if ((pass_cuts[5] == 1) && (nelec_counter == 3) && (ntau_counter == 0) && (nmuon_counter == 0) && (pairParticleTLVDict["elElEl_TLV"][0].Pt() > configDict["elec_pt_cut"]) && (pairParticleTLVDict["elElEl_TLV"][1].Pt() > configDict["elec_pt_cut"]) && (pairParticleTLVDict["elElEl_TLV"][2].Pt() > configDict["elec_pt_cut"]) && (abs(pairParticleTLVDict["elElEl_TLV"][0].Eta()) < configDict["elec_eta_cut"]) && (abs(pairParticleTLVDict["elElEl_TLV"][1].Eta()) < configDict["elec_eta_cut"]) && (abs(pairParticleTLVDict["elElEl_TLV"][2].Eta()) < configDict["elec_eta_cut"]))
    {
      pass_cuts[15] = 1;
    }

    //Fill histograms
    for (Int_t i = 0; i < nDir; i++)
    {

      _hmap_Nevents[i]->Fill(0.0);

      if (pass_cuts[i] == 1)
      {

        _hmap_n_tau[i]->Fill(ntau_counter);
        _hmap_n_jets[i]->Fill(n_jets);
        _hmap_Nevents[i]->Fill(1.0);

        if (jetLeadingVec.Pt() > 1.0)
        {
          _hmap_lead_jet_pT[i]->Fill(jetLeadingVec.Pt());
          _hmap_lead_jet_eta[i]->Fill(jetLeadingVec.Eta());
          _hmap_lead_jet_phi[i]->Fill(jetLeadingVec.Phi());
          _hmap_slead_jet_pT[i]->Fill(jetSleadingVec.Pt());
          _hmap_slead_jet_eta[i]->Fill(jetSleadingVec.Eta());
          _hmap_slead_jet_phi[i]->Fill(jetSleadingVec.Phi());
        }
        if (hadronTausTLV[0].Pt() > 1.0)
        {
          _hmap_tau1_pT[i]->Fill(hadronTausTLV[0].Pt());
          _hmap_tau1_eta[i]->Fill(hadronTausTLV[0].Eta());
          _hmap_tau1_phi[i]->Fill(hadronTausTLV[0].Phi());
        }
        if (hadronTausTLV[1].Pt() > 1.0)
        {
          _hmap_tau2_pT[i]->Fill(hadronTausTLV[1].Pt());
          _hmap_tau2_eta[i]->Fill(hadronTausTLV[1].Eta());
          _hmap_tau2_phi[i]->Fill(hadronTausTLV[1].Phi());
          _hmap_Delta_pT[i]->Fill(abs(hadronTausTLV[0].Pt() - hadronTausTLV[1].Pt()));
        }
        if (hadronMuonsTLV[0].Pt() > 1.0)
        {
          _hmap_muon1_pT[i]->Fill(hadronMuonsTLV[0].Pt());
          _hmap_muon1_eta[i]->Fill(hadronMuonsTLV[0].Eta());
          _hmap_muon1_phi[i]->Fill(hadronMuonsTLV[0].Phi());
        }
        if (hadronMuonsTLV[1].Pt() > 1.0)
        {
          _hmap_muon2_pT[i]->Fill(hadronMuonsTLV[1].Pt());
          _hmap_muon2_eta[i]->Fill(hadronMuonsTLV[1].Eta());
          _hmap_muon2_phi[i]->Fill(hadronMuonsTLV[1].Phi());
        }
        if (hadronElectronsTLV[0].Pt() > 1.0)
        {
          _hmap_elec1_pT[i]->Fill(hadronElectronsTLV[0].Pt());
          _hmap_elec1_eta[i]->Fill(hadronElectronsTLV[0].Eta());
          _hmap_elec1_phi[i]->Fill(hadronElectronsTLV[0].Phi());
        }
        if (hadronElectronsTLV[1].Pt() > 1.0)
        {
          _hmap_elec2_pT[i]->Fill(hadronElectronsTLV[1].Pt());
          _hmap_elec2_eta[i]->Fill(hadronElectronsTLV[1].Eta());
          _hmap_elec2_phi[i]->Fill(hadronElectronsTLV[1].Phi());
        }

        if (tauMass > 0)
        {
          _hmap_tauMass[i]->Fill(tauMass);
        }
        if (MET > 0)
        {
          _hmap_MET[i]->Fill(MET);
        }
        if (tau_transmass > 0)
        {
          _hmap_tau_transmass[i]->Fill(tau_transmass);
        }
        if (muon_transmass > 0)
        {
          _hmap_muon_transmass[i]->Fill(muon_transmass);
        }
        if (elec_transmass > 0)
        {
          _hmap_elec_transmass[i]->Fill(elec_transmass);
        }
        if (ht > 0.0)
        {
          _hmap_ht[i]->Fill(ht);
        }
        if (st > 0.0)
        {
          _hmap_st[i]->Fill(st);
        }
        if (DiJetMass_final > 100.0)
        {
          _hmap_dijet_mass[i]->Fill(DiJetMass_final);
        }
        if (delta_eta_diJet > 0)
        {
          _hmap_dijet_deltaEta[i]->Fill(delta_eta_diJet);
        }
      }
    }
  } // end entry loop for tau channel

  theFile->cd();
  for (int d = 0; d < nDir; d++)
  {
    cdDir[d]->cd();
    _hmap_Nevents[d]->Write();
    _hmap_n_jets[d]->Write();
    _hmap_n_tau[d]->Write();
    _hmap_lead_jet_pT[d]->Write();
    _hmap_lead_jet_eta[d]->Write();
    _hmap_lead_jet_phi[d]->Write();
    _hmap_slead_jet_pT[d]->Write();
    _hmap_slead_jet_eta[d]->Write();
    _hmap_slead_jet_phi[d]->Write();
    _hmap_tau1_pT[d]->Write();
    _hmap_tau1_eta[d]->Write();
    _hmap_tau1_phi[d]->Write();
    _hmap_tau2_pT[d]->Write();
    _hmap_tau2_eta[d]->Write();
    _hmap_tau2_phi[d]->Write();
    _hmap_muon1_pT[d]->Write();
    _hmap_muon1_eta[d]->Write();
    _hmap_muon1_phi[d]->Write();
    _hmap_muon2_pT[d]->Write();
    _hmap_muon2_eta[d]->Write();
    _hmap_muon2_phi[d]->Write();
    _hmap_elec1_pT[d]->Write();
    _hmap_elec1_eta[d]->Write();
    _hmap_elec1_phi[d]->Write();
    _hmap_elec2_pT[d]->Write();
    _hmap_elec2_eta[d]->Write();
    _hmap_elec2_phi[d]->Write();
    _hmap_Delta_pT[d]->Write();
    //Gen/////
    _hmap_Gentau1_pT[d]->Write();
    _hmap_Gentau1_eta[d]->Write();
    _hmap_Gentau1_phi[d]->Write();
    _hmap_Gentau2_pT[d]->Write();
    _hmap_Gentau2_eta[d]->Write();
    _hmap_Gentau2_phi[d]->Write();
    _hmap_st[d]->Write();
    _hmap_GenDelta_pT[d]->Write();
    ////////
    _hmap_tauMass[d]->Write();
    _hmap_MET[d]->Write();
    _hmap_tau_transmass[d]->Write();
    _hmap_muon_transmass[d]->Write();
    _hmap_elec_transmass[d]->Write();
    _hmap_ht[d]->Write();
    _hmap_dijet_mass[d]->Write();
    _hmap_dijet_deltaEta[d]->Write();
  }
  theFile->Close();
}

PhenoAnalysis::~PhenoAnalysis()
{
  // do anything here that needs to be done at desctruction time
}

double PhenoAnalysis::calculateE(double eta, double pt, double mass)
{

  double theta = TMath::ATan(TMath::Exp(-eta));
  double sin_theta = TMath::Sin(2 * theta);
  double p = pt / sin_theta;
  double e = sqrt(pow(p, 2) + pow(mass, 2));

  return e;
}

double PhenoAnalysis::calculate_deltaR(TLorentzVector vector, Track *track)
{

  double eta1 = vector.Eta();
  double phi1 = vector.Phi();
  double eta2 = track->Eta;
  double phi2 = track->Phi;
  double deltaR = sqrt(pow(eta1 - eta2, 2) + pow(phi1 - phi2, 2));
  return deltaR;
}

double PhenoAnalysis::normalizedDphi(double phi)
{
  const double PI = 3.141592653589793238463;
  double twoPI = 2.0 * PI;
  if (phi < -PI)
  {
    phi += twoPI;
  }
  if (phi > PI)
  {
    phi = twoPI - phi;
  }
  else
    phi = TMath::Abs(phi);
  return phi;
}

void PhenoAnalysis::createHistoMaps(int directories)
{
  for (Int_t i = 0; i < directories; i++)
  {
    _hmap_Nevents[i] = new TH1F("Nevents", "Nevents", 3, 0, 3);
    _hmap_lead_jet_pT[i] = new TH1F("jetLeadPt", "p_{T} Leading Jet", 200, 0., 2000.);
    _hmap_lead_jet_eta[i] = new TH1F("jetLeadEta", "#eta Leading Jet", 50, -5.0, 5.0);
    _hmap_lead_jet_phi[i] = new TH1F("jetLeadPhi", "#phi Leading Jet", 70, -3.6, 3.6);
    _hmap_slead_jet_pT[i] = new TH1F("jetSleadPt", "p_{T} Sub-leading Jet", 200, 0., 2000.);
    _hmap_slead_jet_eta[i] = new TH1F("jetSleadEta", "#eta Sub-leading Jet", 50, -5.0, 5.0);
    _hmap_slead_jet_phi[i] = new TH1F("jetSleadPhi", "#phi Sub-leading Jet", 70, -3.6, 3.6);
    _hmap_n_jets[i] = new TH1F("nJets", "N(jet)", 7, 0, 7);
    _hmap_n_tau[i] = new TH1F("nTaus", "N(#tau)", 4, 0, 4);
    _hmap_tau1_pT[i] = new TH1F("tau1Pt", "p_{T}(#tau_{1})", 200, 0., 2000.);
    _hmap_tau1_eta[i] = new TH1F("tau1Eta", "#eta(#tau_{1})", 50, -3.5, 3.5);
    _hmap_tau1_phi[i] = new TH1F("tau1Phi", "#phi(#tau_{1})", 70, -3.6, 3.6);
    _hmap_tau2_pT[i] = new TH1F("tau2Pt", "p_{T}(#tau_{2})", 200, 0., 2000.);
    _hmap_tau2_eta[i] = new TH1F("tau2Eta", "#eta(#tau_{2})", 50, -3.5, 3.5);
    _hmap_tau2_phi[i] = new TH1F("tau2Phi", "#phi(#tau_{2})", 70, -3.6, 3.6);
    _hmap_muon1_pT[i] = new TH1F("muon1Pt", "p_{T}(#muon_{1})", 200, 0., 2000.);
    _hmap_muon1_eta[i] = new TH1F("muon1Eta", "#eta(#muon_{1})", 50, -3.5, 3.5);
    _hmap_muon1_phi[i] = new TH1F("muon1Phi", "#phi(#muon_{1})", 70, -3.6, 3.6);
    _hmap_muon2_pT[i] = new TH1F("muon2Pt", "p_{T}(#muon_{2})", 200, 0., 2000.);
    _hmap_muon2_eta[i] = new TH1F("muon2Eta", "#eta(#muon_{2})", 50, -3.5, 3.5);
    _hmap_muon2_phi[i] = new TH1F("muon2Phi", "#phi(#muon_{2})", 70, -3.6, 3.6);
    _hmap_elec1_pT[i] = new TH1F("elec1Pt", "p_{T}(#elec_{1})", 200, 0., 2000.);
    _hmap_elec1_eta[i] = new TH1F("elec1Eta", "#eta(#elec_{1})", 50, -3.5, 3.5);
    _hmap_elec1_phi[i] = new TH1F("elec1Phi", "#phi(#elec_{1})", 70, -3.6, 3.6);
    _hmap_elec2_pT[i] = new TH1F("elec2Pt", "p_{T}(#elec_{2})", 200, 0., 2000.);
    _hmap_elec2_eta[i] = new TH1F("elec2Eta", "#eta(#elec_{2})", 50, -3.5, 3.5);
    _hmap_elec2_phi[i] = new TH1F("elec2Phi", "#phi(#elec_{2})", 70, -3.6, 3.6);
    _hmap_tauMass[i] = new TH1F("tauMass", "m(#t)", 100, 0, 1000);
    _hmap_MET[i] = new TH1F("MET", "MET", 100, 0, 1000);
    _hmap_tau_transmass[i] = new TH1F("tauTransMass", "Tau Transverse Mass", 100, 0, 1000);
    _hmap_muon_transmass[i] = new TH1F("muonTransMass", "Muon Transverse Mass", 100, 0, 1000);
    _hmap_elec_transmass[i] = new TH1F("elecTransMass", "Electron Transverse Mass", 100, 0, 1000);
    _hmap_ht[i] = new TH1F("HT", "H_{T}", 100, 0, 5000);
    _hmap_st[i] = new TH1F("ST", "S_{T}", 100, 0, 5000);
    _hmap_dijet_mass[i] = new TH1F("diJetMass", "diJet_Mass", 100, 0, 5000);
    _hmap_dijet_deltaEta[i] = new TH1F("diJetDeltaEta", "diJet_deltaEta", 160, 0, 8);
    // Gen Taus
    _hmap_Gentau1_pT[i] = new TH1F("Gentau1Pt", "p_{T}(#tau_{1})", 200, 0., 2000.);
    _hmap_GenDelta_pT[i] = new TH1F("GenDeltaPt", "#Delta p_{T}", 1000, 0., 1000.);
    _hmap_Delta_pT[i] = new TH1F("DeltaPt", "#Delta p_{T}", 1000, 0., 1000.);
    _hmap_Gentau1_eta[i] = new TH1F("Gentau1Eta", "#eta(#tau_{1})", 50, -3.5, 3.5);
    _hmap_Gentau1_phi[i] = new TH1F("Gentau1Phi", "#phi(#tau_{1})", 70, -3.6, 3.6);
    _hmap_Gentau2_pT[i] = new TH1F("Gentau2Pt", "p_{T}(#tau_{2})", 200, 0., 2000.);
    _hmap_Gentau2_eta[i] = new TH1F("Gentau2Eta", "#eta(#tau_{2})", 50, -3.5, 3.5);
    _hmap_Gentau2_phi[i] = new TH1F("Gentau2Phi", "#phi(#tau_{2})", 70, -3.6, 3.6);
  }
}




bool compareTLVPTDescending(TLorentzVector tlv1, TLorentzVector tlv2){
  return (tlv1.Pt() > tlv2.Pt()); 
}

