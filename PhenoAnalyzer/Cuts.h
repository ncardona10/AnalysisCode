/*
                               __
                           .--()°'.'
 Author: Nathalia Cardona '|, . ,'
                           !_-(_\
 Cuts for single particles
*/

#include "./ROOTFunctions.h"
#include "./DelphesFunctions.h"
#include "./Physics.h"
#include "./Overlaps.h"
#include "./LeptonCounter.h"

bool noFilter(ExRootTreeReader *treeReader,
              map<string, TClonesArray *> branchDict,
              int entry,
              vector<int> &cutsArr,
              vector<pair<Jet*,Jet*>> &jetPairs)
{
  return true;
}

bool compareJetsByPT(Jet *j1, Jet *j2)
{
  return (j1->PT > j2->PT);
}

bool met(ExRootTreeReader *treeReader,
         map<string, TClonesArray *> branchDict,
         int entry,
         vector<int> &cutsArr,
              vector<pair<Jet*,Jet*>> &jetPairs)
{
  // met>150

  if (cutsArr[entry])
  {
    treeReader->ReadEntry(entry);

    bool metBool = met(treeReader, branchDict, entry) > 150;

    cutsArr[entry] = metBool;
    return metBool;
  }

  return false;
}

bool bjets(ExRootTreeReader *treeReader,
           map<string, TClonesArray *> branchDict,
           int entry,
           vector<int> &cutsArr,
              vector<pair<Jet*,Jet*>> &jetPairs)
{
  // # bjets = 0

  treeReader->ReadEntry(entry);

  bool bJetsBool = true;

  if (cutsArr[entry])
  {

    for (int leaf = 0; leaf < branchDict["Jet"]->GetEntries(); leaf++)
    {
      Jet *jet = (Jet *)branchDict["Jet"]->At(leaf);
      if (jet->BTag == 1 && jet->PT>30.0 && abs(jet->Eta)<2.5)
      {
        bJetsBool = false;
        break;
      }
    }
  }
  else
  {
    bJetsBool = false;
  }
  cutsArr[entry] = bJetsBool;
  return bJetsBool;
}

bool vbfCut(ExRootTreeReader *treeReader,
            map<string, TClonesArray *> branchDict,
            int entry,
            vector<int> &cutsArr,
              vector<pair<Jet*,Jet*>> &jetPairs)

{
  bool ans = 0;

  // mjj > 0 #
  // eta_j1 * eta_j2 < 0 #
  // |dEta| >5.5

  treeReader->ReadEntry(entry);

  if (cutsArr[entry])
  {
    bool min2JetsBool = min2JetsNotTau(treeReader, branchDict, entry);

    if (min2JetsBool)
    {
      bool mjjBool = mjj(jetPairs[entry].first, jetPairs[entry].second) > 0;
      bool etaMultBool = (jetPairs[entry].first->Eta) * (jetPairs[entry].second->Eta) < 0;
      bool deltaEtaBool = deltaEta(jetPairs[entry].first, jetPairs[entry].second) > 5.5;
      ans = mjjBool && etaMultBool && deltaEtaBool ;
    }
  }
  cutsArr[entry] = ans;
  return ans;
}

bool nParticle(ExRootTreeReader *treeReader,
               map<string, TClonesArray *> branchDict,
               int entry,
               int n_electrons,
               int n_muon,
               int n_tau,
               vector<int> &cutsArr,
               bool checkElecPT,
               bool checkMuonPT)
{

  /*
    must comply with VBF cuts  & cuts
    n_electrons electron with:
      pt>8 
      abs(eta) < 2.4
    n_tau taus with:
      pt>20
      abs(eta)<2.4
    n_muon muons with:
      pt>5
      abs(eta)<2.4

  */

  treeReader->ReadEntry(entry);

  // vbfcut & cuts
  if (cutsArr[entry])
  {
    // verify electron condition
    int nElectrons = 0;
    int i = 0;
    vector<TLorentzVector> particleCharacteristics;

    while (nElectrons <= n_electrons && i < branchDict["Electron"]->GetEntries())
    {
      Electron *elec = (Electron *)branchDict["Electron"]->At(i);
      if (elec->PT >= 8 && abs(elec->Eta) < 2.4)
      {
        if(checkElecPT)
        {
          if(elec->PT <= 40 ){
            nElectrons++;
            particleCharacteristics.push_back(createTLorentzVector(elec->PT, elec->Eta, 0.000510998902, elec->Phi));
          }  
        }

        else
        {
          nElectrons++;
          particleCharacteristics.push_back(createTLorentzVector(elec->PT, elec->Eta, 0.000510998902, elec->Phi));
        }

      }
      i++;
    }
    // && !particleOverlap(particleCharacteristics)
    if (nElectrons == n_electrons)
    {
      //verify number of muons
      int nMuons = 0;
      int i = 0;

      particleCharacteristics.clear();

      while (nMuons <= n_muon && i < branchDict["Muon"]->GetEntries())
      {
        Muon *muon = (Muon *)branchDict["Muon"]->At(i);
        if (muon->PT >= 5 && abs(muon->Eta) < 2.4)
        {

          if(checkMuonPT){
            if(muon->PT <= 40){
              nMuons++;
              particleCharacteristics.push_back(createTLorentzVector(muon->PT, muon->Eta, 0.1056583715, muon->Phi));
            }
          }
          else{
            nMuons++;
            particleCharacteristics.push_back(createTLorentzVector(muon->PT, muon->Eta, 0.1056583715, muon->Phi));
          }
        }
        i++;
      }
      //  && !particleOverlap(particleCharacteristics)
      if (nMuons == n_muon)
      {
        //verify number of taus
        int nTaus = 0;
        int i = 0;
        particleCharacteristics.clear();

        while (nTaus <= n_tau && i < branchDict["Jet"]->GetEntries())
        {
          Jet *jet = (Jet *)branchDict["Jet"]->At(i);
          if (jet->TauTag == 1)
          {
            if (jet->PT >= 20.0 && abs(jet->Eta) < 2.4)
            {
              nTaus++;
              particleCharacteristics.push_back(createTLorentzVector(jet->PT, jet->Eta, jet->Mass, jet->Phi));
            }
          }
          i++;
        }
        return nTaus == n_tau && !particleOverlap(particleCharacteristics);
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }
  else
  {
    return false;
  }
}

/*
SINGLE PARTICLE
*/

bool mono_e(ExRootTreeReader *treeReader,
            map<string, TClonesArray *> branchDict,
            int entry,
            vector<int> &cutsArr,
              vector<pair<Jet*,Jet*>> &jetPairs)
{
  return nParticle(treeReader, branchDict, entry, 1, 0, 0, cutsArr,true,false);
}

bool mono_mu(ExRootTreeReader *treeReader,
             map<string, TClonesArray *> branchDict,
             int entry,
             vector<int> &cutsArr,
              vector<pair<Jet*,Jet*>> &jetPairs)
{
  return nParticle(treeReader, branchDict, entry, 0, 1, 0, cutsArr,false,true);
}

// bool mono_tau(ExRootTreeReader *treeReader,
//               map<string, TClonesArray *> branchDict,
//               int entry,
//               vector<int> &cutsArr,
//               vector<pair<Jet*,Jet*>> &jetPairs)
// {
//   return nParticle(treeReader, branchDict, entry, 0, 0, 1, cutsArr);
// }

/*
DI PARTICLE
*/

// bool di_e(ExRootTreeReader *treeReader,
//           map<string, TClonesArray *> branchDict,
//           int entry,
//           vector<int> &cutsArr,
//               vector<pair<Jet*,Jet*>> &jetPairs)
// {
//   return nParticle(treeReader, branchDict, entry, 2, 0, 0, cutsArr);
// }

// bool di_mu(ExRootTreeReader *treeReader,
//            map<string, TClonesArray *> branchDict,
//            int entry,
//            vector<int> &cutsArr,
//               vector<pair<Jet*,Jet*>> &jetPairs)
// {
//   return nParticle(treeReader, branchDict, entry, 0, 2, 0, cutsArr);
// }

// bool di_tau(ExRootTreeReader *treeReader,
//             map<string, TClonesArray *> branchDict,
//             int entry,
//             vector<int> &cutsArr,
//               vector<pair<Jet*,Jet*>> &jetPairs)
// {
//   return nParticle(treeReader, branchDict, entry, 0, 0, 2, cutsArr);
// }

/*
0 leptons
*/

// bool zero_leptons(ExRootTreeReader *treeReader,
//                   map<string, TClonesArray *> branchDict,
//                   int entry,
//                   vector<int> &cutsArr,
//               vector<pair<Jet*,Jet*>> &jetPairs)
// {
//   return nParticle(treeReader, branchDict, entry, 0, 0, 0, cutsArr);
// }

//filter for the jets pt and |eta|
bool jetsPtEtaFilter(ExRootTreeReader *treeReader,
                     map<string, TClonesArray *> branchDict,
                     int entry,
                     vector<int> &cutsArr,
              vector<pair<Jet*,Jet*>> &jetPairs)
{
  treeReader->ReadEntry(entry);
  if (cutsArr[entry])
  {
    bool ans = false;

    if (jetPairs[entry].first->PT > 30.0 && abs(jetPairs[entry].first->Eta) < 5.0 && 
        jetPairs[entry].second->PT > 30.0 && abs(jetPairs[entry].second->Eta) < 5.0)
    {
      ans = true;
    }
    
    cutsArr[entry] = ans;

    return ans;
  }
  return false;
}


bool jetsLeptonOverlap(ExRootTreeReader *treeReader,
                     map<string, TClonesArray *> branchDict,
                     int entry,
                     vector<int> &cutsArr,
              vector<pair<Jet*,Jet*>> &jetPairs)
{
  treeReader->ReadEntry(entry);
  bool ans = elecOverlap(treeReader,branchDict,jetPairs[entry].first)==-1;
  ans &= muonOverlap(treeReader,branchDict,jetPairs[entry].first)==-1;
  ans &= elecOverlap(treeReader,branchDict,jetPairs[entry].second)==-1;
  ans &= muonOverlap(treeReader,branchDict,jetPairs[entry].second)==-1;

  return ans;

}

bool jetConditions(ExRootTreeReader *treeReader,
                     map<string, TClonesArray *> branchDict,
                     int entry,
                     vector<int> &cutsArr,
              vector<pair<Jet*,Jet*>> &jetPairs)
{
  treeReader->ReadEntry(entry);

  if(cutsArr[entry] && jetPairs[entry].first != NULL && jetPairs[entry].second != NULL){
    bool ans = jetsPtEtaFilter(treeReader,branchDict,entry,cutsArr,jetPairs);
    //jetsLeptonOverlap(treeReader,branchDict,entry,cutsArr,jetPairs) );
    cutsArr[entry] = ans; 
    return ans;
  }
  cutsArr[entry] = false;
  return false;
}



// apply filter, do not generate histograms
void applyFilter(ExRootTreeReader *treeReader,
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
  Long64_t numberOfEntries = treeReader->GetEntries();

  for (Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    if (cutsArr[entry])
    {
      filter(treeReader, branchDict, entry, cutsArr, jetPairs);
    }
  }
}
