/*
                               __
                           .--()°'.'
 Author: Nathalia Cardona '|, . ,'
                           !_-(_\
 Cuts for single particles
*/

#include "../ROOTFunctions.h"
#include "../DelphesFunctions.h"
#include "./Physics.h"
#include "./Overlaps.h"

bool noFilter(ExRootTreeReader *treeReader,
              map<string, TClonesArray *> branchDict,
              int entry,
              vector<int> &vbfCutsArr,
              vector<int> &cutsArr)
{
  bool mjjBool = mjj(treeReader, branchDict, entry) > 0;
  
  cutsArr[entry] = mjjBool;
  return mjjBool;
}

bool compareJetsByPT(Jet *j1, Jet *j2)
{
  return (j1->PT > j2->PT);
}

bool vbfCut(ExRootTreeReader *treeReader,
            map<string, TClonesArray *> branchDict,
            int entry,
            vector<int> &vbfCutsArr,
            vector<int> &cutsArr)

{
  bool ans = 0;

  // mjj > 0 #
  // eta_lead * eta_subLead < 0 #
  // |dEta| > 3.8
  // pt both jets > 30
  // |eta|<5 for both jets

  treeReader->ReadEntry(entry);

  if (cutsArr[entry])
  {

    // MINIMUM 2 JETS!!

    bool min2JetsBool = min2Jets(treeReader, branchDict, entry);

    if (min2JetsBool)
    {
      // Sort by Pt from max to min
      vector<Jet *> daJets;

      for (int i = 0; i < branchDict["Jet"]->GetEntries(); i++)
      {
        daJets.push_back((Jet *)branchDict["Jet"]->At(i));
      }

      sort(daJets.begin(), daJets.end(), compareJetsByPT);

      Jet *leadingJet = daJets[0];
      Jet *subLeadingJet = daJets[1];

      // if there are only overlaps, mjj returns -1
      bool mjjBool = mjj(treeReader, branchDict, entry) > 0;
      bool deltaMultipl = (leadingJet->Eta) * (subLeadingJet->Eta) < 0;
      bool deltaEtaBool = abs(deltaEta(leadingJet, subLeadingJet)) > 3.8;
      bool pTBothBool = leadingJet->PT > 30.0 && subLeadingJet->PT > 30.0;
      bool etaBelow5 = abs(leadingJet->Eta) < 5.0 && abs(subLeadingJet->Eta) < 5.0;

      ans = mjjBool && deltaMultipl && deltaEtaBool && pTBothBool && etaBelow5;
    }
  }
  cutsArr[entry] = ans;
  vbfCutsArr[entry] = ans;
  return ans;
}

bool cuts(ExRootTreeReader *treeReader,
          map<string, TClonesArray *> branchDict,
          int entry,
          vector<int> &vbfCutsArr,
          vector<int> &cutsArr)
{
  // needs to meet vbf
  // met>250
  // # bjets = 0

  treeReader->ReadEntry(entry);
  bool ans = false;

  // if (vbfCutsArr[entry])
  // {
  bool metBool = met(treeReader, branchDict, entry) > 250;

  if (metBool)
  {
    bool bJetsBool = true;

    for (int leaf = 0; leaf < branchDict["Jet"]->GetEntries(); leaf++)
    {
      Jet *jet = (Jet *)branchDict["Jet"]->At(leaf);
      if (jet->BTag == 1)
      {
        bJetsBool = false;
        break;
      }
    }

    ans = bJetsBool;
  }
  // }

  cutsArr[entry] = ans;
  return ans;
}

bool met(ExRootTreeReader *treeReader,
         map<string, TClonesArray *> branchDict,
         int entry,
         vector<int> &vbfCutsArr,
         vector<int> &cutsArr)
{
  // met>250

  if (cutsArr[entry])
  {
    treeReader->ReadEntry(entry);

    bool metBool = met(treeReader, branchDict, entry) > 250;

    cutsArr[entry] = metBool;
    return metBool;
  }

  return false;
}

bool bjets(ExRootTreeReader *treeReader,
           map<string, TClonesArray *> branchDict,
           int entry,
           vector<int> &vbfCutsArr,
           vector<int> &cutsArr)
{
  // # bjets = 0

  treeReader->ReadEntry(entry);

  bool bJetsBool = true;

  if (cutsArr[entry])
  {

    for (int leaf = 0; leaf < branchDict["Jet"]->GetEntries(); leaf++)
    {
      Jet *jet = (Jet *)branchDict["Jet"]->At(leaf);
      if (jet->BTag == 1)
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
bool nParticle(ExRootTreeReader *treeReader,
               map<string, TClonesArray *> branchDict,
               int entry,
               int n_electrons,
               int n_muon,
               int n_tau,
               vector<int> &vbfCutsArr,
               vector<int> &cutsArr)
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

        nElectrons++;
        particleCharacteristics.push_back(createTLorentzVector(elec->PT, elec->Eta, 0.000510998902, elec->Phi));
      }
      i++;
    }

    if (nElectrons == n_electrons && !particleOverlap(particleCharacteristics))
    {
      //verify number of muons
      int nMuons = 0;
      int i = 0;

      particleCharacteristics.clear();

      while (nMuons <= n_muon && i < branchDict["Muon"]->GetEntries())
      {
        Muon *muon = (Muon *)branchDict["Muon"]->At(i);
        if (muon->PT >= 5  && abs(muon->Eta) < 2.4)
        {
          nMuons++;
          particleCharacteristics.push_back(createTLorentzVector(muon->PT, muon->Eta, 0.1056583715, muon->Phi));
        }
        i++;
      }

      if (nMuons == n_muon && !particleOverlap(particleCharacteristics))
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
            vector<int> &vbfCutsArr,
            vector<int> &cutsArr)
{
  return nParticle(treeReader, branchDict, entry, 1, 0, 0, vbfCutsArr, cutsArr);
}

bool mono_mu(ExRootTreeReader *treeReader,
             map<string, TClonesArray *> branchDict,
             int entry,
             vector<int> &vbfCutsArr,
             vector<int> &cutsArr)
{
  return nParticle(treeReader, branchDict, entry, 0, 1, 0, vbfCutsArr, cutsArr);
}

bool mono_tau(ExRootTreeReader *treeReader,
              map<string, TClonesArray *> branchDict,
              int entry,
              vector<int> &vbfCutsArr,
              vector<int> &cutsArr)
{
  return nParticle(treeReader, branchDict, entry, 0, 0, 1, vbfCutsArr, cutsArr);
}

/*
DI PARTICLE
*/

bool di_e(ExRootTreeReader *treeReader,
          map<string, TClonesArray *> branchDict,
          int entry,
          vector<int> &vbfCutsArr,
          vector<int> &cutsArr)
{
  return nParticle(treeReader, branchDict, entry, 2, 0, 0, vbfCutsArr, cutsArr);
}

bool di_mu(ExRootTreeReader *treeReader,
           map<string, TClonesArray *> branchDict,
           int entry,
           vector<int> &vbfCutsArr,
           vector<int> &cutsArr)
{
  return nParticle(treeReader, branchDict, entry, 0, 2, 0, vbfCutsArr, cutsArr);
}

bool di_tau(ExRootTreeReader *treeReader,
            map<string, TClonesArray *> branchDict,
            int entry,
            vector<int> &vbfCutsArr,
            vector<int> &cutsArr)
{
  return nParticle(treeReader, branchDict, entry, 0, 0, 2, vbfCutsArr, cutsArr);
}

/*
0 leptons
*/

bool zero_leptons(ExRootTreeReader *treeReader,
                  map<string, TClonesArray *> branchDict,
                  int entry,
                  vector<int> &vbfCutsArr,
                  vector<int> &cutsArr)
{
  return nParticle(treeReader, branchDict, entry, 0, 0, 0, vbfCutsArr, cutsArr);
}
