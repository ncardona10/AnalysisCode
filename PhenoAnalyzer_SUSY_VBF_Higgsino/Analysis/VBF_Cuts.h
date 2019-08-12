/*
                               __
                           .--()°'.'
 Author: Nathalia Cardona '|, . ,'
                           !_-(_\
 apply VBF cuts
*/

#include "../ROOTFunctions.h"
#include "../DelphesFunctions.h"

bool singleLeptonPT(ExRootTreeReader *treeReader, map<string, TClonesArray *> branchDict, int entry)
{
  float leptonCount = 0;

  //get entry
  treeReader->ReadEntry(entry);
  float PTLowerBound = 8.0;
  float PTUpperBound = 50.0;

  // count number of total leptons

  // electrons
  for (int leaf = 0; leaf < branchDict["Electron"]->GetEntries(); leaf++)
  {
    Electron *lepton = (Electron *)branchDict["Electron"]->At(leaf);
    if (lepton->PT > PTLowerBound && lepton->PT < PTUpperBound)
    {
      leptonCount += 1.0;
    }
  }

  // muons
  for (int leaf = 0; leaf < branchDict["Muon"]->GetEntries(); leaf++)
  {
    Muon *lepton = (Muon *)branchDict["Muon"]->At(leaf);
    if (lepton->PT > PTLowerBound && lepton->PT < PTUpperBound)
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
      if (jet->PT > PTLowerBound && jet->PT < PTUpperBound)
      {
        leptonCount += 1.0;
      }
    }
  }

  return leptonCount == 1;
}

bool met(ExRootTreeReader *treeReader, map<string, TClonesArray *> branchDict, int entry)
{
  //get entry
  treeReader->ReadEntry(entry);

  MissingET *METPointer = (MissingET *)branchDict["MissingET"]->At(0);
  double MET = METPointer->MET;

  return MET > 75;
}

// Calculate delta eta
float deltaEta(Jet *leadingJet, Jet *subLeadingJet)
{
  return abs(leadingJet->Eta - subLeadingJet->Eta);
}

bool min2Jets(ExRootTreeReader *treeReader, map<string, TClonesArray *> branchDict, int entry)
{
  //get entry
  treeReader->ReadEntry(entry);

  return branchDict["Jet"]->GetEntries() >= 2;
}

bool noFilter(ExRootTreeReader *treeReader,
              map<string, TClonesArray *> branchDict,
              int entry)
{
  return true;
}

float mjj(ExRootTreeReader *treeReader,
          map<string, TClonesArray *> branchDict,
          int entry)
{
  // asume that number of jets >=2
  treeReader->ReadEntry(entry);

  Jet *leadingJet = (Jet *)branchDict["Jet"]->At(0);
  Jet *subLeadingJet = (Jet *)branchDict["Jet"]->At(1);

  float dEta = deltaEta(leadingJet, subLeadingJet);

  return TMath::Power(2 * leadingJet->PT * subLeadingJet->PT * TMath::ACosH(dEta), 0.5);
}

bool allEta(ExRootTreeReader *treeReader,
            map<string, TClonesArray *> branchDict,
            int entry)
{
  treeReader->ReadEntry(entry);
  bool ans = true;
  for (int leaf = 0; leaf < branchDict["Jet"]->GetEntries(); leaf++)
  {
    Jet *jet = (Jet *)branchDict["Jet"]->At(leaf);
    if (jet->Eta < 4.0)
    {
      ans = false;
      break;
    }
  }
  return ans;
}

bool cut1(ExRootTreeReader *treeReader,
          map<string, TClonesArray *> branchDict,
          int entry)
{
  // met>75 GeV
  // dEta <5
  // mjj > 500
  // pj > 30
  // bjets = 0

  treeReader->ReadEntry(entry);

  //bjets = 0

  // MINIMUM 2 JETS!!

  bool min2JetsBool = min2Jets(treeReader, branchDict, entry);

  if (min2JetsBool)
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

    
    bool mjjBool = mjj(treeReader, branchDict, entry) > 500;

    Jet *leadingJet = (Jet *)branchDict["Jet"]->At(0);
    Jet *subLeadingJet = (Jet *)branchDict["Jet"]->At(1);

    bool deltaEtaBool = deltaEta(leadingJet, subLeadingJet) < 5.0;

    bool metBool = met(treeReader, branchDict, entry);

    return bJetsBool && mjjBool && deltaEtaBool && metBool;
  }
  else
  {
    return false;
  }
}
	