/*
                               __
                           .--()°'.'
 Author: Nathalia Cardona '|, . ,'
                           !_-(_\
 apply VBF cuts
*/

#include "../ROOTFunctions.h"
#include "../DelphesFunctions.h"

/*
  Returns PT of the lepton in the event if there is exactly 1 lepton. 
  Else returns -1. 
*/

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

bool met(ExRootTreeReader *treeReader, map<string, TClonesArray *> branchDict, int entry){
  //get entry
  treeReader->ReadEntry(entry);

  METpointer mp= (MissingET *)branchDict["MissingET"]->At(0);
  double MET = mp->MET;

  return MET>75; //in GeV?
}

/*
  Calculate delta eta
*/

float deltaEta(Jet leadingJet, Jet subLeadingJet)
{
  return abs(leadingJet->Eta - subLeadingJet->Eta);
}

bool min2Jets(ExRootTreeReader *treeReader, map<string, TClonesArray *> branchDict, int entry){
  //get entry
  treeReader->ReadEntry(entry);

  return branchDict["Jet"]-GetEntries() >= 2;
}

