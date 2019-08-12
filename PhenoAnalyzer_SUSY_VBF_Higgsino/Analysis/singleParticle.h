/*
                               __
                           .--()°'.'
 Author: Nathalia Cardona '|, . ,'
                           !_-(_\
 Cuts for single particles
*/

#include "../ROOTFunctions.h"
#include "../DelphesFunctions.h"

bool cut_e(ExRootTreeReader *treeReader,
           map<string, TClonesArray *> branchDict,
           int entry)
{
  // only get events that have 1 electron between 8 and 50 PT (GeV)

  treeReader->ReadEntry(entry);

  int elecCount = 0;

  int leaf = 0;

  while (elecCount < 2 && leaf < branchDict["Electron"]->GetEntries())
  {
    Electron *elec = (Electron *)branchDict["Electron"]->At(leaf);
    if (elec->PT >= 8 && elec->PT <= 50)
    {
      elecCount += 1;
    }

    leaf++;
  }

  return elecCount == 1;
}

bool cut_mu(ExRootTreeReader *treeReader,
            map<string, TClonesArray *> branchDict,
            int entry)
{
  // only get events that have 1 muon between 8 and 50 PT (GeV)
  // and dont have a single electron between 8 and 50 PT (GeV)

  if (!cut_e(treeReader, branchDict, entry))
  {

    treeReader->ReadEntry(entry);

    int muonCount = 0;

    int leaf = 0;

    while (muonCount < 2 && leaf < branchDict["Muon"]->GetEntries())
    {
      Muon *muon = (Muon *)branchDict["Muon"]->At(leaf);
      if (muon->PT >= 8 && muon->PT <= 50)
      {
        muonCount += 1;
      }

      leaf++;
    }

    return muonCount == 1;
  }
  else
  {
    return false;
  }
}

bool cut_tau(ExRootTreeReader *treeReader,
             map<string, TClonesArray *> branchDict,
             int entry)
{
  // only get events that have 1 tau between 20 and 50 PT (GeV)
  // and dont have a single electron between 8 and 50 PT (GeV)
  // and dont have a single muon between 8 and 50 PT (GeV)

  if (!cut_e(treeReader, branchDict, entry) && !cut_mu(treeReader, branchDict, entry))
  {

    treeReader->ReadEntry(entry);

    int tauCount = 0;

    int leaf = 0;

    while (tauCount < 2 && leaf < branchDict["Jet"]->GetEntries())
    {
      Jet *jet = (Jet *)branchDict["Jet"]->At(leaf);
      if (jet->TauTag == 1)
      {
        if (jet->PT >= 20 && jet->PT <= 50)
        {
          tauCount += 1;
        }
      }

      leaf++;
    }

    return tauCount == 1;
  }
  else
  {
    return false;
  }
}
