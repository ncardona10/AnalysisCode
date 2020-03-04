/*                             __
                           .--()°'.'
 Author: Nathalia Cardona '|, . ,'
                           !_-(_\
 Removes particle overlaps                     
*/
#ifndef OVERLAPS_H
#define OVERLAPS_H


#include "../ROOTFunctions.h"
#include "../DelphesFunctions.h"
#include "./Physics.h"

bool particleOverlap(vector<TLorentzVector> tlvArray)
{
  if (tlvArray.size() > 1)
  {

    bool ans = false;

    for (int i = 0; (unsigned) i < tlvArray.size() - 1 && !ans; i++)
    {
      TLorentzVector tlv1 = tlvArray[i];

      for (int j = i + 1; (unsigned) j < tlvArray.size() && !ans; j++)
      {

        TLorentzVector tlv2 = tlvArray[j];

        double deltaR = deltaR4TLorentzVector(tlv1, tlv2);
        ans = overlap(deltaR);
      }
    }

    return ans;
  }
  else
  {
    return false;
  }
}


int elecOverlap(ExRootTreeReader *treeReader,
                map<string, TClonesArray *> branchDict,
                Jet *jet)
{
  int ans = -1;

  double jetE = calculateE(jet->PT, jet->Eta, jet->Mass);
  TLorentzVector jetTLV(jet->PT, jet->Eta, jet->Phi, jetE);

  int leaf = 0;
  while (ans == -1 && leaf < branchDict["Electron"]->GetEntries())
  {

    Electron *electron = (Electron *)branchDict["Electron"]->At(leaf);

    double electE = calculateE(electron->PT, electron->Eta, 0.000510998902);
    TLorentzVector elecTLV(electron->PT, electron->Eta, electron->Phi, electE);
    double dr = deltaR4TLorentzVector(jetTLV, elecTLV);
    if (overlap(dr))
    {
      ans = leaf;
    }

    leaf++;
  }
  return ans;
}

int muonOverlap(ExRootTreeReader *treeReader,
                map<string, TClonesArray *> branchDict,
                Jet *jet)
{
  int ans = -1;

  double jetE = calculateE(jet->PT, jet->Eta, jet->Mass);
  TLorentzVector jetTLV(jet->PT, jet->Eta, jet->Phi, jetE);

  int leaf = 0;

  while (ans == -1 && leaf < branchDict["Muon"]->GetEntries())
  {

    Muon *muon = (Muon *)branchDict["Muon"]->At(leaf);

    double muonE = calculateE(muon->PT, muon->Eta, 0.1056583715);
    TLorentzVector muonTLV(muon->PT, muon->Eta, muon->Phi, muonE);

    double dr = deltaR4TLorentzVector(jetTLV, muonTLV);

    if (overlap(dr))
    {
      ans = leaf;
    }

    leaf++;
  }
  return ans;
}

int tauOverlap(ExRootTreeReader *treeReader,
                map<string, TClonesArray *> branchDict,
                Jet *jet)
{
  int ans = -1;
  double jetE = calculateE(jet->PT, jet->Eta, jet->Mass);
  TLorentzVector jetTLV(jet->PT, jet->Eta, jet->Phi, jetE);

  int leaf = 0;
  while (ans == -1 && leaf < branchDict["Jet"]->GetEntries())
  {
    
    Jet *tau = (Jet *)branchDict["Jet"]->At(leaf);
    if(tau->TauTag) 
    {
      double tauE = calculateE(tau->PT, tau->Eta, tau->Mass);
      TLorentzVector tauTLV(tau->PT, tau->Eta, tau->Phi, tauE);

      double dr = deltaR4TLorentzVector(jetTLV, tauTLV);

      if (overlap(dr))
      {
        ans = leaf;
      }

    }


    leaf++;
  }
  return ans;
}

#endif
