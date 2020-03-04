#ifndef PHYSICS_H
#define PHYSICS_H

#include "../ROOTFunctions.h"
#include "../DelphesFunctions.h"

double calculateE(double eta, double pt, double mass)
{

  double theta = TMath::ATan(TMath::Exp(-eta));
  double sin_theta = TMath::Sin(2 * theta);
  double p = pt / sin_theta;
  double e = sqrt(pow(p, 2) + pow(mass, 2));

  return e;
}

TLorentzVector createTLorentzVector(double PT, double Eta, double Mass, double Phi)
{
  double E = calculateE(PT, Eta, Mass);
  TLorentzVector TLV(PT, Eta, Phi, E);
  return TLV;
}

double dR(TLorentzVector t1, TLorentzVector t2)
{

  //no need to get abs
  //returns sqrt of the sum of 2 squares
  return t1.DeltaR(t2);
}

bool overlap(double dr)
{
  //well take minimum dr as 0.3
  //TO-DO: make this a constant to take from config file
  return dr < 0.3;
}

double normalizedDphi(double phi)
{
  const double PI = TMath::Pi();
  double twoPI = 2.0 * PI;
  if (phi < -PI)
  {
    phi += twoPI;
  }
  if (phi > PI)
  {
    phi = twoPI - phi;
  }
  // else{
  //   phi = TMath::Abs(phi);
  // }
  return phi;
}

double met(ExRootTreeReader *treeReader, map<string, TClonesArray *> branchDict, int entry)
{
  //get entry
  treeReader->ReadEntry(entry);

  MissingET *METPointer = (MissingET *)branchDict["MissingET"]->At(0);
  double MET = METPointer->MET;

  return MET;
}

double mt(double pt, double met, double deltaPhi)
{
  return TMath::Power(
      2 * pt * abs(met) * (1 - TMath::Cos(deltaPhi)),
      0.5);
}

// Calculate delta eta
float deltaEta(Jet *j1, Jet *j2)
{
  return abs(j1->Eta - j2->Eta);
}

float dEtaLeadSubLead(ExRootTreeReader *treeReader,
                      map<string, TClonesArray *> branchDict,
                      int entry)
{
  Jet *leadingJet = (Jet *)branchDict["Jet"]->At(0);
  Jet *subLeadingJet = (Jet *)branchDict["Jet"]->At(1);

  return deltaEta(leadingJet, subLeadingJet);
}

float dEtaMaxMjj(ExRootTreeReader *treeReader,
                 map<string, TClonesArray *> branchDict,
                 int entry)
{
  // asume that number of jets >=2
  treeReader->ReadEntry(entry);

  float topMjj = -1;
  float tempMjj;
  float returnDEta = -1;
  Jet *j1;
  Jet *j2;

  for (int i = 0; i < branchDict["Jet"]->GetEntries() - 1; i++)
  {
    j1 = (Jet *)branchDict["Jet"]->At(i);

    for (int j = i + 1; j < branchDict["Jet"]->GetEntries(); j++)
    {
      j2 = (Jet *)branchDict["Jet"]->At(j);

      if (j1->PT >= 30 && j2->PT >= 30 && abs(j1->Eta) < 5.0 && abs(j2->Eta) < 5.0 && j1.DeltaR(j2) > 0.3)
      {
        tempMjj = mjjValue(j1, j2);

        if (tempMjj > topMjj)
        {
          topMjj = tempMjj;
          returnDEta = deltaEta(j1, j2);
        }
      }
    }
  }

  return returnDEta;
}

float mjjValue(Jet *j1, Jet *j2)
{

  float dEta = deltaEta(j1, j2);
  return TMath::Power(2 * j1->PT * j2->PT * TMath::ACosH(dEta), 0.5);
}

float mjj(ExRootTreeReader *treeReader,
          map<string, TClonesArray *> branchDict,
          int entry)
{
  // asume that number of jets >=2
  treeReader->ReadEntry(entry);

  float topMjj = -1;
  float tempMjj;
  Jet *j1;
  Jet *j2;

  for (int i = 0; i < branchDict["Jet"]->GetEntries() - 1; i++)
  {
    j1 = (Jet *)branchDict["Jet"]->At(i);

    for (int j = i + 1; j < branchDict["Jet"]->GetEntries(); j++)
    {
      j2 = (Jet *)branchDict["Jet"]->At(j);

      if (j1->PT >= 30 && j2->PT >= 30 && abs(j1->Eta) < 5.0 && abs(j2->Eta) < 5.0 && j1.DeltaR(j2) > 0.3)
      {
        tempMjj = mjjValue(j1, j2);

        if (tempMjj > topMjj)
        {
          topMjj = tempMjj;
        }
      }
    }
  }

  return topMjj;
}

bool min2Jets(ExRootTreeReader *treeReader, map<string, TClonesArray *> branchDict, int entry)
{
  //get entry
  treeReader->ReadEntry(entry);

  // Jet *leadingJet = (Jet *)branchDict["Jet"]->At(0);
  // Jet *subLeadingJet = (Jet *)branchDict["Jet"]->At(1);

  return branchDict["Jet"]->GetEntries() >= 2;
  // bool ans = false;
  // if (branchDict["Jet"]->GetEntries() >= 2)
  // {
  //   if (leadingJet->PT > 30 && abs(leadingJet->Eta) < 5 && subLeadingJet->PT > 30 && abs(subLeadingJet->Eta) < 5)
  //   {
  //     ans = true;
  //   }
  // }
  // return ans;
}

#endif