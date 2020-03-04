#ifndef PHYSICS_H
#define PHYSICS_H

#include "../ROOTFunctions.h"
#include "../DelphesFunctions.h"

//Calculate the energy of the particle
double calculateE(double eta, double pt, double mass)
{

  double theta = TMath::ATan(TMath::Exp(-eta));
  double sin_theta = TMath::Sin(2 * theta);
  double p = pt / sin_theta;
  double e = sqrt(pow(p, 2) + pow(mass, 2));

  return e;
}

//Create a TLorenzt vector of a particle
TLorentzVector createTLorentzVector(double PT, double Eta, double Mass, double Phi)
{
  double E = calculateE(PT, Eta, Mass);
  TLorentzVector TLV(PT, Eta, Phi, E);
  return TLV;
}

//Calculate the delta R between two TLorenzt vectors
double dR(TLorentzVector t1, TLorentzVector t2)
{
  //no need to get abs
  //returns sqrt of the sum of 2 squares
  return t1.DeltaR(t2);
}

//Boolean true or false if there is an overlap
bool overlap(double dr)
{
  //well take minimum dr as 0.3
  //TO-DO: make this a constant to take from config file
  return dr < 0.3;
}

// Normalize variable phi
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
  return phi;
}

// Return the missing energy transverse of an event
double met(ExRootTreeReader *treeReader, map<string, TClonesArray *> branchDict, int entry)
{
  //get entry
  treeReader->ReadEntry(entry);

  MissingET *METPointer = (MissingET *)branchDict["MissingET"]->At(0);
  double MET = METPointer->MET;

  return MET;
}

// Calculate the transverse mass based on a lepton and the missing energy transverse
double mt(double pt, double met, double deltaPhi)
{
  return TMath::Power(
      2 * pt * abs(met) * (1 - TMath::Cos(deltaPhi)),
      0.5);
}

// Calculate delta eta between two jets
float deltaEta(Jet *j1, Jet *j2)
{
  return abs(j1->Eta - j2->Eta);
}

//Calculate the delta eta between the leading and subleading jets
float dEtaLeadSubLead(ExRootTreeReader *treeReader,
                      map<string, TClonesArray *> branchDict,
                      int entry)
{
  Jet *leadingJet = (Jet *)branchDict["Jet"]->At(0);
  Jet *subLeadingJet = (Jet *)branchDict["Jet"]->At(1);

  return deltaEta(leadingJet, subLeadingJet);
}

// Calculate the delta R between two jets
float deltaR(Jet *j1, Jet *j2)
{
  float dphi = abs(j1->Phi - j2->Phi);
  const double PI = TMath::Pi();
  if (dphi > PI)
  {
    dphi = 2 * PI - dphi;
  }
  float DEta = j1->Eta - j2->Eta;
  float dR = sqrt(pow(dphi, 2) + pow(DEta, 2));
  return dR;
}

// Calculate the mjj
float mjjValue(Jet *j1, Jet *j2)
{
  return sqrt(2. * j1->PT * j2->PT * cosh(j1->Eta - j2->Eta));
}

//Calculate the delta eta between the two jets that max mjj
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

      if (deltaR(j1, j2) > 0.3)
      {
        tempMjj = mjjValue(j1, j2);

        if (tempMjj > topMjj)
        {
          topMjj = tempMjj;
          returnDEta = abs(deltaEta(j1, j2));
        }
      }
    }
  }

  return returnDEta;
}

//
float mjj(ExRootTreeReader *treeReader,
          map<string, TClonesArray *> branchDict,
          int entry)
{
  // asume that number of jets >=2
  treeReader->ReadEntry(entry);

  float topMjj = -1.;
  float tempMjj = 0.;
  Jet *j1;
  Jet *j2;

  for (int i = 0; i < branchDict["Jet"]->GetEntries() - 1; i++)
  {
    j1 = (Jet *)branchDict["Jet"]->At(i);

    for (int j = i + 1; j < branchDict["Jet"]->GetEntries(); j++)
    {
      j2 = (Jet *)branchDict["Jet"]->At(j);

      if (deltaR(j1, j2) > 0.3)
      {
        tempMjj = mjjValue(j1, j2);
        //cout << "pT_j1 "<<j1->PT << "pT_j2 "<<j2->PT<< " DeltaEta "<<j1->Eta - j2->Eta << "  mass "<< mass_jj<<endl;
        if (tempMjj > topMjj)
        {
          topMjj = tempMjj;
        }
      }
    }
  }

  return topMjj;
}

//
bool min2Jets(ExRootTreeReader *treeReader, map<string, TClonesArray *> branchDict, int entry)
{
  treeReader->ReadEntry(entry);
  return branchDict["Jet"]->GetEntries() >= 2;
}

#endif
