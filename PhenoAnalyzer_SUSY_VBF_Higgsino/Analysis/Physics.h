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

double dR(TLorentzVector t1, TLorentzVector t2)
{
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
float deltaEta(Jet *leadingJet, Jet *subLeadingJet)
{
  return abs(leadingJet->Eta - subLeadingJet->Eta);
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


bool min2Jets(ExRootTreeReader *treeReader, map<string, TClonesArray *> branchDict, int entry)
{
  //get entry
  treeReader->ReadEntry(entry);
  
  Jet *leadingJet = (Jet *)branchDict["Jet"]->At(0);
  Jet *subLeadingJet = (Jet *)branchDict["Jet"]->At(1);

  bool ans = false;
  if(branchDict["Jet"]->GetEntries() >= 2){
    if(leadingJet->PT>30 && abs(leadingJet->eta)<5 && subLeadingJet->PT>30 && abs(subLeadingJet->eta)<5){
      ans = true; 
    }
  }
  return ans;
}



#endif 