#ifndef PHYSICS_H
#define PHYSICS_H

#include "./ROOTFunctions.h"
#include "./DelphesFunctions.h"

#include <utility>

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

// Calculate the delta R between two jets
float deltaR4TLorentzVector(TLorentzVector t1, TLorentzVector t2)
{
  double dphi = abs(t1.Phi() - t2.Phi());
  const double PI = TMath::Pi();
  if (dphi > PI)
  {
    dphi = 2 * PI - dphi;
  }
  float DEta = t1.Eta() - t2.Eta();
  float dR = sqrt(pow(dphi, 2) + pow(DEta, 2));
  return dR;
}

// Calculate the mjj
float mjj(Jet *j1, Jet *j2)
{
  return sqrt(2. * j1->PT * j2->PT * cosh(j1->Eta - j2->Eta));
}

bool min2JetsNotTau(ExRootTreeReader *treeReader, map<string, TClonesArray *> branchDict, int entry)
{
  treeReader->ReadEntry(entry);

  // cout<<"in min2 jets"<<endl; 
  // cout << "num jets "<< branchDict["Jet"]->GetEntries() << endl;

  Jet *jet;
  vector<Jet *> posJets;

  for (int i = 0; i < branchDict["Jet"]->GetEntries(); i++)
  {
    jet = (Jet *)branchDict["Jet"]->At(i);

    if (!jet->TauTag)
    {
      posJets.push_back(jet);
    }
  }

  // cout<<"after for num array "<< (0 < (posJets.size() - 1)) << " " << ((int) posJets.size() - 1) << endl;

  bool ans = false;
  Jet *j1;
  Jet *j2;

  for(int i = 0; ( i < (int)posJets.size() - 1) && !ans; i++)
  {
    // cout << "entra al loop" << endl;
    j1 = posJets[i];

    for(int j = i + 1; ( j < (int)posJets.size()) && !ans; j++)
    {
      j2 = posJets[j];

      if (deltaR(j1, j2) > 0.3 )
      {
        ans = true;
      }
    }
  }

  return ans;
}

pair<Jet *, Jet *> getMaxMjjJetPair(ExRootTreeReader *treeReader, map<string, TClonesArray *> branchDict, int entry)
{
  treeReader->ReadEntry(entry);

  // cout<<"entered function"<<endl;

  if (min2JetsNotTau(treeReader, branchDict, entry))
  {

    // cout<<"in if"<<endl;

    Jet *j1;
    Jet *j2;
    float topMjj = - 1.;
    float tempMjj = 0.;
    pair <Jet *, Jet *> ans (NULL, NULL);

    for
       (int i = 0; i < branchDict["Jet"]->GetEntries() - 1; i++)
      {
        j1 = (Jet *)branchDict["Jet"]->At(i);
        if (!j1->TauTag)
        {
        for
           (int j = i + 1; j < branchDict["Jet"]->GetEntries(); j++)
          {
            j2 = (Jet *)branchDict["Jet"]->At(j);
            if (!j2->TauTag)
            {
              if
                 (deltaR(j1, j2) > 0.3 )
                {
                  tempMjj = mjj(j1, j2);
                  if
                     (tempMjj > topMjj)
                    {
                      topMjj = tempMjj;
                      ans.first = j1;
                      ans.second = j2;
                    
                    }
                }
            }
          }
        }
      }
    return ans;
  }
  else
  {
    // cout<<"in else"<<endl;

    pair <Jet *, Jet *> ans (NULL, NULL);
    return ans;
  }
}

#endif
