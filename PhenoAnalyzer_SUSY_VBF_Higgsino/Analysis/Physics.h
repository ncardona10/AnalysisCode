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


double dR(TLorentzVector t1, TLorentzVector t2){
  return t1.DeltaR(t2);
}

bool overlap(double dr){
  //well take minimum dr as 0.3
  //TO-DO: make this a constant to take from config file
  return dr<0.3;
}

double normalizedDphi(double phi)
{
  const double PI = 3.141592653589793238463;
  double twoPI = 2.0 * PI;
  if (phi < -PI)
  {
    phi += twoPI;
  }
  if (phi > PI)
  {
    phi = twoPI - phi;
  }
  else
    phi = TMath::Abs(phi);
  return phi;
}