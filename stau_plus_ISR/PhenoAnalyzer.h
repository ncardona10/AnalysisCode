////////////////////////////////////////////////////////////////
//                                                            //
// Author: Andrés Flórez, Universidad de los Andes, Colombia  //
//                                                            //
////////////////////////////////////////////////////////////////

#ifndef PHENOANALYZER_H
#define PHENOANALYZER_H

#include "ROOTFunctions.h"
#include "DelphesFunctions.h"


using namespace std;

class PhenoAnalysis {
public :
   PhenoAnalysis(TChain&, TFile*, TDirectory* dir[], int nDir);
   ~PhenoAnalysis();
   void createHistoMaps (int);
   bool overlapingObjects(double, double, double, double, double);
   double calculateE(double, double, double);
   double normalizedDphi(double);
   double calculate_deltaR(TLorentzVector,  Track*);

   // For Jets
   std::map<unsigned int, TH1*> _hmap_Nevents;
};

#endif
