/*                             __
                           .--()°'.'
 Author: Nathalia Cardona '|, . ,'
                           !_-(_\
*/

#include "PhenoAnalyzer.h"
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <iomanip>
#include <bits/stdc++.h>
#include "Analysis/LeptonCounter.h"

using namespace std;

int main(int argc, char *argv[])
{
  TFile* rootFile = new TFile("/disco1/SIMULACIONES/z+jets/z+jets_1/Events/run_01/m_delphes_events.root"); 
  TTree *MyTree; 
  rootFile->GetObject("MyTree",MyTree);
  // int nentries = MyTree->GetEntries();

  cout<< "number of entries: "<< MyTree->GetEntries() <<endl;


// MyTree->SetBranchAddress("event1",&event1);
// MyTree->SetBranchAddress("event2",&event2);



}