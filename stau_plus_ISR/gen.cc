#include <iostream>
//#include <PhenoAnalyzer.h>

void gen(){
 
TFile* file = TFile::Open("Test.root");
file->cd("No_cuts");
TH1F *jetpt2 = (TH1F*)file->Get("jet_pt");
//TFile *f = new TFile("Test.root", "RECREATE");
//f->GetObject("jet_pt_max", readThis);
//file->Close();
jetpt2->Draw();
//return 0;	

}
