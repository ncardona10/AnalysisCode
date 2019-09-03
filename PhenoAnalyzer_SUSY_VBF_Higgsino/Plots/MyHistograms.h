#include "../ROOTFunctions.h"
#include "../DelphesFunctions.h"
#include <vector>
#include <string>
#include <limits>

/*
Main code taken from: 
http://www.physics.usyd.edu.au/hienergy/index.php/Local_ROOT_Examples#How_to_superimpose_Histograms_and_Histogram_derived_objects
Edited by:
Nathalia Cardona
*/

void Draw_Normalised(TObjArray histos,
                     TPad *pad = 0,
                     bool normalised = false,
                     std::string stacktitle = "",
                     float maxXAxis = std::numeric_limits<double>::quiet_NaN())
{
  // this function draws the histoname from the TObjArray superimposed
  // and normalised if required

  if (histos.GetEntries() == 0)
  {

    return;
  }

  TObjArray RootFiles;
  std::vector<std::string> legends_str;
  std::vector<int> colours = {920,
                              632,
                              416,
                              600,
                              1,
                              616,
                              432,
                              800,
                              30,
                              28,
                              900}; //Colors: https://root.cern/doc/master/classTColor.html#a0f79316b6922be594e6d00e4bc2e11eb

  for (int i = 0; i < histos.GetEntries(); i++)
  {
    TH1F *h = (TH1F *)histos[i];
    legends_str.push_back(h->GetTitle());
  }

  // lets open and draw the canvas

  TCanvas *canvas;
  if (pad == 0)
  {
    canvas = new TCanvas("c5", "TauValidation");
    pad = (TPad *)canvas->cd();
  }
  pad->cd();
  pad->SetTicks(0, 0);

  // lets take the first histoname and see if we can match a title to its which will be HS stack title

  if (stacktitle == "")
    stacktitle = ((TH1F *)histos[0])->GetTitle();
  THStack *Hs = new THStack("hs2", stacktitle.c_str());

  TLegend *legend = new TLegend(0.65, 0.5, 0.85, 0.85); // we need different positions for the legend to not

  for (int i = 0; i < histos.GetEntries(); i++)
  {

    TH1F *h = (TH1F *)histos[i];

    if (normalised)
    {
      double val1 = h->GetSumOfWeights();
      if (fabs(val1) > 0)
        h->Scale(1.0 / val1);
    }

    h->SetLineWidth(2);
    h->SetLineColor(colours[i]);
    h->SetStats(0);
    legend->AddEntry(h, legends_str[i].c_str(), "L");
    Hs->Add(h, "sames");
  }

  int no_error = 0;

  // if the array has more than 0 histograms lets specify the stat boxes and determine whether we should
  // draw errors
  if (histos.GetEntries() > 0)
  {

    // heightboxes = (float)0.5 / (float)histos.GetEntries();
    if ((strcmp(histos.At(0)->GetName(), "hist7132") == 0) || (strcmp(histos.At(0)->GetName(), "hist7032")))
      no_error = 1;
  }

  if (no_error == 1) // do not draw errors
    Hs->Draw("HIST nostack");
  else
    Hs->Draw("HISTE nostack");

  // // limit x axis range
  if (!std::isnan(maxXAxis))
  {
    Hs->GetXaxis()->SetLimits(0, maxXAxis);
  }

  pad->Update();

  legend->Draw("");

  pad->Update();
  pad->Modified(); // so it updates the pad with the new changes
  pad->Draw("");

  return;
}

TH1 *blankHistogram(string title, string filename, int bins, float min_x, float min_y)
{

  char char_array[title.length() + 1];
  strcpy(char_array, title.c_str());

  char charArrayHistoName[filename.length() + 1];
  strcpy(charArrayHistoName, filename.c_str());

  return new TH1F(charArrayHistoName, char_array, bins, min_x, min_y);
}