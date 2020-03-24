#include "./ROOTFunctions.h"
#include "./DelphesFunctions.h"
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
  std::vector<string> legends_str;
  std::vector<int> colours = {
      //1,
      632,
      600,
      820,
      432,
      880,
      800,
      800,

      820,
      28,
      900}; //Colors: https://root.cern/doc/master/classTColor.html#a0f79316b6922be594e6d00e4bc2e11eb

  map<string, string> realNames{
      {"m_n2_100_c1_80_n1_60", "m(#tilde{#chi}^{0}_{2})=100GeV, m(#tilde{#chi}^{#pm}_{1}) = 80GeV, m(#tilde{#chi}^{0}_{1})=60GeV"},
      {"m_n2_100_c1_75_n1_50", "m(#tilde{#chi}^{0}_{2})=100GeV, m(#tilde{#chi}^{#pm}_{1}) = 75GeV, m(#tilde{#chi}^{0}_{1})=50GeV"},
      {"m_n2_400_c1_385_n1_370", "m(#tilde{#chi}^{0}_{2})=400GeV, m(#tilde{#chi}^{#pm}_{1}) = 385GeV, m(#tilde{#chi}^{0}_{1})=370GeV"},
      {"m_n2_200_c1_175_n1_150", "m(#tilde{#chi}^{0}_{2})=200GeV, m(#tilde{#chi}^{#pm}_{1}) = 175GeV, m(#tilde{#chi}^{0}_{1})=150GeV"},
      {"wz", "wz"},
      {"zz", "zz"},
      {"ww", "ww"},
      {
          "w+jets",
          "w+jets",
      },
      {"z+jets", "z+jets"},
      {"ttbar", "ttbar"}

  };

  for (int i = 0; i < histos.GetEntries(); i++)
  {
    TH1F *h = (TH1F *)histos[i];
    legends_str.push_back(realNames[h->GetTitle()]);
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

  // with first histo title
  // THStack *Hs = new THStack("hs2", stacktitle.c_str());
  // no title
  THStack *Hs = new THStack("hs2", "");

  TLegend *legend = new TLegend(0.65, 0.5, 0.9, 0.85); // we need different positions for the legend to not

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
    if (i < 2)
    {
      h->SetLineStyle(5);
      h->SetLineWidth(4);
    }
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

void Draw_Stacked(TObjArray histos,
                  TPad *pad = 0,
                  bool normalised = false,
                  std::string stacktitle = "",
                  float maxXAxis = std::numeric_limits<double>::quiet_NaN())
{ // this function draws histograms stacked and correctly takes into account the
  // stats boxes for each

  if (histos.GetEntries() == 0)
    return; // nothing to do

  // Initial set up
  TObjArray statsboxes;
  std::vector<std::string> legends_str;
  std::vector<int> colours = {
      // 1,
      // 632,
      // 600,
      820,
      432,
      880,
      800,
      800,

      820,
      28,
      900};

  map<string, string> realNames{
      {"m_n2_100_c1_80_n1_60", "m(#tilde{#chi}^{0}_{2})=100GeV, m(#tilde{#chi}^{#pm}_{1}) = 80GeV, m(#tilde{#chi}^{0}_{1})=60GeV"},
      {"m_n2_100_c1_75_n1_50", "m(#tilde{#chi}^{0}_{2})=100GeV, m(#tilde{#chi}^{#pm}_{1}) = 75GeV, m(#tilde{#chi}^{0}_{1})=50GeV"},
      {"m_n2_400_c1_385_n1_370", "m(#tilde{#chi}^{0}_{2})=400GeV, m(#tilde{#chi}^{#pm}_{1}) = 385GeV, m(#tilde{#chi}^{0}_{1})=370GeV"},
      {"m_n2_200_c1_175_n1_150", "m(#tilde{#chi}^{0}_{2})=200GeV, m(#tilde{#chi}^{#pm}_{1}) = 175GeV, m(#tilde{#chi}^{0}_{1})=150GeV"},
      {"wz", "wz"},
      {"zz", "zz"},
      {"ww", "ww"},
      {
          "w+jets",
          "w+jets",
      },
      {"z+jets", "z+jets"},
      {"ttbar", "ttbar"}

  };

  TPaveStats *st1;

  const int n = histos.GetEntries();

  // lets open and draw the canvas
  TCanvas *canvas;
  if (pad == 0)
  {
    canvas = new TCanvas("c5", "Stacked Histograms");
    pad = (TPad *)canvas->cd();
  }
  pad->cd();
  pad->SetTicks(0, 0);
  pad->SetRightMargin(0.20);

  // lets take the first histoname and see if we can match a title to its which will be HS stack title
  if (stacktitle == "")
    stacktitle = ((TH1F *)histos[0])->GetTitle();
  THStack *Hs = new THStack("hs2", stacktitle.c_str());

  // Set Axis Units
  //  Hs->GetXaxis()->SetTitle( ((TH1F*)histos[0])->GetXaxis()->GetTitle() );

  // Set up the LEGEND
  TLegend *legend = new TLegend(0.80, 0.3, 0.995, 0.4); // we need different positions for the legend to not
                                                        // get the plot titles for the legend
  for (int i = n; i > 0; i--)
  {

    TH1F *h = (TH1F *)histos[i - 1];
    legends_str.push_back(realNames[h->GetTitle()]);
  }

  // Add and draw the plots
  for (int i = 0; i < histos.GetEntries(); i++)
  {
    TH1F *h = (TH1F *)histos[i];

    if (normalised)
    {
      double val1 = h->GetSumOfWeights();
      if (fabs(val1) > 0)
        h->Scale(1.0 / val1);
    }

    h->SetLineWidth(1);
    h->SetLineColor(1);
    h->SetFillColor(colours[i]);

    legend->AddEntry(h, legends_str[i].c_str(), "L");
    Hs->Add(h, "sames");
  }

  float heightboxes;
  // the position of the top corner
  float top_corner, deltay;
  // don't draw errors if it is 1
  int no_error = 1;
  top_corner = 0.9;

  // if the array has more than 0 histograms lets specify the stat boxes
  if (histos.GetEntries() > 0)
  {
    heightboxes = (float)0.5 / (float)histos.GetEntries();
  }
  else
    heightboxes = (float)0.5;

  // DRAW not stacked to get correct stats boxes
  if (no_error == 1) // do not draw errors
    Hs->Draw("HIST nostack");
  else
    Hs->Draw("HISTE nostack");

  pad->Update();

  // Work with stats boxes and save copy for later use
  for (int i = histos.GetEntries() - 1; i > -1; i--)
  {
    TH1F *h = (TH1F *)histos.At(n - i - 1);
    if (h != NULL)
    {

      // lets modify the stat boxes
      deltay = i * heightboxes;
      st1 = (TPaveStats *)h->GetListOfFunctions()->FindObject("stats");

      if (st1 != NULL)
      {
        //st1->SetOptStat(1111);
        st1->SetOptFit(0011);
        st1->SetStatFormat("2.3e");
        st1->SetFitFormat("2.3e");
        st1->SetFillColor(19);

        st1->SetY1NDC(top_corner - deltay - heightboxes);
        st1->SetY2NDC(top_corner - deltay);
        st1->SetX1NDC(0.80);
        st1->SetX2NDC(.995);
        st1->SetTextColor(colours[i]);

        // Copy them for later
        TPaveStats *temp = (TPaveStats *)st1->Clone();
        statsboxes.AddLast(temp);
      }
      else
      {
        printf("NULL\n");
      }
    }
  }

  // // limit x axis range
  // if (!std::isnan(maxXAxis))
  // {
  //   Hs->GetXaxis()->SetLimits(0, maxXAxis);
  // }

  pad->Update();
  pad->Modified(); // so it updates the pad with the new changes
  pad->Draw("");

  // DRAW finally stacked
  if (no_error == 1) // do not draw errors
    Hs->Draw("HIST");
  else
    Hs->Draw("HISTE");

  // Redraw correct stats boxes
  for (int i = histos.GetEntries() - 1; i > -1; i--)
  {
    ((TPaveStats *)statsboxes.At(n - i - 1))->Draw();
  }
  legend->Draw("");

  return;
}
