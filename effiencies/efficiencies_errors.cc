/*                             __
                           .--()°"."
 Author: Nathalia Cardona "|, . ,"
                           !_-(_\
*/

#include "PhenoAnalyzer.h"
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
// #include "error_bayesian.h"
// #include "error_binomial.h"

using namespace std;

float calcBinEffError(float numerator, float denominator)
{
  if (denominator > 0)
  {
    float efficiency = numerator / denominator;
    float efferror = sqrt(efficiency * (1.0 - efficiency) / denominator);
    return efferror;
  }
  return -1;
}

float calcBayEffError(float numerator, float denominator)
{
  float effError;

  TH1F *theNumHisto = new TH1F("theNumHisto", "theNumHisto", 1, 0, 1);
  theNumHisto->SetBinContent(1, numerator);
  theNumHisto->Sumw2();

  TH1F *theDenHisto = new TH1F("theDenHisto", "theDenHisto", 1, 0, 1);
  theDenHisto->SetBinContent(1, denominator);
  theDenHisto->Sumw2();

  TGraphAsymmErrors *bayesEff = new TGraphAsymmErrors();
  bayesEff->BayesDivide(theNumHisto, theDenHisto, "b");

  if (bayesEff->GetErrorYhigh(0) > bayesEff->GetErrorYlow(0))
  {
    effError = bayesEff->GetErrorYhigh(0);
  }
  else
  {
    effError = bayesEff->GetErrorYlow(0);
  }

  delete theNumHisto;
  delete theDenHisto;
  delete bayesEff;

  return effError;
}

int main(int argc, char *argv[])
{

  if (argc != 3)
  {

    cout << "Wrong number of parameters. Pleas use ./exec inputFile outputFile." << endl;
    return -1;
  }

  ifstream theFile(argv[1]);
  ofstream outputFile(argv[2]);

  if (!theFile.is_open())
    throw runtime_error("Could not open file");

  string line, colname;

  float val;

  vector<float> data;

  // Read the column names
  if (theFile.good())
  {
    // Extract the first line in the file
    getline(theFile, line);

    // Create a stringstream from line
    stringstream ss(line);

    cout << "headers: " << line << endl;
  }

  int dataIndex = 0;
  // Read data, line by line
  while (getline(theFile, line))
  {
    // Create a stringstream of the current line
    stringstream ss(line);
    data.clear();
    // cout << "the line: " << line << endl;

    // Extract each integer
    while (ss >> val)
    {

      // // Add the current integer to the 'colIdx' column's values vector
      // result.at(colIdx).second.push_back(val);

      data.push_back(val);

      // cout << val << " ";
      // If the next token is a comma, ignore it and move on
      if (ss.peek() == ',')
        ss.ignore();
    }

    // process data
    for (int i = 0; (unsigned)i < data.size() - 1; i++)
    {
      float numerator = data[i + 1];
      float denominator = i == data.size() - 2 ? data[i - 1] : data[i];

      float error;
      float eff = numerator / denominator;

      if (eff > 1)
        throw runtime_error("Efficiency is greater than 1. Check the order of the columns.");

      if (eff < 0.25 || eff > 0.75)
      {
        error = calcBayEffError(numerator, denominator);
      }
      else
      {
        error = calcBinEffError(numerator, denominator);
      }
      // cout << "num,denom,eff,error " << numerator << " " << denominator << " " << eff << " " << error << endl;

      outputFile << error;
      if (i != data.size() - 2)
      {
        outputFile << ",";
      }
    }

    outputFile << endl;

    // cout << endl;
    dataIndex++;
  }

  outputFile.close();
  return 0;
}
