#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TGaxis.h>
#include <TLatex.h>
#include <string> 
#include "TMath.h"
#include "TGraph.h"

double GetNDC(double x) 
{
  gPad->Update();//this is necessary!
  return (x - gPad->GetX1())/(gPad->GetX2()-gPad->GetX1());
}   