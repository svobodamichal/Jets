#ifndef __StJetBuildResponseMatrix__hh
#define __StJetBuildResponseMatrix__hh

#include "TString.h"

class TH1D;
class TH2D;
class TFile;

class StJetBuildResponseMatrix
{
 public:
  StJetBuildResponseMatrix();
  StJetBuildResponseMatrix(TString path, Float_t R, Float_t pTleading, TString RMtype="BG_sp");
  virtual ~StJetBuildResponseMatrix();
  
  void BuildDeltaPtResponseMatrix();
  void BuildGaussianResponseMatrix();
  
 private:
  Double_t SmearWithDeltaPt(Double_t pT);

 protected:
  TFile *fout;
  TFile *finput;
  TFile *fv2;

  TH1D *hntrue;
  TH1D *hdpT[20];
  TH2D *hResponse;

  Bool_t kCentral;
  Bool_t kv2corr;
	
  Int_t nbins;
  Int_t nevts;
  
  Double_t pTmin;
  Double_t pTmax;

  TString str;
  TString v2path;

  //static const Int_t nEmb=12;
  //static const Int_t nEmb=5;
	static const Int_t nEmb=4;

  ClassDef(StJetBuildResponseMatrix, 1)
};

#endif
