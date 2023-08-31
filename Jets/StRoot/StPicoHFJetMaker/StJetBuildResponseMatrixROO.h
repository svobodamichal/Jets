#ifndef __StJetBuildResponseMatrixROO__hh
#define __StJetBuildResponseMatrixROO__hh

#include "TString.h"

class TF1;
class TH1D;
class TH2D;
class TFile;

class StJetBuildResponseMatrixROO
{
 public:
  StJetBuildResponseMatrixROO(Float_t R, Float_t pTthresh, TString type, Int_t nEvents=1E9, TH1D* hpriorhisto=NULL);
  virtual ~StJetBuildResponseMatrixROO();
  
  void BuildDeltaPtResponseMatrix();
  
 private:
  Double_t SmearWithDeltaPt(Double_t pT);

 protected:
  TFile *fout;
  TFile *finput;
  TFile *fepsilon;

  TH2D *hRMin;
  TH1D *hdpT[500];
  TH2D *hResponse[20];
  TH1D *hMCtrue[20];
  TH1D *hMCreco[20];
  TH1D *hprior0;
  TF1 *fprior[20];
  TF1 *fpythia;

  /*
  TH1D *hepsilon;
  TH1D *hMCreco_eff;
  TH2D *hResponse_eff;
  TH2D *hResponse_15;
  TH2D *hResponse_15_eff;
  TH1D *hMCtrue_15;
  TH1D *hMCreco_15;
  TH1D *hMCreco_15_eff;
*/

  Int_t nbins;
  Int_t nevts;
  Int_t priorStart;
  Int_t priorStop;
  Int_t priorSkip1;
  Int_t priorSkip2;
  Int_t priorSkip3;
  
  Double_t pTmin;
  Double_t pTmax;
  Double_t pTcutoff;
  Double_t pTlead;

  TString mtype;
  TString str;
  TString path_in;
  TString path_out;
  TString prior_path;
  TString pyEmb_path;

  ClassDef(StJetBuildResponseMatrixROO, 1)
};

#endif
