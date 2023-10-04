#include "TF1.h"
#include "TH2D.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TString.h"
#include "TSystem.h"
#include "Riostream.h"
#include "TDirectoryFile.h"
 
#include "StJetBuildResponseMatrix.h"

//static const  Float_t fEmbPt[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 20.0, 40.0, 90.0};
//static const  Float_t fEmbPt[] = {3.0, 5.0, 7.0, 10.0, 20.0};
static const  Float_t fEmbPt[] = {5.0, 7.0, 10.0, 20.0}; //get rid of discontinuity below 5 GeV

ClassImp(StJetBuildResponseMatrix)

//==============================================================================
StJetBuildResponseMatrix::StJetBuildResponseMatrix()
{
  nbins = 800;
  nevts = 1E9;

  pTmin = 0.0;
  pTmax = 100.;

  hResponse = 0x0;
  
  fout = new TFile("response_matrix_gaussian.root", "RECREATE");
}

//=============================================================================
StJetBuildResponseMatrix::StJetBuildResponseMatrix(TString path, Float_t R, Float_t pTleading, TString RMtype)
{
  kCentral=atoi(gSystem->Getenv("CENTRAL"));
  kv2corr=atoi(gSystem->Getenv("V2CORR"));
  v2path=gSystem->Getenv("V2PATH");

  nbins = 800;
  nevts = 1E9;

  //pTmin = pTleading; //changed to 0
	pTmin = 0;
  pTmax = 100.;

  //prepare histograms
  hResponse = 0x0;
  hntrue=new TH1D("hntrue", "n true generated;p_{T}^{true};entries", nbins, -pTmax, +pTmax);
  
  //str=Form("%s/histos_embeddedjet_R%.1lf.root", path.Data(),R);
  str=Form("%s/deltapt_sp_centrality.root", path.Data()); //single particle file
	//cout<<"opening file:"<<str<<endl;
  finput = new TFile(str.Data(), "OPEN");
	//TList* list1 = (TList *)finput->Get("stPicoHFJetMaker"); 

  str=Form("%s/v2corr.root", v2path.Data()); //file with v2 correction factors for delta-pT
  fv2=new TFile(str.Data(), "OPEN");
	float Rv2=R;
	if(R>0.4) Rv2=0.4; //we have corrections only up to R=0.4
  str=Form("dpt_v2corrVSuncorr_R0%.0lf",Rv2*10);
  if(kCentral) str+="_central";
	else str+="_peripheral";
  TH2D* hv2=(TH2D*)fv2->Get(str.Data());
  TString suff="_normal";
  if(kv2corr)suff="_v2";
	  //std::cout << Form("%s/rmatrix%s/response_matrix_%s_R%.1lf_pTlead%.0lf.root", path.Data(), suff.Data(),RMtype.Data(), R, pTleading) << std::endl;
  fout = new TFile(Form("%s/rmatrix%s/response_matrix_%s_R%.1lf_pTlead%.0lf.root", path.Data(), suff.Data(),RMtype.Data(), R, pTleading), "recreate");

//the following lines were inside the emb loop
   //TString name = Form("delta_pt_%s_%.0lf_R0%.0lf", RMtype.Data(),pTleading,R*10);
		TString name = Form("delta_pt_%s_%.0lf_R0%.0lf_central", RMtype.Data(),pTleading,R*10);
    if (!kCentral) name = Form("delta_pt_%s_%.0lf_R0%.0lf_peripheral", RMtype.Data(),pTleading,R*10);
		//cout<<"loading histogram:"<<name.Data()<<endl;
		//TH2D* htmp= (TH2D*)list1->FindObject(name.Data());
		TH2D* htmp= (TH2D*)finput->Get(name.Data());
  for(Int_t idist = 0; idist < nEmb; idist++)
    {
  
		name=Form("dpt_pTl%.0lf_emb%i",pTleading,idist);
//cout<<"htmp->GetName()"<<htmp->GetName()<<endl;
		Int_t firstbin=htmp->GetXaxis()->FindBin(fEmbPt[idist]);
		Int_t lastbin=firstbin;
      hdpT[idist] = (TH1D*)htmp->ProjectionY(name,firstbin,lastbin);

      // removing 1 +/- 1
      for(Int_t bin = 1; bin <= 2000; bin++)
	if(hdpT[idist]->GetBinContent(bin) == hdpT[idist]->GetBinError(bin))
	  {
	    hdpT[idist]->SetBinContent(bin, 0);
	    hdpT[idist]->SetBinError(bin, 0);
	  }

//v2 correction from Alex
		if(kv2corr){
		hdpT[idist]->Write(Form("hdeltapT_%.0lf_v2uncorr",fEmbPt[idist]));
      for(Int_t bin = 1; bin <= 2000; bin++)
		{
			Float_t pTreco=hdpT[idist]->GetBinCenter(bin);
			Float_t prob=hdpT[idist]->GetBinContent(bin);
			if(!prob>0)continue;
			Float_t prob_err=hdpT[idist]->GetBinError(bin);
			Int_t binx=hv2->GetXaxis()->FindBin(fEmbPt[idist]);
			Int_t biny=hv2->GetYaxis()->FindBin(pTreco);			
			Float_t scaler=hv2->GetBinContent(binx,biny);
			if(!scaler>0)scaler=1;
			hdpT[idist]->SetBinContent(bin,prob*scaler);
			hdpT[idist]->SetBinError(bin,prob_err*scaler);
		}
		}//v2 corr
	hdpT[idist]->Write(Form("hdeltapT_%.0lf",fEmbPt[idist]));
    }//pTemb loop

		delete htmp;

}

//==============================================================================
StJetBuildResponseMatrix::~StJetBuildResponseMatrix()
{
  fout->Close();
  delete fout;
}

//==============================================================================
void StJetBuildResponseMatrix::BuildGaussianResponseMatrix()
{
  TString name;
  Double_t sigma = 4.;
  Int_t save_evt = 1;

  hResponse = new TH2D("hResponse", "hResponse;p_{T}^{meas};p_{T}^{true};entries", nbins, -pTmax, +pTmax, nbins, -pTmax, +pTmax);

  for(Int_t ievt = 0; ievt < nevts; ievt++)
    {
      Double_t pT = gRandom->Uniform(pTmin, pTmax);
      Double_t dpT = gRandom->Gaus(0, sigma);
      hResponse->Fill(pT + dpT, pT);

      if(ievt != TMath::Power(10, save_evt) - 1) continue;
      name = Form("hResponse_1E%d", save_evt);
      fout->cd();
      hResponse->Write(name.Data());
      cout << Form("%s saved!", name.Data()) << endl;
      save_evt++;
    }

  fout->cd();
  hResponse->Write();
  delete hResponse;
}

//==============================================================================
void StJetBuildResponseMatrix::BuildDeltaPtResponseMatrix()
{
  TString name;
  Int_t save_evt = 0;

  hResponse = new TH2D("hResponse", "hResponse;p_{T}^{meas};p_{T}^{true};entries", nbins, -pTmax, +pTmax, nbins, -pTmax, +pTmax);
  gRandom = new TRandom3(0);

  for(Int_t ievt = 0; ievt < nevts; ievt++)
    {
      Double_t pT = gRandom->Uniform(pTmin, pTmax);
      Double_t dpT = SmearWithDeltaPt(pT);
      hResponse->Fill(pT + dpT, pT);
		hntrue->Fill(pT);

		if(ievt+1==TMath::Power(10, save_evt) )
		{
			cout <<"doing event 1E"<<save_evt<<endl;
			save_evt++;
		}
/*
      if(ievt != TMath::Power(10, save_evt) - 1) continue;
      name = Form("hResponse_1E%d", save_evt);
      fout->cd();
      hResponse->Write(name.Data());
      cout << Form("%s saved!", name.Data()) << endl;
      save_evt++;*/
    }//event loop

//normalize matrix
for(int j=1; j<nbins+1; j++){
      double ntot=hntrue->GetBinContent(j);
   for(int i=1; i<nbins+1; i++){
      double binc = hResponse->GetBinContent(i,j);
		if(binc<1 || ntot<1)continue;
		double prob=binc/ntot;
      hResponse->SetBinContent(i,j,prob);
   }//j - true
}//i - measured

  fout->cd();
  hResponse->Write("hResponse_1E9");
  hntrue->Write("hntrue");
 
  delete hResponse;
}

//==============================================================================
Double_t StJetBuildResponseMatrix::SmearWithDeltaPt(Double_t pT)
{
   Double_t dpT = 0;
  if (pT <= fEmbPt[0]){
    dpT = hdpT[0]->GetRandom();
   }
  for(int i=0; i<(nEmb-1); i++){
   if(pT>fEmbPt[i] && pT<=fEmbPt[i+1]){
   Double_t rnd=gRandom->Uniform(fEmbPt[i],fEmbPt[i+1]);
   if(pT>rnd) dpT = hdpT[i+1]->GetRandom();
   else dpT = hdpT[i]->GetRandom();
   }
  }
  if(pT>fEmbPt[nEmb-1]) {
   dpT = hdpT[nEmb-1]->GetRandom();
}
  return dpT;

/*
//-------------------------
  Double_t dpT = 0.0;
	  
  if(pT < 1.0)
    dpT = hdpT[0]->GetRandom();
  else
  if(pT < 2.0)
    dpT = hdpT[1]->GetRandom();
  else
  if(pT < 3.0)
    dpT = hdpT[2]->GetRandom();
  else
  if(pT < 5.0)
    dpT = hdpT[3]->GetRandom();
  else
  if(pT < 10.0)
	 dpT = hdpT[4]->GetRandom();
  else
  if(pT < 15.0)
    dpT = hdpT[5]->GetRandom();
  else
    dpT = hdpT[6]->GetRandom();
  return dpT;
  */
}
