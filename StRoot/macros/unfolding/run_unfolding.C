#include "unfold.C"

void run_unfolding()
{

  Int_t priorMin= atoi(gSystem->Getenv("PRIOR_MIN"));
  Int_t priorMax= atoi(gSystem->Getenv("PRIOR_MAX"));
  Int_t priorSkip1= atoi(gSystem->Getenv("PRIOR_SKIP1"));
  Int_t priorSkip2= atoi(gSystem->Getenv("PRIOR_SKIP2"));
  Int_t priorSkip3= atoi(gSystem->Getenv("PRIOR_SKIP3"));
  int pTleadMin = atoi(gSystem->Getenv("PTLEAD_MIN"));
  int pTleadMax = atoi(gSystem->Getenv("PTLEAD_MAX"));

  Int_t efficorr= atoi(gSystem->Getenv("EFFICORR"));//correct the result for jet reconstruction efficiency
  TString epsilon_path = gSystem->Getenv("EPSILON_PATH"); //path to efficiency files
  Double_t R = atof(gSystem->Getenv("RPARAM"));
  TString centrality = gSystem->Getenv("SUFFIX"); //centrality 

	//Efficiency correction  
	TH1D* hepsilon;
	if(efficorr){
		//TString str = Form("%s/pythia_emb_R%.1lf.root",epsilon_path.Data(),R);
		TString str = Form("%s/pythia6_effi.root",epsilon_path.Data());
		TFile *fepsilon= new TFile(str.Data(), "OPEN");
	}


  for(int prior=priorMin;prior<=priorMax;prior++){
	  if(prior==priorSkip1 || prior==priorSkip2 || prior==priorSkip3) continue;
	  for(int pTlead=pTleadMin;pTlead<=pTleadMax;pTlead++){
			if (pTlead == 8 || pTlead ==6) continue; //we don't have rm with pTlead 6 or 8
			//if(efficorr) hepsilon=(TH1D*) fepsilon->Get(Form("hepsilon_pTl%i",pTlead));
			if(efficorr) hepsilon=(TH1D*) fepsilon->Get(Form("heffi_pTl%i_R0%.0f_%s",pTlead, R*10,centrality.Data()));
		  unfold(prior,(double) pTlead,hepsilon);
	  }
  }
if(efficorr)fepsilon->Close();
}
