#include "unfold.C"

void run_unfolding()
{

  int pTleadMin = atoi(gSystem->Getenv("PTLEAD_MIN"));
  int pTleadMax = atoi(gSystem->Getenv("PTLEAD_MAX"));

  Int_t efficorr= atoi(gSystem->Getenv("EFFICORR"));//correct the result for jet reconstruction efficiency
  TString epsilon_path = gSystem->Getenv("EPSILON_PATH"); //path to efficiency files
  Double_t R = atof(gSystem->Getenv("RPARAM"));

	//Efficiency correction  
	TH1D* hepsilon;
	if(efficorr){
		TString str = Form("%s/pythia_emb_R%.1lf.root",epsilon_path.Data(),R);
		TFile *fepsilon= new TFile(str.Data(), "OPEN");
	}


	  for(int pTlead=pTleadMin;pTlead<=pTleadMax;pTlead++){
	  	if (pTlead == 8 || pTlead == 6) continue;
			if(efficorr) hepsilon=(TH1D*) fepsilon->Get(Form("hepsilon_pTl%i",pTlead));
		  unfold(0,(double) pTlead,hepsilon);
	  }
	  
if(efficorr)fepsilon->Close();
}
