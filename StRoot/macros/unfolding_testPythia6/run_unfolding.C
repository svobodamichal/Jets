#include "unfold.C"

void run_unfolding()
{

  int pTleadMin = atoi(gSystem->Getenv("PTLEAD_MIN"));
  int pTleadMax = atoi(gSystem->Getenv("PTLEAD_MAX"));

  Int_t efficorr= atoi(gSystem->Getenv("EFFICORR"));//correct the result for jet reconstruction efficiency
  TString epsilon_path = gSystem->Getenv("EPSILON_PATH"); //path to efficiency files
  TString systematics = gSystem->Getenv("SYS"); //systematic uncertainty version
  Double_t R = atof(gSystem->Getenv("RPARAM"));
  TString centrality = gSystem->Getenv("SUFFIX"); //centrality 

	//Efficiency correction  
	TH1D* hepsilon;
	if(efficorr){
		TString str = Form("%s/pythia6_effi_%s.root",epsilon_path.Data(),systematics.Data());
		TFile *fepsilon= new TFile(str.Data(), "OPEN");
	}


	  for(int pTlead=pTleadMin;pTlead<=pTleadMax;pTlead++){
	  	if (pTlead == 1 || pTlead == 2 || pTlead == 8 || pTlead == 6) continue;
			if(efficorr) hepsilon=(TH1D*) fepsilon->Get(Form("heffi_pTl%i_R0%.0f_%s",pTlead, R*10,centrality.Data()));
		  unfold(0,(double) pTlead,hepsilon);
	  }
	  
if(efficorr)fepsilon->Close();
}
