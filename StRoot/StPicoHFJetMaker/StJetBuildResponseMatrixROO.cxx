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

#include "StJetBuildResponseMatrixROO.h"

/*

	REMEMBER, MODIFIED FOR PRIOR TESTING

*/

ClassImp(StJetBuildResponseMatrixROO)

	TString prior_type[]={"","flat","pythia","powlaw4","powlaw45","powlaw5","powlaw55","tsalis_1","tsalis_2","tsalis_3","tsalis_4","tsalis_5","tsalis_6","tsalis_7","tsalis_8","tsalis_9"};
	//TString prior_type[]={"","mtsalis_1","mtsalis_2","mtsalis_3","mtsalis_4","mtsalis_5","mtsalis_6", "gauspol", "expol", "mpowlaw4","mpowlaw45","mpowlaw5","mpowlaw55","mpowlaw6","gauspol_1","gauspol_2","gauspol_3","gauspol_4","gauspol_5"}; //TEST - modified tsalis functions	
	//TString prior_type[]={"","mtsalis_1","mtsalis_2","mtsalis_3","mtsalis_4","mtsalis_5","mtsalis_6", "gauspol", "expol", "gammapol_1","gammapol_2","gammapol_3","gammapol_4","gammapol_5","gauspol_1","gauspol_2","gauspol_3","gauspol_4","gauspol_5"}; //TEST - modified tsalis functions	
	//=============================================================================
StJetBuildResponseMatrixROO::StJetBuildResponseMatrixROO(Float_t R, Float_t pTthresh, TString type, Int_t nEvents, TH1D* hpriorhisto)
{
	path_in=gSystem->Getenv("RM_PATH");
	path_out=gSystem->Getenv("OUT_PATH");
	//path_out=Form("%s/ROOformat/", path_in.Data());
	pyEmb_path = gSystem->Getenv("PYEMB_PATH");
	prior_path= gSystem->Getenv("PRIOR_PATH");
	priorStart=atoi(gSystem->Getenv("PRIOR_START"));
	priorStop=atoi(gSystem->Getenv("PRIOR_STOP"));
	priorSkip1=atoi(gSystem->Getenv("PRIOR_SKIP1"));
	priorSkip2=atoi(gSystem->Getenv("PRIOR_SKIP2"));
	priorSkip3=atoi(gSystem->Getenv("PRIOR_SKIP3"));
	TString centrality = gSystem->Getenv("CENTRALITY");

	pTlead=pTthresh;

	nbins = 800;
	nevts = nEvents;

	pTmax = 100.;
	pTmin = -pTmax;
	if(type=="effi"){
		pTmin=0;
		nbins=400;
	}

	//hResponse = 0x0;
	mtype=type;


	//for priorNo=0 take prior from input histogram
	if(priorStart==0) hprior0=hpriorhisto;

	//TString prior_type[]={"","flat","pythia","powlaw3","powlaw45","powlaw5","powlaw55","levy_jan","levy_alex"};

	//str=Form("%s/response_matrix_%s_R%.1lf_pTlead%.0lf.root", path_in.Data(),mtype.Data(),R,pTthresh); //input matrix file classical
	//str=Form("%s/pythia6_normalized_combined_response_embed.root", path_in.Data()); //input matrix file embedding
	str=Form("%s/pythia6_normalized_combined_response_HT2.root", path_in.Data()); //input matrix file embedding with HT2 cut
	if(mtype=="effi" || mtype=="dete")str=Form("%s/pythia_emb_R%.1lf.root",pyEmb_path.Data(),R);
	finput = new TFile(str, "OPEN");
	cout<<"reading input file "<<str.Data()<<endl;

//	TString name = "hResponse_1E9"; //use "classical" RM as input
	TString name = Form("hResponseMatrix_pTl%.0lf_R0%.0f_%s", pTthresh, R*10, centrality.Data()); //used for embedding RM
	if(mtype=="effi" || mtype=="dete") name = Form("hresponse_pTl%.0lf",pTthresh);
 	//cout << "getting histogram " << name.Data() << endl; 
	hRMin = (TH2D*)finput->Get(name.Data());


	//get delta-pT histograms for all positive probe pT values
	Int_t binzero=hRMin->GetYaxis()->FindBin(0.0);
	for(Int_t bin=binzero; bin<hRMin->GetNbinsY()+1; bin++){
		TString name=Form("hdpT_%i",bin);
		hdpT[bin-binzero]=(TH1D*)hRMin->ProjectionX(name, bin, bin);
		//cout<<"bin:"<<bin<<", pT:"<<hRMin->GetYaxis()->GetBinCenter(bin)<<", integral:"<<hdpT[bin-binzero]->Integral()<<endl;
	}

	//create output
	//str = Form("%s/response_matrix_%s_R%.1lf_pTlead%.1lf_allPriors.root", path_out.Data(),mtype.Data(),R,pTthresh); //classic RM
	str = Form("%s/response_matrix_emb_R%.1lf_pTlead%.1lf_allPriors.root", path_out.Data(),R,pTthresh); //embedding RM
	
	//str = Form("%s/response_matrix_emb_R%.1lf_pTlead%.1lf_newPriors.root", path_out.Data(),R,pTthresh); //TEST - modified tsalis priors	
	
	if(mtype=="effi")str= Form("%s/rmatrix/response_matrix_%s_R%.1lf_pTlead%.1lf_allPriors.root", pyEmb_path.Data(),mtype.Data(),R,pTthresh);
	fout = new TFile(str, "RECREATE");


	//LOAD PRIOR FUNCTION
	str = Form("%s/priors_default.root", prior_path.Data());
	//str = Form("%s/priors_mtsalis.root", prior_path.Data()); //TEST - modified tsalis priors	
	TFile *fpriorfile = new TFile(str.Data(), "OPEN");
	cout<<"reading prior file: "<<str.Data()<<endl;

	for (int priorNo=priorStart; priorNo<=priorStop; priorNo++)
	{	
		if(priorNo==priorSkip1 || priorNo==priorSkip2 || priorNo==priorSkip3) continue;
		str=Form("%s",prior_type[priorNo].Data());
		cout << "priorL: " << str.Data() << endl;
		if(priorNo==2) str=Form("%s_R%.0lf_pTlead%.0lf",prior_type[priorNo].Data(),R*10,pTthresh); //TEST - UNCOMMENT FOR CORRECT PYTHIA USE AFTERWARDS
		//if(priorNo==2) str=Form("%s_R%.0lf_pTlead5",prior_type[priorNo].Data(),R*10);
		fprior[priorNo] = (TF1*)fpriorfile->Get(str.Data());
		
		cout << "prior: " << fprior[priorNo]->GetName() << endl;
	}


}

//==============================================================================
StJetBuildResponseMatrixROO::~StJetBuildResponseMatrixROO()
{
	fout->Close();
	delete fout;
}

//==============================================================================
void StJetBuildResponseMatrixROO::BuildDeltaPtResponseMatrix()
{
	TString name;
	TString title;

	for (int prior=priorStart; prior<=priorStop; prior++)
	{
		if(prior==priorSkip1 || prior==priorSkip2 || prior==priorSkip3) continue;
		name=Form("hResponse_%s",prior_type[prior].Data());
		title=Form("hResponse_%s; p_{T}^{meas};p_{T}^{true};entries",prior_type[prior].Data());
		hResponse[prior] = new TH2D(name,title, nbins, pTmin, pTmax, nbins, pTmin, pTmax);
		name=Form("hMCtrue_%s",prior_type[prior].Data());
		title=Form("hMCtrue_%s; p_{T}^{true};entries",prior_type[prior].Data());
		hMCtrue[prior] = new TH1D(name,title, nbins, pTmin, pTmax);
		name=Form("hMCreco_%s",prior_type[prior].Data());
		title=Form("hMCreco_%s; p_{T}^{meas};entries",prior_type[prior].Data());
		hMCreco[prior] = new TH1D(name, title, nbins, pTmin, pTmax);
		hResponse[prior]->Sumw2();
		hMCtrue[prior]->Sumw2();
		hMCreco[prior]->Sumw2();
	}
	gRandom = new TRandom3(0);
	for(Int_t ievt = 0; ievt < nevts; ievt++){
		if(ievt%1000000==0)	cout<<"EVENT:"<<ievt<<endl;
		Double_t pT=gRandom->Uniform(0,pTmax); 
		//if(pT<=pTlead) continue; //DON'T USE THIS CUT NOW
		//get random pT-measured value using delta-pT distribution for given pT-generated
		Double_t dpT = SmearWithDeltaPt(pT);
		if(TMath::Abs(dpT)<1E-9) continue; //dpT=0 - dpT is set to 0 for invalid pTvalues
		//fill values for all priors
		for (int prior=priorStart; prior<=priorStop; prior++)
		{
			if(prior==priorSkip1 || prior==priorSkip2 || prior==priorSkip3) continue;
			Double_t scale=0;
			if(prior==0)
			{
				cout<<"creating prior function from input histogram"<<endl;
				int bn=hprior0->FindBin(pT);
				float yield=hprior0->GetBinContent(bn);
				float intg=hprior0->Integral();
				if(intg>0)scale=yield/intg;		
			}
			else	
			{
				scale = fprior[prior]->Eval(pT);
			}
			if(scale!=scale) continue; //fprior is not defined for this pT
			//weight the output by the prior distribution - necessary for ROOUnfold
			hMCtrue[prior]->Fill(pT,scale);
			hResponse[prior]->Fill(dpT, pT,scale);
			hMCreco[prior]->Fill(dpT, scale);
		}//prior loop
	}//event loop

	//save output histograms
	fout->cd();
	for (int prior=priorStart; prior<=priorStop; prior++)
	{
		if(prior==priorSkip1 || prior==priorSkip2 || prior==priorSkip3) continue;
		name = Form("hResponse_%s", prior_type[prior].Data());
		hResponse[prior]->Write(name.Data());
		name = Form("hMCtrue_%s", prior_type[prior].Data());
		hMCtrue[prior]->Write(name.Data());
		name = Form("hMCreco_%s", prior_type[prior].Data());
		hMCreco[prior]->Write(name.Data());
		cout << "All events saved!"<< endl;

		delete hResponse[prior];
		delete hMCtrue[prior];
		delete hMCreco[prior];
	}//prior loop
	delete gRandom;
}

//==============================================================================
Double_t StJetBuildResponseMatrixROO::SmearWithDeltaPt(Double_t pT)
{
	if(pT>pTmax-0.001)pT=pTmax-0.001;  //to avoid pT=pTmax
	Double_t dpT = 0;
	Int_t  bin = hRMin->GetYaxis()->FindBin(pT);
	Int_t binzero=hRMin->GetYaxis()->FindBin(0.0);
	float integral=hdpT[bin-binzero]->Integral();
	if (pT < 0 || integral<1E-9)dpT=0; //the input response matrix does not have values for this pT
	else dpT = hdpT[bin-binzero]->GetRandom();
	return dpT;

}
