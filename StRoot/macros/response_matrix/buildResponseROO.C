#include <TSystem>

void buildResponseROO()
{
	//load libraries
  // gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	//loadSharedLibraries();
	
	//gSystem->Load("$FASTJETDIR/lib/libfastjet.so");
  // gSystem->Load("$FASTJETDIR/lib/libfastjettools.so");
	//gSystem->Load("$ANALYSISDIR/lib/StRefMultCorr");
	//gSystem->Load("$ANALYSISDIR/lib/StPicoEvent");
	//gSystem->Load("$ANALYSISDIR/lib/StJetAna");

	Float_t R=atof(gSystem->Getenv("RPARAM"));
	Float_t pTlead=atof(gSystem->Getenv("PTTHRESH"));
	int nevents=atoi(gSystem->Getenv("NEVENTS"));
	TString rmatrix_type = gSystem->Getenv("RMATRIX_TYPE");

	//measured distribution histogram
/*
	TString data_path = gSystem->Getenv("DATA_PATH");
	TString str = Form("%s/histos_inclusivejet_R%.1lf.root", data_path.Data(),R);
	TFile *finput = new TFile(str.Data(), "OPEN");

	TH2D *h2 = (TH2D*)finput->Get(Form("hpT_pTlead_R0%.0lf",R*10));
	Int_t firstbin = h2->GetYaxis()->FindBin(pTlead);
	Int_t lastbin = h2->GetNbinsY();
	TH1D *hmeasured = h2->ProjectionX("hprmeasured", firstbin, lastbin);
*/
		StJetBuildResponseMatrixROO *rmatrix = new StJetBuildResponseMatrixROO(R, pTlead, rmatrix_type, nevents/*,hmeasured*/);
		rmatrix->BuildDeltaPtResponseMatrix();
		delete rmatrix;
}
