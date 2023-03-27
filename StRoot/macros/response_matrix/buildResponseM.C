#include <TSystem>

void buildResponseM()
{
	//load libraries
  // gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	//loadSharedLibraries();

	//gSystem->Load("$FASTJETDIR/lib/libfastjet.so");
  //gSystem->Load("$FASTJETDIR/lib/libfastjettools.so");
	//gSystem->Load("$ANALYSISDIR/lib/StRefMultCorr");
	//gSystem->Load("$ANALYSISDIR/lib/StPicoEvent");
	//gSystem->Load("$ANALYSISDIR/lib/StJetAna");


	Float_t R=atof(gSystem->Getenv("RPARAM"));
	Float_t pTleading=atof(gSystem->Getenv("PTLEAD"));
	TString path = gSystem->Getenv("PATH_TO_DELTA_PT_HISTOGRAMS");
	TString RMtype = gSystem->Getenv("RMTYPE");
	
	StJetBuildResponseMatrix *rmatrix = new StJetBuildResponseMatrix(path, R, pTleading,RMtype);
	rmatrix->BuildDeltaPtResponseMatrix();
}
