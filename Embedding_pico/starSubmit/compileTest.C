#include "TString.h"
#include "TROOT.h"

void compileTest(const char* name = "runPicoHFJetMaker.C"){ //orig. runPicoDpmAnaMaker.C
  Long_t ret= gROOT->ProcessLine(".L StRoot/macros/loadSharedHFLibraries.C");
  cout << ret << endl;
  loadSharedHFLibraries();
  ret = gROOT->ProcessLine(Form(".L StRoot/macros/%s++", name));
}
