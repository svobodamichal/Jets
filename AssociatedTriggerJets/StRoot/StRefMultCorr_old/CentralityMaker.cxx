//----------------------------------------------------------------------------------------------------
// $Id: CentralityMaker.cxx,v 1.5 2015/05/22 06:51:56 hmasui Exp $
// $Log: CentralityMaker.cxx,v $
// Revision 1.5  2015/05/22 06:51:56  hmasui
// Add grefmult for Run14 Au+Au 200 GeV
//
// Revision 1.4  2013/05/10 18:33:33  hmasui
// Add TOF tray mult, preliminary update for Run12 U+U
//
// Revision 1.3  2012/05/19 00:49:16  hmasui
// Update refmult3
//
// Revision 1.2  2012/05/08 03:19:36  hmasui
// Move parameters to Centrality_def_refmult.txt
//
// Revision 1.1  2012/04/23 21:32:12  hmasui
// Interface for future extention of centrality correction maker to other centrality measures, like refmult2
//
//
//----------------------------------------------------------------------------------------------------

#include <iostream>
#include "StRefMultCorr.h"
#include "CentralityMaker.h"

using std::cout ;
using std::endl ;

ClassImp(CentralityMaker)

  CentralityMaker* CentralityMaker::fInstance = 0 ;

//____________________________________________________________________________________________________
CentralityMaker::CentralityMaker()
{
  // Create instance for centrality classes
  fRefMultCorr  = new StRefMultCorr("refmult") ;
  fRefMult2Corr = new StRefMultCorr("refmult2") ;
  fRefMult3Corr = new StRefMultCorr("refmult3") ;
  fTofTrayMultCorr = new StRefMultCorr("toftray") ;
  fgRefMultCorr  = new StRefMultCorr("grefmult") ;
  fgRefMultCorr_P16id  = new StRefMultCorr("grefmult_P16id") ;
  fgRefMultCorr_P18ih  = new StRefMultCorr("grefmult_P18ih") ;
  fgRefMultCorr_VpdMB30 = new StRefMultCorr("grefmult_VpdMB30") ;
  fgRefMultCorr_VpdMBnoVtx = new StRefMultCorr("grefmult_VpdMBnoVtx") ;
}

//____________________________________________________________________________________________________
CentralityMaker::~CentralityMaker()
{ }

//____________________________________________________________________________________________________
CentralityMaker* CentralityMaker::instance()
{
  if ( !fInstance ) {
    // Initialize StRefMultCorr only once
    fInstance = new CentralityMaker() ;
  }

  return fInstance ;
}

//____________________________________________________________________________________________________
StRefMultCorr* CentralityMaker::getRefMultCorr()
{
  return fRefMultCorr ;
}

//____________________________________________________________________________________________________
StRefMultCorr* CentralityMaker::getRefMult2Corr()
{
  return fRefMult2Corr ;
}

//____________________________________________________________________________________________________
StRefMultCorr* CentralityMaker::getRefMult3Corr()
{
  return fRefMult3Corr ;
}

//____________________________________________________________________________________________________
StRefMultCorr* CentralityMaker::getTofTrayMultCorr()
{
  return fTofTrayMultCorr ;
}

//____________________________________________________________________________________________________
StRefMultCorr* CentralityMaker::getgRefMultCorr()
{
  return fgRefMultCorr ;
}

//____________________________________________________________________________________________________
StRefMultCorr* CentralityMaker::getgRefMultCorr_P16id()
{
    return fgRefMultCorr_P16id ;
}

//____________________________________________________________________________________________________
StRefMultCorr* CentralityMaker::getgRefMultCorr_VpdMB30()
{
    return fgRefMultCorr_VpdMB30 ;
}

//____________________________________________________________________________________________________
StRefMultCorr* CentralityMaker::getgRefMultCorr_VpdMBnoVtx()
{
    return fgRefMultCorr_VpdMBnoVtx ;
}

//____________________________________________________________________________________________________
void CentralityMaker::help() const
{
  cout << endl;
  cout << "//------------------------------------------------------------------------------" << endl;
  cout << "How to get centrality bins by CentralityMaker ?" << endl;
  cout << "  (Please also take a look at StRoot/StRefMultCorr/macros/getCentralityBins.C" << endl;
  cout << endl;
  cout << "1. Initialize run index to read proper data base" << endl;
  cout << "  Suppose we want to get centrality from refmult at run index 11078000" << endl;
  cout << endl;
  cout << "  // NOTE:" << endl;
  cout << "  //  Use BBC coincidence rate (NOT ZDC coincidence rate) for refmult2)" << endl;
  cout << "  StRefMultCorr* refmultCorr = CentralityMaker::instance()->getRefMultCorr();" << endl;
  cout << "  refmultCorr->init(11078000);" << endl;
  cout << endl;
  cout << "2. Initialize relevant variables event-by-event" << endl;
  cout << endl;
  cout << "  // NOTE:" << endl;
  cout << "  //  1st argument is original multiplicity" << endl;
  cout << "  //  If one wants to have centrality from refmult2, you have to put refmult2" << endl;
  cout << "  refmultCorr->initEvent(refmult, vz, zdcCoincidenceRate);" << endl;
  cout << endl;
  cout << "3. Get centrality bins" << endl;
  cout << endl;
  cout << "  const Int_t cent16 = refmultCorr->getCentralityBin16() ;" << endl;
  cout << "  const Int_t cent9  = refmultCorr->getCentralityBin9() ;" << endl;
  cout << "//------------------------------------------------------------------------------" << endl;
  cout << endl;
}


