#include "StPicoHFJetMaker.h"
#include "JetInfo.h"

#include "StEmcADCtoEMaker/StBemcData.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StEmcRawMaker/StBemcRaw.h"
#include "StEmcRawMaker/defines.h"
#include "StEmcRawMaker/StBemcTables.h"

#include "BemcNewCalib.h"

//FastJet 3
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <../../fastjet/config.h>
#include <../../fastjet/PseudoJet.hh>
#include <../../fastjet/JetDefinition.hh>
#include <../../fastjet/ClusterSequence.hh>
#include <../../fastjet/ClusterSequenceArea.hh>
//#ifdef FASTJET_VERSION
#include <../../fastjet/Selector.hh>
#include <../../fastjet/tools/Subtractor.hh>
#include <../../fastjet/tools/JetMedianBackgroundEstimator.hh>
//#include </gpfs01/star/pwg/licenrob/jets/fastjet/contrib/SoftDrop.hh>	//robotmon 
//#include </gpfs01/star/pwg/licenrob/jets/fastjet/contrib/RecursiveSoftDrop.hh>
//#include <fastjet/internal/base.hh> //robotmon
//#endif


using namespace std;  //robotmon
using namespace fastjet;
//using namespace contrib;	//robotmon

ClassImp(StPicoHFJetMaker)

// was
//
//StPicoHFJetMaker::StPicoHFJetMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName) : StPicoJetMaker(name, picoMaker, outputBaseFileName) {

//}

//jana chagned to

StPicoHFJetMaker::StPicoHFJetMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName, char const* inputHFListHFtree = "") :
StPicoJetMaker(name, picoMaker, outputBaseFileName, inputHFListHFtree),  mRefmultCorrUtil(NULL) {


  // constructor
}

// _________________________________________________________
StPicoHFJetMaker::~StPicoHFJetMaker() {
  // destructor
}

// _________________________________________________________
int StPicoHFJetMaker::InitJets() {


  mADCtoEMaker = dynamic_cast<StEmcADCtoEMaker*>(GetMaker("Eread")); 
	assert(mADCtoEMaker);
  mTables = mADCtoEMaker->getBemcData()->getTables();

  // -- INITIALIZE USER HISTOGRAMS ETC HERE -------------------
  //    add them to the output list mOutList which is automatically written

  // EXAMPLE //  mOutList->Add(new TH1F(...));
  // EXAMPLE //  TH1F* hist = static_cast<TH1F*>(mOutList->Last());

    int npTleadbins = 25;
    float pTleadmin = 0;
    float pTleadmax = 25;
    int nptbins = 1000;
    float ptminbin = -40;
    float ptmaxbin = 60;
    int netabins = 100*2;
    float etaminbin = -1;
    float etamaxbin = 1;
    int nphibins = 120;
    float phiminbin = 0;//-TMath::Pi();
    float phimaxbin = TMath::TwoPi();
    int zbins = 50;
    float zmin = -50;
    float zmax = 50;
    int npttrackbins = 120;
    float pttrackmin = 0;
    float pttrackmax = 30;
		
    int ncetabins = 100*2;
    float cetaminbin = -1;
    float cetamaxbin = 1;
    int ncphibins = 120;
    float cphiminbin = 0;//-TMath::Pi();
    float cphimaxbin = TMath::TwoPi();

		int nptembbins = 200;
    float ptembminbin = 0;
    float ptembmaxbin = 25;
    float deltaptembminbin = -30;
    float deltaptembmaxbin = 50;
    
    TH1::SetDefaultSumw2();

	mOutList->Add(new TH1D("hweight", "weight", 135, 0, 3));
	mOutList->Add(new TH1D("hcent", "centrality", 10, -1, 9));
		
	//mOutList->Add(new TH2D("hrunIdcent", "runId vs centrality", 90913, 15076101, 15167014, 10, -1, 9)); //not used

    //General track QA
	mOutList->Add(new TH2D("heta_phi_tr", "track phi vs. eta;#phi [-];#eta [-]", nphibins, phiminbin, phimaxbin, netabins, etaminbin, etamaxbin));
	mOutList->Add(new TH1D("hphi_tr", "track phi;#phi", nphibins, phiminbin, phimaxbin));
	mOutList->Add(new TH1D("heta_tr", "track eta;#eta", netabins, etaminbin, etamaxbin));    
	mOutList->Add(new TH1D("hpT_tr", "track pT; p_{T} [GeV/c]", npttrackbins, pttrackmin, pttrackmax));
	//	mOutList->Add(new TH1D("hE_tr", "track E from BEMC; E [GeV]", nEtrackbins, Etrackmin, Etrackmax));
	//	mOutList->Add(new TH1D("hpTBEMC_tr", "track pT from BEMC; p_{T} [GeV/c]", nptBEMCtrackbins, ptBEMCtrackmin, ptBEMCtrackmax));
	mOutList->Add(new TH2D("hdca_z_tr", "track DCA vs z-vertex", 90, 0, 3, zbins , zmin, zmax));
	mOutList->Add(new TH1D("hdca_tr", "track DCA", 90, 0, 3));
	mOutList->Add(new TH2D("hdca_pT", "track DCA vs. p_{T}", 90, 0, 3, npttrackbins, pttrackmin, pttrackmax));
	mOutList->Add(new TH1D("hcharged_tr", "track charge", 90, 0, 3));	

		//mOutList->Add(new TH1D("hdca_XY_tr", "global track DCA_XY; DCA_{xy} [cm]", 200, -20, 20));	
    //mOutList->Add(new TH1D("hdca_Z_tr", "global track DCA_Z; DCA_{z} [cm]", 200, -20, 20));	
    //mOutList->Add(new TH2D("hdca_XYZ_tr", "global track DCA_Z vs DCA_XY; DCA_{z} [cm]; DCA_{xy} [cm]", 200, -20, 20, 200, -20, 20));	
   	//for (int lumi = 0; lumi < 5; lumi++) 
		//mOutList->Add(new TH2D(Form("hdca_XYZ_tr_lumi%i", lumi), Form("global track DCA_Z vs DCA_XY, %i < lumi < %i kHz; DCA_{z} [cm]; DCA_{xy} [cm]", lumi*15, (lumi+1)*15), 200, -20, 20, 200, -20, 20));	


		//mOutList->Add(new TH1D("hEcorr", "E from BTow after correction; E [GeV]", 340, -4.5, 29.5));

    //mOutList->Add(new TH2D("hE_tow", "energy vs tower ID;tower ID;E [GeV]", 4800, 1, 4800, 300, -0.5, 29.5));
    //mOutList->Add(new TH2D("hp_tow", "track momentum vs tower ID;tower ID;p [GeV/c]", 4800, 1, 4800, 300, -0.5, 29.5));
  //mOutList->Add(new TH2D("heta_phi_tow", "tower eta vs phi; #eta [-]; #phi [-]", netabins, etaminbin, etamaxbin, nphibins, phiminbin, phimaxbin));		
	  mOutList->Add(new TH2D("heta_phi_tow", "tower eta vs phi; #eta [-]; #phi [-]", netabins/5, etaminbin, etamaxbin, nphibins, phiminbin, phimaxbin));		

    if (isMcMode()) {cout << "MC mode activated" << endl;
        //TODO: Add required histograms for MC Mode

				for(unsigned int r=0; r<fR.size(); r++)
				{
				TString hname=Form("hjetpTembArea_R0%.0lf",fR[r]*10);
				mOutList->Add(new TH2D(hname,"jet pTemb vs area; A [-]; p_{T} [GeV/c]", 100, 0, 1, nptembbins, ptembminbin, ptembmaxbin));
				TString htitle = "delta pT for BG corrections, using sp probe; p_{T}^{emb} [GeV/c]; #delta p_{T} [GeV/c]";
				for (int centbin = -1; centbin < 9; centbin++) {
				for(Int_t pTlcut = 0; pTlcut<npTlead; pTlcut++)
					{
						hname=Form("delta_pt_BG_sp_%i_R0%.0lf_centbin%i", pTlcut, fR[r]*10, centbin);
						if (kEmbPythia == 1) {htitle = "delta pT for BG corrections, using pythia probe; p_{T}^{emb} [GeV/c]; #delta p_{T} [GeV/c]"; hname=Form("delta_pt_BG_py_%i_R0%.0lf_centbin%i", pTlcut, fR[r]*10, centbin);}
						//cout << hname << endl;
 	   				mOutList->Add(new TH2D(hname, htitle, nptembbins, ptembminbin, ptembmaxbin, nptembbins, deltaptembminbin, deltaptembmaxbin));
					}	
				}
				}


		//initialize TPythia6 for pythia jet embedding
		//Setting random seed
		TDatime dt;
	  UInt_t curtime = dt.Get();
		UInt_t procid = gSystem->GetPid();
  	UInt_t seed = curtime - procid;
  	gRandom->SetSeed(seed);

		//PYTHIA
		fpythia = new TPythia6();
		fpythia->SetMRPY(1, seed);

    }

   
						TString hname ="hrho";
            mOutList->Add(new TH1D(hname,"median energy density; #rho_{ch} [GeV/sr]",50,0,50));
					
            hname = "hrho_mult";
            mOutList->Add(new TH2D(hname,"median energy density vs multiplicity; RefMult [-] ; #rho_{ch} [GeV/sr]",600, 0, 600, 50,0,50));

            hname = "hfrho";
            mOutList->Add(new TH1D(hname,"median energy density; #rho [GeV/sr]",100,0,100));

            hname = "hfrho_mult";
            mOutList->Add(new TH2D(hname,"median energy density vs multiplicity ; RefMult [-] ; #rho [GeV/sr]",600, 0, 600, 100,0,100));

 			if (!isMcMode()) {
        for(unsigned int r = 0; r < fR.size(); r++) {
            //TString hname = Form("hpT_pTlead_R0%.0lf",fR[r]*10);
            //mOutList->Add(new TH2D(hname, "jet pTcorr vs pTleading; p_{T} [GeV/c]; p_{T}^{lead} [GeV/c]", nptbins, ptminbin, ptmaxbin, npTleadbins, pTleadmin, pTleadmax));
	    			hname = Form("hpT_R0%.0lf",fR[r]*10);
	    			mOutList->Add(new TH1D(hname, "pT; p_{T} [GeV/c]", nptbins, ptminbin, ptmaxbin)); 

            //hname = Form("heta_phi_R0%.0lf",fR[r]*10);
            //mOutList->Add(new TH2D(hname, "jet eta vs phi;#eta;#phi", netabins, etaminbin, etamaxbin, nphibins, phiminbin, phimaxbin));

	    			hname = Form("heta_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH1D(hname, "jet eta;#eta", netabins, etaminbin, etamaxbin));

	    			hname = Form("hphi_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH1D(hname, "jet phi;#phi", nphibins, phiminbin, phimaxbin));

            //hname = Form("hjetarea_cut_R0%.0lf",fR[r]*10);
            //mOutList->Add(new TH1D(hname,"jet area after cut", 100, 0, 1));

            hname = Form("hjetarea_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH1D(hname,"jet area", 100, 0, 1));

            hname = Form("hjetpTarea_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH2D(hname,"jet pTmeasured vs area; A [-]; p_{T} [GeV/c]" ,100, 0, 1, nptbins, ptminbin, ptmaxbin));

            hname = Form("hjetpTcorrArea_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH2D(hname,"jet pTreco vs area; A [-]; p_{T}^{corr} [GeV/c]",100, 0, 1, nptbins, ptminbin, ptmaxbin));

            //hname = Form("hjetstructure_R0%.0lf",fR[r]*10);
            //mOutList->Add(new TH2D(hname,"jet constituents pT vs pTlead ; p_{T}^{part} [GeV/c]; p_{T}^{lead} [GeV/c]", 50, 0, 25, npTleadbins, pTleadmin, pTleadmax));

            hname = Form("hnparticlesinjet_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH2D(hname,"#particles in jet vs jet pT; # of particles; p_{T}^{lead} [GeV/c]", 300*fR[r], 0, 300*fR[r], npTleadbins, pTleadmin, pTleadmax));

				    hname = Form("hceta_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH1D(hname, "jet constituents eta;#eta", ncetabins, cetaminbin, cetamaxbin));

				    hname = Form("hcphi_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH1D(hname, "jet constituents phi;#phi", ncphibins, cphiminbin, cphimaxbin));

						//full jet histos
 						hname = Form("hfceta_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH1D(hname, "full jet constituents eta;#eta", ncetabins, cetaminbin, cetamaxbin));

				    hname = Form("hfcphi_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH1D(hname, "full jet constituents phi;#phi", ncphibins, cphiminbin, cphimaxbin));	

 						hname = Form("hfjetarea_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH1D(hname,"full jet area", 100, 0, 1));

            hname = Form("hfjetpTarea_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH2D(hname,"full jet pTmeasured vs area; A [-]; p_{T} [GeV/c]", 100, 0, 1, nptbins, ptminbin, ptmaxbin));

            hname = Form("hfjetpTcorrArea_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH2D(hname,"full jet pTreco vs area; A [-]; p_{T}^{corr} [GeV/c]", 100, 0, 1, nptbins, ptminbin, ptmaxbin));

						//hname = Form("hfrho_rho_R0%.0lf",fR[r]*10);
            //mOutList->Add(new TH2D(hname,"mean energy density charged vs full ; #rho ; #rho_f ",100, 0, 100, 100,0,100));

            //hname = Form("hfjetstructure_R0%.0lf",fR[r]*10);
            //mOutList->Add(new TH2D(hname,"full jet constituents pT vs pTlead ; p_{T}^{part} [GeV/c]; p_{T}^{lead} [GeV/c]", 50, 0, 25, npTleadbins, pTleadmin, pTleadmax));

            hname = Form("hfnparticlesinjet_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH2D(hname,"#particles in full jet vs full jet pT; # of particles; p_{T}^{lead} [GeV/c]", 500*fR[r], 0, 500*fR[r], npTleadbins, pTleadmin, pTleadmax));

						hname = Form("hfpT_R0%.0lf",fR[r]*10);
	    			mOutList->Add(new TH1D(hname, "full jet p_{T}; p_{T} [GeV/c]", nptbins, ptminbin, ptmaxbin)); 

           	hname = Form("hfeta_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH1D(hname, "full jet eta;#eta", netabins, etaminbin, etamaxbin));

	    			hname = Form("hfphi_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH1D(hname, "full jet phi;#phi", nphibins, phiminbin, phimaxbin));

						hname = Form("hNF_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH1D(hname, "jet neutral energy fraction; NEF", 100, 0, 1));
	

            for (int centbin = -1; centbin < 9; centbin++) {
						
						hname = Form("hjetpT_R0%.0lf_centbin%i",fR[r]*10, centbin);
            mOutList->Add(new TH1D(hname, "jet p_{T}; p_{T} [GeV/c]", nptbins, 0, ptmaxbin));

            //hname = Form("hjetpT_R0%.0lf_centbin%i_corr",fR[r]*10, centbin);
            //mOutList->Add(new TH1D(hname, "Corrected jet pT; p_{T} [GeV/c]", nptbins, ptminbin, ptmaxbin));
						//full jet histos
            hname = Form("hfjetpT_R0%.0lf_centbin%i",fR[r]*10, centbin);
            mOutList->Add(new TH1D(hname, "full jet p_{T}; p_{T} [GeV/c]", nptbins, 0, ptmaxbin));



						hname = Form("hjetpTlead_R0%.0lf_centbin%i",fR[r]*10, centbin);
            mOutList->Add(new TH1D(hname, "jet p_{T}^{lead}; p_{T} [GeV/c]", 2*npTleadbins, pTleadmin, pTleadmax));


            hname = Form("hfjetpTlead_R0%.0lf_centbin%i",fR[r]*10, centbin);
            mOutList->Add(new TH1D(hname, "full jet p_{T}^{lead}; p_{T} [GeV/c]", 2*npTleadbins, pTleadmin, pTleadmax));

 						hname = Form("hfjetpTleadNeutral_R0%.0lf_centbin%i",fR[r]*10, centbin);
            mOutList->Add(new TH1D(hname, "full jet p_{T}^{lead, N}; p_{T} [GeV/c]", 2*npTleadbins, pTleadmin, pTleadmax));
  					hname = Form("hfjetpTleadCharged_R0%.0lf_centbin%i",fR[r]*10, centbin);
            mOutList->Add(new TH1D(hname, "full jet p_{T}^{lead, Ch}; p_{T} [GeV/c]", 2*npTleadbins, pTleadmin, pTleadmax));
            
            hname = Form("hfNjets_R0%.0lf_centbin%i",fR[r]*10, centbin);
            mOutList->Add(new TH1D(hname, "number of full jets; N", 100, 0, 100));

            //hname = Form("hfjetpT_R0%.0lf_centbin%i_corr",fR[r]*10, centbin);
            //mOutList->Add(new TH1D(hname, "Corrected full jet pT; p_{T} [GeV/c]", nptbins, ptminbin, ptmaxbin));

            		for(Int_t pTl = 0; pTl < npTlead; pTl++) {
                hname = Form("hpT_pTl%i_R0%.0lf_centbin%i",pTl,fR[r]*10,centbin);
                TString hdesc = Form("jet p_{T} for p_{T}lead>%i ; p_{T}^{corr} [GeV/c]",pTl);
                mOutList->Add(new TH1D(hname, hdesc, nptbins, ptminbin, ptmaxbin));

								hname = Form("hfpT_pTl%i_R0%.0lf_centbin%i",pTl,fR[r]*10,centbin);
                hdesc = Form("full jet p_{T} for p_{T}lead>%i ; p_{T}^{corr} [GeV/c]",pTl);
                mOutList->Add(new TH1D(hname, hdesc, nptbins, ptminbin, ptmaxbin));
            		}
        		}
							
         }
				/*for (int centbin = -1; centbin < 9; centbin++) {
            mOutList->Add(new TH1D(Form("hEcorr_centbin%i", centbin), Form("E from BTow after correction for centbin %i; E [GeV]", centbin), 340, -4.5, 29.5));
						}*/

   		}
		

    if (isMakerMode() == StPicoJetMaker::kWrite) {
        //TODO: Fill trees with required variables
    }

   // mRefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
   // mRefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14.txt");
  return kStOK;
}

// _________________________________________________________
void StPicoHFJetMaker::ClearJets(Option_t *opt="") {
  
  return;
}

// _________________________________________________________
int StPicoHFJetMaker::FinishJets() {
  return kStOK;
}

// _________________________________________________________
int StPicoHFJetMaker::MakeJets() {

    vector<PseudoJet> jetTracks;
		vector<PseudoJet> jetTracks_emb; //tmp for embedding
		vector<PseudoJet> neutraljetTracks; //from bemc clusters only
		vector<PseudoJet> fullTracks;

//very rough centrality estimation - until proper centrality definition available
		//int centrality = -1;

	/*	if (refMult > 10 && refMult <= 21) centrality = 0;
		if (refMult > 21 && refMult <= 40) centrality = 1;
		if (refMult > 40 && refMult <= 71) centrality = 2;
		if (refMult > 71 && refMult <= 115) centrality = 3;
		if (refMult > 115 && refMult <= 176) centrality = 4;
		if (refMult > 176 && refMult <= 257) centrality = 5;
		if (refMult > 257 && refMult <= 364) centrality = 6;
		if (refMult > 364 && refMult <= 430) centrality = 7;
		if (refMult > 430) centrality = 8;
	*/

		fRunNumber = mPicoDst->event()->runId();
		int eventId = mPicoDst->event()->eventId(); //eventID
		int refMult = mPicoDst->event()->refMult();
		//get centrality						
		double vz = mPrimVtx.z();
		mRefmultCorrUtil->setEvent(fRunNumber, refMult, mPicoDst->event()->ZDCx(), vz);
		int centrality = mRefmultCorrUtil->centrality9(); //0 = 0-5 %,..., 8 = 70-80 %
		float weight = 1.0;
		weight = mRefmultCorrUtil->weight();
		static_cast<TH1D*>(mOutList->FindObject("hweight"))->Fill(weight);
		static_cast<TH1D*>(mOutList->FindObject("hcent"))->Fill(centrality, weight);
		
		//static_cast<TH2D*>(mOutList->FindObject("hrunIdcent"))->Fill(fRunNumber,centrality,weight); //not used
		
		//if (centrality > 1) return kStOk; //REMEMBER NOW ONLY CENTRAL

//		cout << "centrality " << centrality << " weight " << weight << " runID " << fRunNumber <<  endl;  

						
						
						//USE reweight from RefmultCorr!! FILL HISTOS WITH WEIGHT

		GetCaloTrackMomentum(mPicoDst,mPrimVtx); //fill array Sump with momenta of tracks which are matched to BEMC
	

		for (int iTow = 0; iTow < 4800; iTow++){ //get btow info
		StPicoBTowHit *towHit = mPicoDst->btowHit(iTow);
		vector<int> ids = {0,0,0,0,0,0,0,0,0}; 
		if (!towHit || towHit->isBad()) {/*cout << "bad tower information" << endl; TowArr.push_back(-1111); TowEta.push_back(-1111); TowPhi.push_back(-1111); Clusters.push_back(ids); */continue;} //if the tower is bad or missing info
		int realtowID = towHit->numericIndex2SoftId(iTow);
		//if (any_of(begin(BadTowerArr), end(BadTowerArr), [=](int n){return n == realtowID;})) {
            cout <<iTow<< " bad tower " << realtowID << endl;
            if (BadTowerMap[realtowID]) {
			//TowArr.push_back(-1111); TowEta.push_back(-1111); TowPhi.push_back(-1111); Clusters.push_back(ids);//if the tower is bad
			continue;} //exclude bad towers

		//double towE = towHit->energy(); //in GeV
		double towE = GetTowerCalibEnergy(iTow+1);
		//cout << towE << endl;
		towE-= fHadronCorr*Sump[iTow];
		//if (Sump[iTow] > 0)cout << "after hadroncorr " << towE << endl;
		if (towE < 0) towE = 0;
				
		StEmcGeom* mEmcGeom;
		mEmcGeom = StEmcGeom::getEmcGeom("bemc");
		float Toweta_tmp = 0, Towphi = 0;
		mEmcGeom->getEtaPhi(realtowID,Toweta_tmp,Towphi);
		float Toweta = vertexCorrectedEta(Toweta_tmp, vz); //max eta 1.05258 max difference: ET = 0.124452 for E = 0.2, if we cut on |Vz| < 30 cm
		static_cast<TH2D*>(mOutList->FindObject("heta_phi_tow"))->Fill(Toweta, Towphi+TMath::Pi(), weight);
		double ET = towE/cosh(Toweta);
		if (ET > 30) {/*cout << towE << endl;*/ TowArr.clear();TowEta.clear();TowPhi.clear();Clusters.clear(); return kStOK;} //discard events with E > 30 GeV towers 
		//no clustering
		double px,py,pz;
		//px = towE*cos(Towphi)/cosh(Toweta);
		//py = towE*sin(Towphi)/cosh(Toweta);
		px = ET*cos(Towphi);
		py = ET*sin(Towphi);
		pz = towE*tanh(Toweta);

		PseudoJet inputTower(px, py, pz, towE);
		if (inputTower.perp() > fETmincut){
		inputTower.set_user_index(0); //default index is -1, 0 means neutral particle
		neutraljetTracks.push_back(inputTower);}
		
		
		//double etacluster = -1111, phicluster = -1111, ecluster = -1111; 

		//mEmcGeom->getEtaPhi(realtowID, eta, phi);
			//if good tower (30 GeV > energy > 0.2 GeV and ADC > 4 default in StPicoCuts) make a 3x3 cluster around this tower (see JetInfo.h)
			/*if (mPicoCuts->isGoodTowHit(towHit) && makeRealCluster(mPicoDst, mEmcGeom, realtowID, &ecluster, &etacluster, &phicluster, &ids)) {


			static_cast<TH2D*>(mOutList->FindObject("heta_phi_tow"))->Fill(etacluster, phicluster+TMath::Pi(), weight);

			TowArr.push_back(ecluster); //fill energy of good towers
			Clusters.push_back(ids); //fill REAL ids of towers in a cluster
			TowEta.push_back(etacluster); TowPhi.push_back(phicluster); //cluster location

			} else {TowArr.push_back(-9999); TowEta.push_back(-9999); TowPhi.push_back(-9999); Clusters.push_back(ids);} //basically if the energy is too low
			*/
		} // get btow info



		//cout << "event " << mPicoDst->event()->eventId() << " with centrality " << centrality << endl;
		//check which towers are in multiple clusters
		/*for (int iTow = 0; iTow < 4800; iTow++){ 
		int towerId = iTow+1;
		vector<unsigned int> clusterids;
		//StPicoBTowHit *towHit = mPicoDst->btowHit(iTow);
		//double towene = towHit->energy();
		double towene = GetTowerCalibEnergy(iTow+1);
		if (TowArr[towerId-1] > 0 && doesExistInMultiple(Clusters, towerId, &clusterids) && clusterids.size() > 1) 
			{
			int n = 0; //n of clusters with the same tower
			int totalene = 0;
			for (unsigned int i = 0; i < clusterids.size(); i++) {if (TowArr[clusterids[i]] > 0) {n++; totalene+=TowArr[clusterids[i]];}}
			if (n > 1) { 
			for (unsigned int i = 0; i < clusterids.size(); i++) 
			TowArr[clusterids[i]]-=(towene*(1-1.0/n)); //split the shared tower's energy equally among the clusters
			//TowArr[clusterids[i]]-=(towene*(1-TowArr[clusterids[i]]/totalene)); //divide the shared tower's energy according to the clusters' energy 
			}
			}
		}*/

		
		//EnergyCl = TowArr;

		//loop over primary tracks
    for (unsigned int i = 0; i < mIdxPicoParticles.size(); i++) {
        StPicoTrack *trk = mPicoDst->track(mIdxPicoParticles[i]);
		
			//using primary tracks
			  double pT = trk->pMom().Perp();
				if(pT != pT) continue; // NaN test. 
				
        float eta = trk->pMom().PseudoRapidity();
        float phi = trk->pMom().Phi();

        float dca = (mPrimVtx - trk->origin()).Mag();
				float charged = trk->charge();


				
        static_cast<TH1D*>(mOutList->FindObject("hpT_tr"))->Fill(pT, weight);
        static_cast<TH2D*>(mOutList->FindObject("heta_phi_tr"))->Fill(phi + TMath::Pi(), eta,  weight);
				static_cast<TH1D*>(mOutList->FindObject("heta_tr"))->Fill(eta, weight);
				static_cast<TH1D*>(mOutList->FindObject("hphi_tr"))->Fill(phi + TMath::Pi(), weight); //to shift by pi
        static_cast<TH2D*>(mOutList->FindObject("hdca_z_tr"))->Fill(dca, vz, weight);
        static_cast<TH2D*>(mOutList->FindObject("hdca_pT"))->Fill(dca, pT, weight);
        static_cast<TH1D*>(mOutList->FindObject("hdca_tr"))->Fill(dca, weight);
				static_cast<TH1D*>(mOutList->FindObject("hcharged_tr"))->Fill(charged, weight);


        //PseudoJet inputParticle(trk->gMom().x(), trk->gMom().y(), trk->gMom().z(), trk->gMom().Mag());
				//primary tracks        
				PseudoJet inputParticle(trk->pMom().x(), trk->pMom().y(), trk->pMom().z(), trk->pMom().Mag());
        jetTracks.push_back(inputParticle);
 
			} 		//end loop over primary tracks


//BEGIN hadronic correction - loop over global tracks  - towers

/*				UInt_t nTracks = mPicoDst->numberOfTracks();
				//cout << nTracks << endl;
				for (unsigned int itrack = 0; itrack < nTracks; itrack++) {

        StPicoTrack *trk = mPicoDst->track(itrack);


				TVector3 gMom = trk->gMom();

				//using global tracks
			  double pT = gMom.Perp();
				if(pT != pT || pT < 0.2) continue;
				
        float eta = gMom.PseudoRapidity();
				if (fabs(eta) > 1) continue;
        float phi = gMom.Phi();

        //float dca = (mPrimVtx - trk->origin()).Mag();
				//if (dca > 10) continue;

				float nHitsFit = trk->nHitsFit();
				float nHitsMax = trk->nHitsMax();
	
				if (nHitsFit < 15 || nHitsFit/nHitsMax < 0.52) continue; //some basic QA cuts

				double Bfield = mPicoDst->event()->bField();

				StPicoPhysicalHelix trkhelix = trk->helix(Bfield);
//begin qa code
				float vtx_x = mPrimVtx.x();
				float vtx_y = mPrimVtx.y();
				float vtx_z = mPrimVtx.z();

				//float dca_xy = float(trkhelix.geometricSignedDistance(vtx_x, vtx_y)); //dca_xy
				TVector3 dcaPoint = trkhelix.at(trkhelix.pathLength(vtx_x, vtx_y));
        float dca_z = dcaPoint.z() - vtx_z; //check

				//if (fabs(dca_xy) > 1) continue;

				//cout << "dca_xy " << dca_xy << " dca_z " << dca_z << endl;
				//static_cast<TH1D*>(mOutList->FindObject("hdca_XY_tr"))->Fill(dca_xy, weight);
				//static_cast<TH1D*>(mOutList->FindObject("hdca_Z_tr"))->Fill(dca_z, weight);
				//static_cast<TH2D*>(mOutList->FindObject("hdca_XYZ_tr"))->Fill(dca_z, dca_xy, weight);
				
				//int lumibin = floor(mPicoDst->event()->ZDCx()/15000);	
				//cout << mPicoDst->event()->ZDCx() << " " << lumibin<< endl;
				//static_cast<TH2D*>(mOutList->FindObject(Form("hdca_XYZ_tr_lumi%i", lumibin)))->Fill(dca_z, dca_xy, weight);
//end qa code

				if (fabs(dca_z) > maxdcazhadroncorr) continue; //PUT IT HERE JUST SO HISTOGRAMS ARE FILLED CORRECTLY, DON'T FORGET
				int TowIndex = -99999;
				TowIndex = trk->bemcTowerIndex();
				cout << TowIndex << endl;


				}// END global track loop


*/


	//BEGIN hadronic correction - loop over global tracks - clusters

				/*UInt_t nTracks = mPicoDst->numberOfTracks();
				//cout << nTracks << endl;
				for (unsigned int itrack = 0; itrack < nTracks; itrack++) {

        StPicoTrack *trk = mPicoDst->track(itrack);


				TVector3 gMom = trk->gMom();

				//using global tracks
			  double pT = gMom.Perp();
				if(pT != pT || pT < 0.2) continue; // NaN test. 
				
        float eta = gMom.PseudoRapidity();
				if (fabs(eta) > 1) continue;
        float phi = gMom.Phi();

        //float dca = (mPrimVtx - trk->origin()).Mag();
				//if (dca > 10) continue;

				float nHitsFit = trk->nHitsFit();
				float nHitsMax = trk->nHitsMax();
	
				if (nHitsFit < 15 || nHitsFit/nHitsMax < 0.52) continue; //some basic QA cuts

				double Bfield = mPicoDst->event()->bField();

				StPicoPhysicalHelix trkhelix = trk->helix(Bfield);
//begin qa code
				float vtx_x = mPrimVtx.x();
				float vtx_y = mPrimVtx.y();
				float vtx_z = mPrimVtx.z();

				//float dca_xy = float(trkhelix.geometricSignedDistance(vtx_x, vtx_y)); //dca_xy
				TVector3 dcaPoint = trkhelix.at(trkhelix.pathLength(vtx_x, vtx_y));
        float dca_z = dcaPoint.z() - vtx_z; //check

				//if (fabs(dca_xy) > 1) continue;

				//cout << "dca_xy " << dca_xy << " dca_z " << dca_z << endl;
				//static_cast<TH1D*>(mOutList->FindObject("hdca_XY_tr"))->Fill(dca_xy, weight);
				//static_cast<TH1D*>(mOutList->FindObject("hdca_Z_tr"))->Fill(dca_z, weight);
				//static_cast<TH2D*>(mOutList->FindObject("hdca_XYZ_tr"))->Fill(dca_z, dca_xy, weight);
				
				//int lumibin = floor(mPicoDst->event()->ZDCx()/15000);	
				//cout << mPicoDst->event()->ZDCx() << " " << lumibin<< endl;
				//static_cast<TH2D*>(mOutList->FindObject(Form("hdca_XYZ_tr_lumi%i", lumibin)))->Fill(dca_z, dca_xy, weight);
//end qa code

				if (fabs(dca_z) > maxdcazhadroncorr) continue; //PUT IT HERE JUST SO HISTOGRAMS ARE FILLED CORRECTLY, DON'T FORGET


				//StEmcPosition *emcposition = new StEmcPosition();
				//if (emcposition) cout << "emcposition" << endl;
				TVector3 *finalposition = new TVector3(0,0,0);
				TVector3 *finalmomentum = new TVector3(0,0,0);

				StPicoPhysicalHelix *helix;
				StPicoPhysicalHelix tmphelix = trk->helix(Bfield);
				helix = &tmphelix;
				//const double bemcradius = 225.405;
				projectTrack(finalposition, finalmomentum, helix, Bfield, 229, 1); //229 cm is the projection radius - so that we project inside the BEMC, not on the inner/outer edge, see JetInfo.h
				double eta_proj = finalposition->PseudoRapidity();
				double phi_proj = finalposition->Phi();
				
				StEmcGeom* mEmcGeom;
				mEmcGeom = StEmcGeom::getEmcGeom("bemc");
		  	Int_t m = 0, e = 0, s = 0, towerId = -1;

  			mEmcGeom->getBin(phi_proj, eta_proj, m, e, s);
  			if (m == 0) continue;
  			if (s < 0) s = 1;
  			mEmcGeom->getId(m, e, s, towerId);
				//int clusterid = findClusterId(towerId); // obsolete 2x2 clustering
				//int clusterid = -1;
				vector<unsigned int> clusterids;
				//cout << "new track with pT " << pT << " in event " << mPicoDst->event()->eventId() << endl;
				//if (doesExist(Clusters, towerId, &clusterid)) 
				if (EnergyCl[towerId-1] > 0 && doesExistInMultiple(Clusters, towerId, &clusterids)) 
				{
				//be carefull, clusterids[i] = position in the Clusters vector, real seed tower ID = clusterids[i]+1
				//cout << "track with pT = " << pT << " matched to tower " << towerId  << " of cluster " << clusterid+1 << " with energy " << EnergyCl[clusterid] << endl;

				int n = 0; //n of clusters with the same tower
				for (unsigned int i = 0; i < clusterids.size(); i++) {if (EnergyCl[clusterids[i]] > 0) n++;} //count good clusters that match
							if (n > 0) { //cout << n << endl;
							for (unsigned int i = 0; i < clusterids.size(); i++) {if (EnergyCl[clusterids[i]] > 0) {Nmatched[clusterids[i]]++; SumpCl[clusterids[i]] += pT*cosh(eta)/n;}}
												 }//if matched to n good clusters, count the track hit and 1/n of its momentum
				} //find clusters with a track match 

				}// END global track loop

*/


				/*for (unsigned int i = 0; i < EnergyCl.size(); i++){
				if (EnergyCl[i] < 0) continue;				
				//
			  //if (Ecorr[i] > 0.5 && (Sump[i] < 0.001)) static_cast<TH2D*>(mOutList->FindObject("hEtrigvsEnext"))->Fill(Ecorr[i], Ecorr[i]+E1[i]+E2[i],  weight);
				
				//hadronic correction 
				//Ecorr[i]-=Sump[i];
				//if (Ecorr[i] < 0) Ecorr[i] = 0;				
	
				//cluster based
				//double enecl = EneCl[i];
				//EneCl[i]-=Sump[i];
				//if (EneCl[i] < 0) EneCl[i] = 0;				

				//real cluster based hadronic correction
				double enecl = EnergyCl[i];
				EnergyCl[i]-=SumpCl[i];
				if (EnergyCl[i] < 0) EnergyCl[i] = 0;	

				//double deltaEcorr = TowArr[i] - Ecorr[i];
				//double deltaEcorr = enecl - EneCl[i];
				double deltaEcorr = enecl - EnergyCl[i];
				//cout << deltaEcorr << endl;
				
				//if (Nmatch[i] > 0) cout << i+1 << " " << Nmatch[i] << endl;

				//int ptbin = ceil(SumpT[i]);
				//int pbin = ceil(Sump[i]);
				//double Rcorr = deltaEcorr/Sump[i];
				//int ptbin = ceil(SumpCl[i]); //!!!!!!! remember this is wrong (asking pT from p)
				int pbin = ceil(SumpCl[i]);
				double Rcorr = deltaEcorr/SumpCl[i];
				//if (Rcorr > 0)    static_cast<TH1D*>(mOutList->FindObject("hRcorr"))->Fill(Rcorr, weight);
				//if (Nmatch[i] > 0) cout << i+1 << " " << Nmatch[i] << endl;
				//if (Rcorr > 0  && ptbin < 11 && Nmatch[i] > 0 && Nmatch[i] < 6)  static_cast<TH1D*>(mOutList->FindObject(Form("hRcorr_nmatch%i_pT%i", Nmatch[i], ptbin)))->Fill(Rcorr, weight);
				//if (Rcorr > 0  && pbin < 11 && Nmatch[i] > 0 && Nmatch[i] < 6)  static_cast<TH1D*>(mOutList->FindObject(Form("hRcorr_nmatch%i_p%i", Nmatch[i], pbin)))->Fill(Rcorr, weight);

				//if (Rcorr > 0 && Nmatch[i] > 0 && Nmatch[i] < 6) static_cast<TH1I*>(mOutList->FindObject("hnmatch"))->Fill(Nmatch[i], weight);

				//if (Rcorr > 0  && ptbin < 11 && Nmatched[i] > 0 && Nmatched[i] < 6)  static_cast<TH1D*>(mOutList->FindObject(Form("hRcorr_nmatch%i_pT%i", Nmatched[i], ptbin)))->Fill(Rcorr, weight);
				//if (Rcorr > 0  && pbin < 11 && Nmatched[i] > 0 && Nmatched[i] < 6)  static_cast<TH1D*>(mOutList->FindObject(Form("hRcorr_nmatch%i_p%i", Nmatched[i], pbin)))->Fill(Rcorr, weight);

				//if (Rcorr > 0 && Nmatched[i] > 0 && Nmatched[i] < 6) static_cast<TH1I*>(mOutList->FindObject("hnmatch"))->Fill(Nmatched[i], weight);
				
				//calculate momentum vector from energy, phi and eta of cluster (assuming m=0)
				double px,py,pz;
				px = EnergyCl[i]*cos(TowPhi[i])/cosh(TowEta[i]);
				py = EnergyCl[i]*sin(TowPhi[i])/cosh(TowEta[i]);
				pz = EnergyCl[i]*tanh(TowEta[i]);

				PseudoJet inputTower(px, py, pz, EnergyCl[i]);
				if (inputTower.perp() > fJetPtMin && inputTower.perp() < 30) 
				{ if (Nmatched[i] == 0) {inputTower.set_user_index(99);} else {inputTower.set_user_index(0);} //default index is -1, distinguish between fully neutral (99), cluster with tracks (0) and charged leading particle
				neutraljetTracks.push_back(inputTower);}
				} //end clusters loop
*/

				fullTracks = neutraljetTracks;
				fullTracks.insert(fullTracks.end(), jetTracks.begin(), jetTracks.end()); //commenting this line will cause only neutral jets, MAX NEUTRAL FRACTION HAS TO BE TURNED OFF

			//make jets

      // background estimation - charged jets
       JetDefinition jet_def_bkgd(kt_algorithm, fRBg);
       AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(fGhostMaxrap, 1, 0.01));
			 if (centrality == 0 || centrality == 1) nJetsRemove = 2;//remove two hardest jets in central collisions, one in others
       Selector selector = SelectorAbsEtaMax(1.0) * (!SelectorNHardest(nJetsRemove)) * SelectorPtMin(0.01);
       JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
       bkgd_estimator.set_particles(jetTracks);

       float rho = bkgd_estimator.rho();		
       float rho_sigma = bkgd_estimator.sigma();
				//cout << "mean area = " << bkgd_estimator.mean_area() << endl;
				//static_cast<TH1D*>(mOutList->FindObject(Form("hrho_R0%.0lf",fR[i]*10)))->Fill(rho_man, weight);
				//static_cast<TH2D*>(mOutList->FindObject(Form("hrho_R0%.0lf_mult",fR[i]*10)))->Fill(refMult, rho_man, weight);
			 static_cast<TH1D*>(mOutList->FindObject("hrho"))->Fill(rho, weight);
			 static_cast<TH2D*>(mOutList->FindObject("hrho_mult"))->Fill(refMult, rho, weight);

       // full jets
       JetDefinition fjet_def_bkgd(kt_algorithm, fRBg);
       AreaDefinition farea_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(fGhostMaxrap, 1, 0.01));
			 if (centrality == 0 || centrality == 1) nJetsRemove = 2;//remove two hardest jets in central collisions, one in others
       Selector fselector = SelectorAbsEtaMax(1.0) * (!SelectorNHardest(nJetsRemove)) * SelectorPtMin(0.01);
       JetMedianBackgroundEstimator fbkgd_estimator(fselector, fjet_def_bkgd, farea_def_bkgd);
       fbkgd_estimator.set_particles(fullTracks);

			 float frho   = fbkgd_estimator.rho();
       float frho_sigma = fbkgd_estimator.sigma();

			 static_cast<TH1D*>(mOutList->FindObject("hfrho"))->Fill(frho, weight);
			 static_cast<TH2D*>(mOutList->FindObject("hfrho_mult"))->Fill(refMult, frho, weight);


    for (unsigned int i = 0; i < fR.size(); i++) {
        float maxRapJet = mPicoCuts->getCutEta() - fR[i];


        if (isMcMode()) {	
		//TODO: Add MC code here
		//embedding of simulated jet into real event
							//vector<PseudoJet> input_vector;
							float pT_emb=0;
							//float eta_emb=0;
							//float phi_emb=0;
	//float maxRapJet=fMaxRap - 0.3; //fiducial jet acceptance

						for (unsigned iemb = 0; iemb < fEmbPt.size(); iemb++)
							{
							//JetReco(input_vector, R, nJetsRemove, weight,iemb);
						
							pT_emb=fEmbPt[iemb];
							//jetTracks_emb = EmbedJet(kEmbPythia,pT_emb,maxRapJet,jetTracks/*,&eta_emb,&phi_emb*/,"u"/*sParton*/,fpythia); //charged jets
							jetTracks_emb = EmbedJet(kEmbPythia,pT_emb,maxRapJet,fullTracks/*,&eta_emb,&phi_emb*/,"u"/*sParton*/,fpythia); //full jets

							//setup fastjet
	            JetDefinition jet_def(antikt_algorithm, fR[i]);
	            // jet area definition
	            //GhostedAreaSpec area_spec(fGhostMaxrap);
	            //AreaDefinition area_def(active_area, are0a_spec);
	            AreaDefinition area_def(active_area_explicit_ghosts, GhostedAreaSpec(fGhostMaxrap, 1, 0.01));
	
	            //run jet reconstruction
	            ClusterSequenceArea clust_seq_hard(jetTracks_emb, jet_def, area_def);
	            vector<PseudoJet> jets_all = sorted_by_pt(clust_seq_hard.inclusive_jets(fJetPtMin));
							//cout << "n all anti-kT jets: " << jets_all.size() << endl;
	            Selector Fiducial_cut_selector = SelectorAbsEtaMax(maxRapJet) * SelectorPtMin(0.01); // Fiducial cut for jets
	            vector<PseudoJet> jets = Fiducial_cut_selector(jets_all);
							//cout << "n accepted anti-kT jets: " << jets.size() << endl;


	            // background estimation
	            JetDefinition jet_def_bkgd(kt_algorithm, fRBg);
	            AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(fGhostMaxrap, 1, 0.01));
							if (centrality == 0 || centrality == 1) nJetsRemove = 2; //remove two hardest jets in central collisions, one in others
	            Selector selector = SelectorAbsEtaMax(1.0) * (!SelectorNHardest(nJetsRemove)) * SelectorPtMin(0.01);
	            JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
	            bkgd_estimator.set_particles(jetTracks_emb);
						

	   					float rho   = bkgd_estimator.rho();
	   					float rho_sigma = bkgd_estimator.sigma();
	
						
		 					for(unsigned int pjet = 0; pjet < jets.size(); pjet++) {
                float phi_jet = jets[pjet].phi();
                float eta_jet = jets[pjet].eta();
                float pT_jet = jets[pjet].perp();
                float area_jet = jets[pjet].area();
                vector<PseudoJet> constituents = sorted_by_pt(jets[pjet].constituents());
							
                float pTlead = constituents[0].perp();
                float pTcorr_jet = pT_jet - area_jet*rho;
							
                //set acceptance
                float etaMinCut = -(maxRapJet);
                float etaMaxCut = (maxRapJet);

                if(eta_jet < etaMinCut || eta_jet > etaMaxCut) continue; // fiducial acceptance

		//*********************
      //FILLING HISTOGRAMS
		//*********************
								bool found=false;
								//is it embedded jet?
								found=FindEmbeddedJet(constituents,pT_emb);
								if(found) {
									static_cast<TH2D*>(mOutList->FindObject(Form("hjetpTembArea_R0%.0lf",fR[i]*10)))->Fill(area_jet, pT_emb, weight);
									if(area_jet < fAcuts[i]) continue;
									double dpT=	pTcorr_jet-pT_emb;
									//if (fabs(fR[i]-0.4) < 0.001) cout << "pTcorr " << pTcorr_jet << " pTemb " << pT_emb << " dpT " << dpT << endl;
      						for(Int_t pTl=0; pTl<npTlead; pTl++)
										{
										if(pTlead >= pTl){
											if (kEmbPythia == 1) {static_cast<TH2D*>(mOutList->FindObject(Form("delta_pt_BG_py_%i_R0%.0lf_centbin%i", pTl, fR[i]*10, centrality)))->Fill(pT_emb, dpT, weight);}
 											else {static_cast<TH2D*>(mOutList->FindObject(Form("delta_pt_BG_sp_%i_R0%.0lf_centbin%i", pTl, fR[i]*10, centrality)))->Fill(pT_emb, dpT, weight);}
											}
										} //pTlead loop
									} //if(found)
								}//jet loop
							}	//emb loop








        } else { //!isMcMode()
						//setup fastjet
            JetDefinition jet_def(antikt_algorithm, fR[i]);
            // jet area definition
            //GhostedAreaSpec area_spec(fGhostMaxrap);
            //AreaDefinition area_def(active_area, are0a_spec);
            AreaDefinition area_def(active_area_explicit_ghosts, GhostedAreaSpec(fGhostMaxrap, 1, 0.01));

            //run jet reconstruction
            ClusterSequenceArea clust_seq_hard(jetTracks, jet_def, area_def);
            vector<PseudoJet> jets_all = sorted_by_pt(clust_seq_hard.inclusive_jets(fJetPtMin));
						//cout << "n all anti-kT jets: " << jets_all.size() << endl;
            Selector Fiducial_cut_selector = SelectorAbsEtaMax(maxRapJet); // Fiducial cut for jets
            vector<PseudoJet> jets = Fiducial_cut_selector(jets_all);
						//cout << "n accepted anti-kT jets: " << jets.size() << endl;


            for(unsigned int pjet = 0; pjet < jets.size(); pjet++) {
                float phi_jet = jets[pjet].phi();
                float eta_jet = jets[pjet].eta();
                float pT_jet = jets[pjet].perp();
                float area_jet = jets[pjet].area();
                vector<PseudoJet> constituents = sorted_by_pt(jets[pjet].constituents());

							
                float pTlead = constituents[0].perp();
                float pTcorr_jet = pT_jet - area_jet*rho;
							
                //set acceptance
                float etaMinCut = -(maxRapJet);
                float etaMaxCut = (maxRapJet);

                if(eta_jet < etaMinCut || eta_jet > etaMaxCut) continue; // fiducial acceptance

                static_cast<TH1D*>(mOutList->FindObject(Form("hjetarea_R0%.0lf",fR[i]*10)))->Fill(area_jet, weight);
                static_cast<TH2D*>(mOutList->FindObject(Form("hjetpTarea_R0%.0lf",fR[i]*10)))->Fill(area_jet, pT_jet, weight);
                static_cast<TH2D*>(mOutList->FindObject(Form("hjetpTcorrArea_R0%.0lf",fR[i]*10)))->Fill(area_jet, pTcorr_jet, weight);

                    if(area_jet < fAcuts[i]) continue;
		
                int nparticles = constituents.size();
								//cout << nparticles << endl;
								for(unsigned int ic = 0; ic < constituents.size(); ++ic) {
									float ceta = constituents[ic].eta();
									float cphi = constituents[ic].phi();
									static_cast<TH1D*>(mOutList->FindObject(Form("hceta_R0%.0lf",fR[i]*10)))->Fill(ceta, weight);
									static_cast<TH1D*>(mOutList->FindObject(Form("hcphi_R0%.0lf",fR[i]*10)))->Fill(cphi, weight);
								}		




                static_cast<TH1D*>(mOutList->FindObject(Form("hjetpT_R0%.0lf_centbin%d",fR[i]*10, centrality)))->Fill(pT_jet, weight);
                static_cast<TH1D*>(mOutList->FindObject(Form("hjetpTlead_R0%.0lf_centbin%d",fR[i]*10, centrality)))->Fill(pTlead, weight);

                    //static_cast<TH1D*>(mOutList->FindObject(Form("hjetpT_R0%.0lf_centbin%d_corr",fR[i]*10, centrality)))->Fill(pTcorr_jet);
                    //static_cast<TH2D*>(mOutList->FindObject(Form("hpT_pTlead_R0%.0lf",fR[i]*10)))->Fill(pTcorr_jet, pTlead);
		    						static_cast<TH1D*>(mOutList->FindObject(Form("hpT_R0%.0lf",fR[i]*10)))->Fill(pT_jet, weight);
                    //static_cast<TH2D*>(mOutList->FindObject(Form("heta_phi_R0%.0lf", fR[i]*10)))->Fill(eta_jet, phi_jet);
		    						static_cast<TH1D*>(mOutList->FindObject(Form("heta_R0%.0lf", fR[i]*10)))->Fill(eta_jet, weight);
		    						static_cast<TH1D*>(mOutList->FindObject(Form("hphi_R0%.0lf", fR[i]*10)))->Fill(phi_jet, weight);
                    static_cast<TH2D*>(mOutList->FindObject(Form("hnparticlesinjet_R0%.0lf",fR[i]*10)))->Fill(nparticles, pTlead); //this includes ghosts and is probably not correct!

                    for(Int_t pTl = 0; pTl < npTlead; pTl++) {
                        if(pTl < pTlead) {
                            static_cast<TH1D*>(mOutList->FindObject(Form("hpT_pTl%i_R0%.0lf_centbin%d",pTl,fR[i]*10, centrality)))->Fill(pTcorr_jet, weight);
                        }
                    }
               

                /*for(int pr = 0; pr < nparticles; pr++) {
                    for (unsigned int j = 0; j < fAcuts.size(); j++) {
                        if(area_jet < fAcuts[j]) continue;

                        float part_pT = constituents[pr].perp();
                        static_cast<TH2D*>(mOutList->FindObject(Form("hjetstructure_R0%.0lf",fR[i]*10)))->Fill(part_pT, pTlead);
                    }
                }*/


            } // for(unsigned int pjet = 0; pjet < jets.size(); pjet++)


					//full jet reconstruction




				    JetDefinition fjet_def(antikt_algorithm, fR[i]);
            // jet area definition
            //GhostedAreaSpec farea_spec(fGhostMaxrap);
            //AreaDefinition area_def(active_area, are0a_spec);
            AreaDefinition farea_def(active_area_explicit_ghosts, GhostedAreaSpec(fGhostMaxrap, 1, 0.01));

            //run full jet reconstruction
					
		
            ClusterSequenceArea fclust_seq_hard(fullTracks, fjet_def, farea_def);

            vector<PseudoJet> fjets_all = sorted_by_pt(fclust_seq_hard.inclusive_jets(fJetPtMin));
            Selector fFiducial_cut_selector = SelectorAbsEtaMax(maxRapJet); // Fiducial cut for jets
            vector<PseudoJet> fjets = fFiducial_cut_selector(fjets_all);


		int naccJets = 0;

						//static_cast<TH2D*>(mOutList->FindObject(Form("hfrho_rho_R0%.0lf",fR[i]*10)))->Fill(rho, frho);
						//cout << "number of jets: " << fjets.size() << endl;
						//double totaljetE = 0;
						for(unsigned int pjet = 0; pjet < fjets.size(); pjet++) {
                float phi_jet = fjets[pjet].phi();
                float eta_jet = fjets[pjet].eta();
                float pT_jet = fjets[pjet].perp();
                float area_jet = fjets[pjet].area();
                vector<PseudoJet> constituents = sorted_by_pt(fjets[pjet].constituents());
                int nparticles = constituents.size();
	//totaljetE += fjets[pjet].E();

               float pTcorr_jet = pT_jet - area_jet*frho;
							
                //set acceptance
                float etaMinCut = -(maxRapJet);
                float etaMaxCut = (maxRapJet);

                if(eta_jet < etaMinCut || eta_jet > etaMaxCut) continue; // fiducial acceptance

                static_cast<TH1D*>(mOutList->FindObject(Form("hfjetarea_R0%.0lf",fR[i]*10)))->Fill(area_jet, weight);
                static_cast<TH2D*>(mOutList->FindObject(Form("hfjetpTarea_R0%.0lf",fR[i]*10)))->Fill(area_jet, pT_jet, weight);
                static_cast<TH2D*>(mOutList->FindObject(Form("hfjetpTcorrArea_R0%.0lf",fR[i]*10)))->Fill(area_jet, pTcorr_jet, weight);

                if(area_jet < fAcuts[i]) continue;


								//is the leading particle charged or neutral?

								float neutralpT = 0;

								for(unsigned int ic = 0; ic < constituents.size(); ++ic) {

									if (constituents[ic].user_index() >= 0) neutralpT+=constituents[ic].pt(); //select clusters
									float ceta = constituents[ic].eta();
									float cphi = constituents[ic].phi();
									static_cast<TH1D*>(mOutList->FindObject(Form("hfceta_R0%.0lf",fR[i]*10)))->Fill(ceta, weight);
									static_cast<TH1D*>(mOutList->FindObject(Form("hfcphi_R0%.0lf",fR[i]*10)))->Fill(cphi, weight);
								}		
								float nfraction = neutralpT/pT_jet;
								//cout << "neutral fraction " << nfraction << endl;
								static_cast<TH1D*>(mOutList->FindObject(Form("hNF_R0%.0lf",fR[i]*10)))->Fill(nfraction, weight);
								if (nfraction > maxneutralfrac) continue; //keep only jets with reasonable cluster fraction of jet pT (default 0.95)	
								
		naccJets++;		
							
                float pTlead = constituents[0].perp();
								
								float leadindex = constituents[0].user_index();
								//cout << leadindex << endl;

 
				
								

								if (leadindex >= 0) { static_cast<TH1D*>(mOutList->FindObject(Form("hfjetpTleadNeutral_R0%.0lf_centbin%d",fR[i]*10, centrality)))->Fill(pTlead, weight);} 
								else { static_cast<TH1D*>(mOutList->FindObject(Form("hfjetpTleadCharged_R0%.0lf_centbin%d",fR[i]*10, centrality)))->Fill(pTlead, weight);}
									
                static_cast<TH1D*>(mOutList->FindObject(Form("hfjetpT_R0%.0lf_centbin%d",fR[i]*10, centrality)))->Fill(pT_jet, weight);
                static_cast<TH1D*>(mOutList->FindObject(Form("hfjetpTlead_R0%.0lf_centbin%d",fR[i]*10, centrality)))->Fill(pTlead, weight);

                    //static_cast<TH1D*>(mOutList->FindObject(Form("hfjetpT_R0%.0lf_centbin%d_corr",fR[i]*10, centrality)))->Fill(pTcorr_jet);
                    //static_cast<TH2D*>(mOutList->FindObject(Form("hpT_pTlead_R0%.0lf",fR[i]*10)))->Fill(pTcorr_jet, pTlead);
		    						static_cast<TH1D*>(mOutList->FindObject(Form("hfpT_R0%.0lf",fR[i]*10)))->Fill(pT_jet, weight);
                    //static_cast<TH2D*>(mOutList->FindObject(Form("heta_phi_R0%.0lf", fR[i]*10)))->Fill(eta_jet, phi_jet);

		    						static_cast<TH1D*>(mOutList->FindObject(Form("hfeta_R0%.0lf", fR[i]*10)))->Fill(eta_jet, weight);
		    						static_cast<TH1D*>(mOutList->FindObject(Form("hfphi_R0%.0lf", fR[i]*10)))->Fill(phi_jet, weight);
                    static_cast<TH2D*>(mOutList->FindObject(Form("hfnparticlesinjet_R0%.0lf",fR[i]*10)))->Fill(nparticles, pTlead);

                    for(Int_t pTl = 0; pTl < npTlead; pTl++) {
                        if(pTl < pTlead) {
                            static_cast<TH1D*>(mOutList->FindObject(Form("hfpT_pTl%i_R0%.0lf_centbin%d",pTl,fR[i]*10, centrality)))->Fill(pTcorr_jet, weight);
                        }
                    }


                /*for(int pr = 0; pr < nparticles; pr++) {
                    for (unsigned int j = 0; j < fAcuts.size(); j++) {
                        if(area_jet < fAcuts[j]) continue;

                        float part_pT = constituents[pr].perp();
                        static_cast<TH2D*>(mOutList->FindObject(Form("hfjetstructure_R0%.0lf",fR[i]*10)))->Fill(part_pT, pTlead);
                    }
                }*/

            } // for(unsigned int pjet = 0; pjet < jets.size(); pjet++)
            //cout << fR[i]<< " " << fjets.size() <<  " " << naccJets << endl;
            static_cast<TH1D*>(mOutList->FindObject(Form("hfNjets_R0%.0lf_centbin%d", fR[i]*10,centrality)))->Fill(naccJets, weight);
           
			//cout << "total jet energy in this event: " << totaljetE << endl;
        } //if (!isMcMode())
    } //for (unsigned int i = 0; i < fR.size(); i++)

	//clear vectors/arrays
	TowArr.clear();
	Ecorr.clear();
	TowEta.clear();
	TowPhi.clear();
	//Nmatch.fill(0);
	//SumpT.fill(0);
	Sump.fill(0);
	//EneCl.fill(0);
	//E1.fill(0);
	//E2.fill(0);

	Clusters.clear();
	Nmatched.fill(0);
	SumpCl.fill(0);
	EnergyCl.clear();
	
  return kStOK;
}

Double_t StPicoHFJetMaker::GetTowerCalibEnergy(Int_t TowerId)
{
  StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(TowerId-1));

  Float_t pedestal, rms;
  Int_t status;
  
  mTables->getPedestal(BTOW, TowerId, 0, pedestal, rms);
  mTables->getStatus(BTOW, TowerId, status);	

  Double_t *TowerCoeff;
  if(fRunNumber <= 15094020) TowerCoeff = CPre;
  else TowerCoeff = CLowMidHigh;
  
  Double_t calibEnergy = TowerCoeff[TowerId-1]*(tower->adc() - pedestal);

  return calibEnergy;
}


//----------------------------------------------------------------------------- 
//Find energy and position of a cluster surrounding known tower - withiout track hit position
//----------------------------------------------------------------------------- 

Bool_t StPicoHFJetMaker::makeRealCluster(StPicoDst *mPicoDst, StEmcGeom *mEmcGeom, int towerId, double *ecluster, double *etacluster, double *phicluster, vector<int> *ids)
{

		ids->resize(9);
   int localTowerId = -1;
   //double energy1 = 0;
   double energyTemp = 0;
   //double dist1 = 1000, dist2 = 1000;
   //double distTemp = 0;
   Float_t etaTemp = 0, phiTemp = 0;
		double cenergy = 0;
		array<double, 9> e;
		array<double, 9> phi;
		array<double, 9> eta;
		e.fill(0);
		phi.fill(0);
		eta.fill(0);
	const float pi = 3.141592;

	//find energy of seed tower
	StPicoBTowHit *emcHit = mPicoDst->btowHit(towerId-1);
		//double seedEne = emcHit->energy();
			double seedEne = GetTowerCalibEnergy(towerId);
	Float_t seedPhi = 0;
	mEmcGeom->getPhi(towerId, seedPhi);		
	//cout << seedPhi << endl;
		//cout << "seed energy " << seedEne << endl;
   //if (mEmcGeom->getId(position.phi(), position.pseudoRapidity(), towerId) == 1) return kTRUE;
  	StEmcPosition *emcpos = new StEmcPosition();

    for (int ieta = -1; ieta < 2; ++ieta) {
      for (int iphi = -1; iphi < 2; ++iphi) {
        localTowerId++;//loops from 0 to 8

        int nextTowerId = emcpos->getNextTowerId(towerId, ieta, iphi);	
				/*bool bad = false;
				for (int i = 0; i<nBadTowers; i++) {if (BadTowerArr[i] == nextTowerId) bad = true;}
				if (bad) {
					//cout << "found bad tower " << nextTowerId <<" , throw cluster away" << endl; 
					for (int i=0; i<9; i++) ids->at(i) = -1; 
					return false;
					}*/
        if (nextTowerId < 1 || nextTowerId > 4800 || BadTowerMap[nextTowerId]) continue;

				ids->at(localTowerId) = nextTowerId;
				//emcHit = mPicoDst->btowHit(nextTowerId-1);
        //if (emcHit == 0) continue;
				//energyTemp = emcHit->energy();
				energyTemp = GetTowerCalibEnergy(nextTowerId);
				//cout << "tower energy " << energyTemp << endl;
				//if (ieta != 0 && iphi !=0 && energyTemp > 0.2 && fabs(energyTemp - seedEne)< 0.001) cout << "same energy in neighbor towers" << nextTowerId << " " << energyTemp << " " << towerId << " " << seedEne << endl;
				if (energyTemp > seedEne) {delete emcpos; emcpos = NULL;return false;}
        if (energyTemp < 0.2) continue; // don't include any noise tower

			
				cenergy+=energyTemp;
				
          mEmcGeom->getEta(nextTowerId, etaTemp);
          mEmcGeom->getPhi(nextTowerId, phiTemp);
          if (fabs(phiTemp - seedPhi) > 0.10){ 
          	if (phiTemp < 0) phiTemp+=2*pi;
          	else phiTemp-=2*pi;
		}
			e[localTowerId] = energyTemp;
			//e[localTowerId] = fabs(energyTemp);
			//cout << towerId << " next tower ID " << nextTowerId << " energy " << e[localTowerId] << " ieta " << ieta << " eta " << etaTemp << " iphi " << iphi << " phi " << phiTemp  << endl;
		//weighing phi and eta of individual towers by their energy
			phi[localTowerId] = phiTemp*e[localTowerId];
			eta[localTowerId] = etaTemp*e[localTowerId];
				//cout << towerId << " energy " << e[localTowerId] << " eta " << eta[localTowerId] << " phi " << phi[localTowerId] << endl;
      } //for (int iphi = -1; iphi < 2; ++iphi)
    }

		*ecluster = cenergy;
		if (cenergy > 0)
		{
		double sum = accumulate(phi.begin(), phi.end(), 0.0);

		*phicluster = sum/cenergy;
		if (fabs(*phicluster)>pi){ //should happen extremely rarely, if ever
			if (*phicluster < -pi) *phicluster += 2*pi;
			else *phicluster -= 2*pi;
			
		}
			//cout << "phi cluster " << *phicluster << endl;
		//cout << accumulate(phi.begin(), phi.end(), sum) << " " << cenergy << endl;
		sum = accumulate(eta.begin(), eta.end(), 0.0);
		*etacluster = sum/cenergy;
		//cout << sum << " " << cenergy << endl;
		
		}

	e.fill(0);
	phi.fill(0);
	eta.fill(0);

	delete emcpos; 
	emcpos = NULL;

	return true;
}



//----------------------------------------------------------------------------- 
//Correct tower eta for Vz position
//----------------------------------------------------------------------------- 
Double_t StPicoHFJetMaker::vertexCorrectedEta(double eta, double vz) {
    double tower_theta = 2.0 * atan(exp(-eta));
    double z = 0.0;
    if (eta != 0.0) z = mBarrelRadius / tan(tower_theta);
    double z_diff = z - vz;
    double theta_corr = atan2(mBarrelRadius, z_diff);
    double eta_corr = -log(tan(theta_corr / 2.0));
    return eta_corr;
}

//----------------------------------------------------------------------------- 
//Fill array with momentum of BEMC-matched tracks
//----------------------------------------------------------------------------- 
Bool_t StPicoHFJetMaker::GetCaloTrackMomentum(StPicoDst *mPicoDst, TVector3 mPrimVtx) {
	//loop over global tracks  - towers

	UInt_t nTracks = mPicoDst->numberOfTracks();
	//cout << nTracks << endl;
	for (unsigned int itrack = 0; itrack < nTracks; itrack++) {
        StPicoTrack *trk = mPicoDst->track(itrack);
	TVector3 gMom = trk->gMom();
	//using global tracks
	double pT = gMom.Perp();
	if(pT != pT || pT < 0.2) continue;
        float eta = gMom.PseudoRapidity();
	if (fabs(eta) > 1) continue;
        float phi = gMom.Phi();

	float nHitsFit = trk->nHitsFit();
	float nHitsMax = trk->nHitsMax();
	if (nHitsFit < 15 || nHitsFit/nHitsMax < 0.52) continue; //some basic QA cuts
	double Bfield = mPicoDst->event()->bField();

	StPicoPhysicalHelix trkhelix = trk->helix(Bfield);
	float vtx_x = mPrimVtx.x();
	float vtx_y = mPrimVtx.y();
	float vtx_z = mPrimVtx.z();


	TVector3 dcaPoint = trkhelix.at(trkhelix.pathLength(vtx_x, vtx_y));
        float dca_z = dcaPoint.z() - vtx_z; //check
        if (fabs(dca_z) > maxdcazhadroncorr) continue; 
	int TowIndex = -99999;
	TowIndex = trk->bemcTowerIndex();
	//cout << TowIndex << endl;
	float p = 0;
	if (TowIndex > 0) {
		p = gMom.Mag();
		Sump[TowIndex-1] += p;
		//cout << p << " " << Sump[TowIndex-1] << endl;
		}
	}// END global track loop
	
	return true;
}
