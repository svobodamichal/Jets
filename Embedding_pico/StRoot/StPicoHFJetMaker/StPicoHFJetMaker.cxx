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



//_________________________________match MC tracks to TPC tracks...has to be here, fastjet cannot be included in StPicoHFJetMaker.h file___________________________________________________________
bool MatchTracks(vector<PseudoJet> &McTracks,vector<PseudoJet> &RcTracks){
	bool found = false;
	for (unsigned int i = 0; i < McTracks.size(); i++){
		double mceta = McTracks[i].eta();	
		double mcphi = McTracks[i].phi();	
		double mcpt = McTracks[i].pt();	
		double cut = 0; //delta r cut
		if (mcpt < 1) cut = 0.0475; //cut from deltar distribution minimum - pT dependent
		else if (mcpt < 3) cut = 0.0285;
		else cut = 0.014;
		int nmatch = 0; //number of matched tracks
		int ind = 0; //index of matched track
		vector<int> jvec; //indices of matched track candidates

		for (unsigned int j = 0; j < RcTracks.size(); j++){
			if (RcTracks[j].user_index() > 1000000) continue; //skip already-matched tracks
			double rceta = RcTracks[j].eta();
			double rcphi = RcTracks[j].phi();	
			double deltaeta = mceta-rceta;
			double deltaphi = mcphi-rcphi;
			if (deltaphi > TMath::Pi()) deltaphi = 2.0*TMath::Pi() - deltaphi;
			if (deltaphi < -1.0*TMath::Pi()) deltaphi = 2.0*TMath::Pi() + deltaphi;
			double deltar = sqrt(pow(deltaeta,2)+pow(deltaphi,2));
			if (deltar < cut){ //cout << mcpt << " track match with deltar " << deltar << endl;
					nmatch++;
					jvec.push_back(j);		
					found = true;
			}
		}
		if (nmatch == 0) continue;
		if (nmatch > 1) {//cout << nmatch <<  " matches" << endl; //if more than 1 match, select random match
				srand (time(NULL)); //random seed
				ind = rand() % nmatch; //integer between 0 and nmatch-1 (hopefully)
				//cout << nmatch << " random match number " << ind << endl;
			}

		McTracks[i].set_user_index(1000000+i);
		//cout << McTracks[i].user_index() << endl;
		RcTracks[jvec[ind]].set_user_index(1000000+i);
		//cout << "next mc track" << endl;
	}


	return found;
}

//_________________________________rewritten jet matching using matched tracks...has to be here, fastjet cannot be included in StEffEmbed.h file___________________________________________________________
bool MatchJets(vector<PseudoJet> McJets, vector<PseudoJet> Rcjets, vector<double> McPtLeads, vector<double> Rcleads, vector<pair<PseudoJet,PseudoJet>>* matched, vector<pair<double,double>>* matchedPtLead, 
		vector<pair<double,double>>* matchedNeutralFraction, /*vector<pair<int,int>>* matchedNNeutral, vector<pair<int,int>>* matchedNCharged, vector<pair<int,int>>* matchedNTot,*/ double R){
	bool found = false;
	//double deltaeta = 0, deltaphi = 0;
	//double r = 0; 
	//double pTfracmin = 0.15, pTfracmax = 150; //changed upper limit to a larger number from 1.5
	//double pT, deltapT, deltaR;
	vector<PseudoJet> RcJets = Rcjets; //copy RC jets, so we can remove after match
	vector<double> RcPtLeads = Rcleads; //copy RC leading particles, so we can remove after match
	vector<double> matchpT,matchR,matchpTfrac,matchR2,matchpTfrac2,matchdeltapT; //must be cleared at the end (of first loop)
	vector<pair<PseudoJet,PseudoJet>> matchedTmp, matchedTmp2;	//must be cleared at the end
	vector<pair<double,double>> matchedPtLeadTmp;	//must be cleared at the end
	vector<pair<double,double>> matchedNeutralFractionTmp;//must be cleared at the end
	vector<int> jvec; //indices of matched jet candidates
	vector<int> mcidx; //user indices of MC tracks, must be cleared at the end
	vector<double> matchtrackpT; //pT from matched tracks, must be cleared at the end
	for (unsigned int i = 0; i < McJets.size(); i++) {
		found = false;
		vector<PseudoJet> constituentsMc = sorted_by_pt(McJets[i].constituents());
		//double pTleadMc = constituentsMc[0].perp();
		double nfractionMc = 0;
		//int nneutralMc = 0, nchargedMc = 0;
		double neutralpTMc = 0;
		double pT_jetMc = McJets[i].perp();
		//double pTfrac = 0; //pT fraction det/part
		for(unsigned int ic = 0; ic < constituentsMc.size(); ++ic) {
			int uidx = constituentsMc[ic].user_index();
			if (uidx >= 1000000) {mcidx.push_back(uidx);}//select matched mc tracks
			else if (uidx >= 0) {neutralpTMc+=constituentsMc[ic].perp(); /*nneutralMc++;*/}//select towers
			}	
		nfractionMc = neutralpTMc/pT_jetMc;
       		for (unsigned int j = 0; j < RcJets.size(); j++){
			double pT_jetRc = RcJets[j].perp();
			//if (pT_jetMc > 0) pTfrac = pT_jetRc/pT_jetMc;
			//deltaphi = fabs(McJets[i].phi() - RcJets[j].phi());
			//deltaeta = fabs(McJets[i].eta() - RcJets[j].eta());
			//if (deltaeta < R && deltaphi > (2*TMath::Pi()-R)) deltaphi = 2*TMath::Pi() - deltaphi;
			//r = sqrt(pow(deltaeta,2)+pow(deltaphi,2));
			/*if (r < R && (pTfrac > pTfracmin && pTfrac < pTfracmax)) 
			{
			//matchR.push_back(r);
								matchpTfrac.push_back(pTfrac);
								matchpT.push_back(pT_jetMc);
								matchdeltapT.push_back(pT_jetMc - RcJets[j].perp());
								matchedTmp.push_back(make_pair(McJets[i], RcJets[j]));
								matchedPtLeadTmp.push_back(make_pair(McPtLeads[i], RcPtLeads[j]));
								jvec.push_back(j);
								found = true; 
							*/ //remember, last } of this snippet is missing

			vector<PseudoJet> constituentsRc = sorted_by_pt(RcJets[j].constituents());
			//double pTleadRc = constituentsRc[0].perp();
			double nfractionRc = 0;
			//int nneutralRc = 0, nchargedRc = 0;
			double neutralpTRc = 0;
			double pTmatch = 0;
			for(unsigned int irc = 0; irc < constituentsRc.size(); ++irc) {
				int uidx = constituentsRc[irc].user_index();
				if (uidx >= 1000000) {if ( std::find(mcidx.begin(), mcidx.end(), uidx) != mcidx.end()) pTmatch+=constituentsRc[irc].perp();} //search for RC track user index in the mcidx vector 
				if (uidx >= 0 && uidx < 1000000) {neutralpTRc+=constituentsRc[irc].perp(); /*nneutralRc++;*/}//select towers
				}	
			if (pTmatch > 0) {found = true; // found at least one RC track which matches track at MC level
				matchtrackpT.push_back(pTmatch);
				nfractionRc = neutralpTRc/pT_jetRc;
				matchedNeutralFractionTmp.push_back(make_pair(nfractionMc, nfractionRc));		
				matchedTmp.push_back(make_pair(McJets[i], RcJets[j]));
				matchedPtLeadTmp.push_back(make_pair(McPtLeads[i], RcPtLeads[j]));	
				jvec.push_back(j);		
				} //end of if (pTmatch > 0) - track match candidate found 
			} //end of RC jets loop
			if (!found) continue; //no match
			//int index = std::min_element(matchR.begin(),matchR.end()) - matchR.begin(); //closest jet in R
			//int index = std::max_element(matchpTfrac.begin(),matchpTfrac.end()) - matchpTfrac.begin(); //closest jet in pT
			int index = std::max_element(matchtrackpT.begin(),matchtrackpT.end()) - matchtrackpT.begin(); //largest matched track pT
		
		
			//loop the other way for matched detector level jet
			/*for (uint j = 0; j < McJets.size(); j++) {
				deltaphi = fabs(McJets[j].phi() - matchedTmp[index].second.phi());
				deltaeta = fabs(McJets[j].eta() - matchedTmp[index].second.eta());
				if (deltaeta < R && deltaphi > (2*TMath::Pi()-R)) deltaphi = 2*TMath::Pi() - deltaphi;
							r = sqrt(pow(deltaeta,2)+pow(deltaphi,2));

							if (McJets[j].perp() > 0) pTfrac = matchedTmp[index].second.perp()/McJets[j].perp();
							if (r < R && (pTfrac > pTfracmin && pTfrac < pTfracmax))
							{
								matchR2.push_back(r);
								matchpTfrac2.push_back(pTfrac);
								matchedTmp2.push_back(make_pair(McJets[j], matchedTmp[index].second));
							}

				}*/
			//int index2 = std::min_element(matchR2.begin(),matchR2.end()) - matchR2.begin();

			//int index2 = std::max_element(matchpTfrac2.begin(),matchpTfrac2.end()) - matchpTfrac2.begin(); //closest jet in pT
			//cout << matchedTmp2[index2].first.perp() << " " << matchedTmp[index].first.perp() << endl;
			//if (fabs(matchedTmp2[index2].first.perp() - matchedTmp[index].first.perp()) < 0.00001) {
				//cout << "found 2way match between MC jet " << matchedTmp[index].first.perp() << " and RC jets " << matchedTmp[index].second.perp() << " which should correspond to " << RcJets[jvec[index]].perp() << endl;
				//apply area and eta cut
			RcJets.erase(RcJets.begin()+jvec[index]); //remove already-matched det-lvl jet
			RcPtLeads.erase(RcPtLeads.begin()+jvec[index]); //and corresponding pTlead
			//cout << "removed jet no. " << jvec[index] << endl;
			double area = matchedTmp[index].second.area();
			double Area_cuts[3] = {0.07,0.2,0.4}; //stupid way to "access" fAcuts
			int ac = R*10-2; //find index R=0.2 -> 0, R=0.3 -> 1, R=0.4 -> 2
			if (area > Area_cuts[ac] && fabs(matchedTmp[index].second.eta()) < 1-R && fabs(matchedTmp[index].first.eta()) < 1-R) { 
				matched->push_back((pair<PseudoJet,PseudoJet>)matchedTmp[index]);
				matchedPtLead->push_back(matchedPtLeadTmp[index]); 
				matchedNeutralFraction->push_back(matchedNeutralFractionTmp[index]);
			}
			//}
	
		jvec.clear();
		mcidx.clear();
		matchtrackpT.clear();
		//matchR.clear();
		//matchR2.clear();
		//matchpTfrac.clear();
		//matchpTfrac2.clear();
		//matchpT.clear();
		//matchdeltapT.clear();
		matchedTmp.clear();
		//matchedTmp2.clear();
		matchedPtLeadTmp.clear();
		matchedNeutralFractionTmp.clear();	
    }

//for (unsigned int i = 0; i < matched->size(); i++) cout << matched->at(i).first.perp() << " " << matched->at(i).second.perp() << endl;

		
    return found;
}

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

	mOutList->Add(new TH1D("hweight", "weight", 135, 0, 3));
	mOutList->Add(new TH1D("hcent", "centrality", 10, -1, 9));

    //General track QA
	mOutList->Add(new TH2D("heta_phi_tr", "track phi vs. eta;#phi [-];#eta [-]", nphibins, phiminbin, phimaxbin, netabins, etaminbin, etamaxbin));
	mOutList->Add(new TH1D("hphi_tr", "track phi;#phi", nphibins, phiminbin, phimaxbin));
	mOutList->Add(new TH1D("heta_tr", "track eta;#eta", netabins, etaminbin, etamaxbin));    
	mOutList->Add(new TH1D("hpT_tr", "track pT; p_{T} (GeV/c)", npttrackbins, pttrackmin, pttrackmax));
//	mOutList->Add(new TH1D("hE_tr", "track E from BEMC; E [GeV]", nEtrackbins, Etrackmin, Etrackmax));
//	mOutList->Add(new TH1D("hpTBEMC_tr", "track pT from BEMC; p_{T} (GeV/c)", nptBEMCtrackbins, ptBEMCtrackmin, ptBEMCtrackmax));
	mOutList->Add(new TH2D("hdca_z_tr", "track DCA vs z-vertex", 90, 0, 3, zbins , zmin, zmax));
    	mOutList->Add(new TH1D("hdca_tr", "track DCA", 90, 0, 3));
    	mOutList->Add(new TH2D("hdca_pT", "track DCA vs. p_{T}", 90, 0, 3, npttrackbins, pttrackmin, pttrackmax));
	mOutList->Add(new TH1D("hcharged_tr", "track charge", 90, 0, 3));	
	
	//MC tracks
	mOutList->Add(new TH2D("hMcEtaPhi", "MC particle eta vs phi; #eta (-); #phi (-)", netabins, etaminbin, etamaxbin,nphibins, phiminbin, phimaxbin));	
		
		

		//mOutList->Add(new TH1D("hdca_XY_tr", "global track DCA_XY; DCA_{xy} [cm]", 200, -20, 20));	
    //mOutList->Add(new TH1D("hdca_Z_tr", "global track DCA_Z; DCA_{z} [cm]", 200, -20, 20));	
    //mOutList->Add(new TH2D("hdca_XYZ_tr", "global track DCA_Z vs DCA_XY; DCA_{z} [cm]; DCA_{xy} [cm]", 200, -20, 20, 200, -20, 20));	
   	//for (int lumi = 0; lumi < 5; lumi++) 
		//mOutList->Add(new TH2D(Form("hdca_XYZ_tr_lumi%i", lumi), Form("global track DCA_Z vs DCA_XY, %i < lumi < %i kHz; DCA_{z} [cm]; DCA_{xy} [cm]", lumi*15, (lumi+1)*15), 200, -20, 20, 200, -20, 20));	


		//mOutList->Add(new TH1D("hEcorr", "E from BTow after correction; E [GeV]", 340, -4.5, 29.5));

    //mOutList->Add(new TH2D("hE_tow", "energy vs tower ID;tower ID;E [GeV]", 4800, 1, 4800, 300, -0.5, 29.5));
    //mOutList->Add(new TH2D("hp_tow", "track momentum vs tower ID;tower ID;p (GeV/c)", 4800, 1, 4800, 300, -0.5, 29.5));
  //mOutList->Add(new TH2D("heta_phi_tow", "tower eta vs phi; #eta [-]; #phi [-]", netabins, etaminbin, etamaxbin, nphibins, phiminbin, phimaxbin));		
	  mOutList->Add(new TH2D("heta_phi_tow", "tower eta vs phi; #eta [-]; #phi [-]", netabins/5, etaminbin, etamaxbin, nphibins, phiminbin, phimaxbin));		

    if (isMcMode()) {cout << "MC mode activated" << endl;
        //TODO: Add required histograms for MC Mode

	for(unsigned int r=0; r<fR.size(); r++)
	{
		TString hname=Form("hjetpTembArea_R0%.0lf",fR[r]*10);
		mOutList->Add(new TH2D(hname,"jet pTemb vs area; A [-]; p_{T} (GeV/c)", 100, 0, 1, nptembbins, ptembminbin, ptembmaxbin));
		TString htitle = "delta pT for BG corrections, using sp probe; p_{T}^{emb} (GeV/c); #delta p_{T} (GeV/c)";
		for (int centbin = -1; centbin < 9; centbin++) {
			for(Int_t pTlcut = 0; pTlcut<npTlead; pTlcut++)	{
				hname=Form("delta_pt_BG_sp_%i_R0%.0lf_centbin%i", pTlcut, fR[r]*10, centbin);
					if (kEmbPythia == 1) {htitle = "delta pT for BG corrections, using pythia probe; p_{T}^{emb} (GeV/c); #delta p_{T} (GeV/c)"; hname=Form("delta_pt_BG_py_%i_R0%.0lf_centbin%i", pTlcut, fR[r]*10, centbin);}
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
	//mOutList->Add(new TH1D(hname,"median energy density; #rho_{ch} [GeV/sr]",50,0,50));		
        hname = "hrho_mult";
        //mOutList->Add(new TH2D(hname,"median energy density vs multiplicity; RefMult [-] ; #rho_{ch} [GeV/sr]",600, 0, 600, 50,0,50));
        hname = "hfrho";
        mOutList->Add(new TH1D(hname,"median energy density; #rho [GeV/sr]",100,0,100));
        hname = "hfrho_mult";
        mOutList->Add(new TH2D(hname,"median energy density vs multiplicity ; RefMult [-] ; #rho [GeV/sr]",600, 0, 600, 100,0,100));

 	if (!isMcMode()) {
        for(unsigned int r = 0; r < fR.size(); r++) {
            //TString hname = Form("hpT_pTlead_R0%.0lf",fR[r]*10);
            //CHARGED JETS
            //mOutList->Add(new TH2D(hname, "jet pTcorr vs pTleading; p_{T} (GeV/c); p_{T}^{lead} (GeV/c)", nptbins, ptminbin, ptmaxbin, npTleadbins, pTleadmin, pTleadmax));
	/*	hname = Form("hpT_R0%.0lf",fR[r]*10);
		mOutList->Add(new TH1D(hname, "pT; p_{T} (GeV/c)", nptbins, ptminbin, ptmaxbin)); 
		hname = Form("heta_R0%.0lf",fR[r]*10);
        	mOutList->Add(new TH1D(hname, "jet eta;#eta", netabins, etaminbin, etamaxbin));
	    	hname = Form("hphi_R0%.0lf",fR[r]*10);
            	mOutList->Add(new TH1D(hname, "jet phi;#phi", nphibins, phiminbin, phimaxbin));
          	hname = Form("hjetarea_R0%.0lf",fR[r]*10);
	        mOutList->Add(new TH1D(hname,"jet area", 100, 0, 1));
	        hname = Form("hjetpTarea_R0%.0lf",fR[r]*10);
        	mOutList->Add(new TH2D(hname,"jet pTmeasured vs area; A [-]; p_{T} (GeV/c)" ,100, 0, 1, nptbins, ptminbin, ptmaxbin));
	        hname = Form("hjetpTcorrArea_R0%.0lf",fR[r]*10);
        	mOutList->Add(new TH2D(hname,"jet pTreco vs area; A [-]; p_{T}^{reco} (GeV/c)",100, 0, 1, nptbins, ptminbin, ptmaxbin));
	        hname = Form("hnparticlesinjet_R0%.0lf",fR[r]*10);
        	mOutList->Add(new TH2D(hname,"#particles in jet vs jet pT; # of particles; p_{T}^{lead} (GeV/c)", 300*fR[r], 0, 300*fR[r], npTleadbins, pTleadmin, pTleadmax));
		hname = Form("hceta_R0%.0lf",fR[r]*10);
		mOutList->Add(new TH1D(hname, "jet constituents eta;#eta", ncetabins, cetaminbin, cetamaxbin));
		hname = Form("hcphi_R0%.0lf",fR[r]*10);
        	mOutList->Add(new TH1D(hname, "jet constituents phi;#phi", ncphibins, cphiminbin, cphimaxbin));
*/
		//FULL JETS
		hname = Form("hfceta_R0%.0lf",fR[r]*10);
	        mOutList->Add(new TH1D(hname, "full jet constituents eta;#eta", ncetabins, cetaminbin, cetamaxbin));
		hname = Form("hfcphi_R0%.0lf",fR[r]*10);
            	mOutList->Add(new TH1D(hname, "full jet constituents phi;#phi", ncphibins, cphiminbin, cphimaxbin));	
		hname = Form("hfjetarea_R0%.0lf",fR[r]*10);
	        mOutList->Add(new TH1D(hname,"full jet area", 100, 0, 1));
	        hname = Form("hfjetpTarea_R0%.0lf",fR[r]*10);
        	mOutList->Add(new TH2D(hname,"full jet pTmeasured vs area; A [-]; p_{T} (GeV/c)", 100, 0, 1, nptbins, ptminbin, ptmaxbin));
	        hname = Form("hfjetpTcorrArea_R0%.0lf",fR[r]*10);
        	mOutList->Add(new TH2D(hname,"full jet pTreco vs area; A [-]; p_{T}^{reco} (GeV/c)", 100, 0, 1, nptbins, ptminbin, ptmaxbin));
	        hname = Form("hfnparticlesinjet_R0%.0lf",fR[r]*10);
        	mOutList->Add(new TH2D(hname,"#particles in full jet vs full jet pT; # of particles; p_{T}^{lead} (GeV/c)", 500*fR[r], 0, 500*fR[r], npTleadbins, pTleadmin, pTleadmax));
		hname = Form("hfpT_R0%.0lf",fR[r]*10);
		mOutList->Add(new TH1D(hname, "full jet p_{T}; p_{T} (GeV/c)", nptbins, ptminbin, ptmaxbin)); 
           	hname = Form("hfeta_R0%.0lf",fR[r]*10);
	        mOutList->Add(new TH1D(hname, "full jet eta;#eta", netabins, etaminbin, etamaxbin));
		hname = Form("hfphi_R0%.0lf",fR[r]*10);
	        mOutList->Add(new TH1D(hname, "full jet phi;#phi", nphibins, phiminbin, phimaxbin));
		hname = Form("hNF_R0%.0lf",fR[r]*10);
	        mOutList->Add(new TH1D(hname, "jet neutral energy fraction; NEF", 100, 0, 1));
	
	        for (int centbin = -1; centbin < 9; centbin++) {
			//hname = Form("hjetpT_R0%.0lf_centbin%i",fR[r]*10, centbin);
            		//mOutList->Add(new TH1D(hname, "jet p_{T}; p_{T} (GeV/c)", nptbins, 0, ptmaxbin));
			//full jet histos
            		hname = Form("hfjetpT_R0%.0lf_centbin%i",fR[r]*10, centbin);
		        mOutList->Add(new TH1D(hname, "full jet p_{T}; p_{T} (GeV/c)", nptbins, 0, ptmaxbin));
			hname = Form("hjetpTlead_R0%.0lf_centbin%i",fR[r]*10, centbin);
		      //  mOutList->Add(new TH1D(hname, "jet p_{T}^{lead}; p_{T} (GeV/c)", 2*npTleadbins, pTleadmin, pTleadmax));
		        hname = Form("hfjetpTlead_R0%.0lf_centbin%i",fR[r]*10, centbin);
            		mOutList->Add(new TH1D(hname, "full jet p_{T}^{lead}; p_{T} (GeV/c)", 2*npTleadbins, pTleadmin, pTleadmax));
			hname = Form("hfjetpTleadNeutral_R0%.0lf_centbin%i",fR[r]*10, centbin);
		        mOutList->Add(new TH1D(hname, "full jet p_{T}^{lead, N}; p_{T} (GeV/c)", 2*npTleadbins, pTleadmin, pTleadmax));
 			hname = Form("hfjetpTleadCharged_R0%.0lf_centbin%i",fR[r]*10, centbin);
            		mOutList->Add(new TH1D(hname, "full jet p_{T}^{lead, Ch}; p_{T} (GeV/c)", 2*npTleadbins, pTleadmin, pTleadmax));
                        hname = Form("hfNjets_R0%.0lf_centbin%i",fR[r]*10, centbin);
            		mOutList->Add(new TH1D(hname, "number of full jets; N", 100, 0, 100));
            		hname = Form("hetaphi_MCmatched_R0%.0lf_centbin%i",fR[r]*10, centbin);
            		mOutList->Add(new TH2D(hname, ";eta^{MC}_{jet}; phi^{MC}_{jet}",netabins, etaminbin, etamaxbin,nphibins,phiminbin,phimaxbin));
                       	hname = Form("hetaphi_RCmatched_R0%.0lf_centbin%i",fR[r]*10, centbin);
            		mOutList->Add(new TH2D(hname, ";eta_{jet}; phi_{jet}",netabins, etaminbin, etamaxbin,nphibins,phiminbin,phimaxbin)); 		
            		hname = Form("hpTleads_R0%.0lf_centbin%i",fR[r]*10, centbin);
            		mOutList->Add(new TH2D(hname, "; p^{det}_{T,lead} (GeV/c); p^{true}_{T,lead} (GeV/c)",2*npTleadbins, pTleadmin, pTleadmax,2*npTleadbins, pTleadmin, pTleadmax)); 		
            		
            		for(Int_t pTl = 0; pTl < npTlead; pTl++) {
                		//hname = Form("hpT_pTl%i_R0%.0lf_centbin%i",pTl,fR[r]*10,centbin);
                		TString hdesc = Form("jet p_{T} for p_{T}lead>%i ; p_{T}^{reco} (GeV/c)",pTl);
                		//mOutList->Add(new TH1D(hname, hdesc, nptbins, ptminbin, ptmaxbin));
                		
				hname = Form("hfpT_pTl%i_R0%.0lf_centbin%i",pTl,fR[r]*10,centbin);
               			hdesc = Form("full jet p_{T} for p_{T}lead>%i ; p_{T}^{reco} (GeV/c)",pTl);
            			mOutList->Add(new TH1D(hname, hdesc, nptbins, ptminbin, ptmaxbin));
                
              			hname = Form("hMCpT_pTl%i_R0%.0lf_centbin%i",pTl,fR[r]*10,centbin);
              			hdesc = Form("MC jet p_{T} for p_{T}lead>%i ; p^{MC}_{T} (GeV/c)",pTl);
              			mOutList->Add(new TH1D(hname, hdesc, nptbins, ptminbin, ptmaxbin));
              			
              			hname = Form("hResponseMatrix_pTl%i_R0%.0lf_centbin%i",pTl,fR[r]*10,centbin);
              			hdesc = "; p^{det}_{T}; p^{true} (GeV/c)";
              			mOutList->Add(new TH2D(hname, hdesc, 420, -20, 100, 420, -20, 100));
              						 	
              			hname = Form("hMCmatchedpT_pTl%i_R0%.0lf_centbin%d",pTl,fR[r]*10, centbin);
              			hdesc = "MC matched jets (p_{T,lead} on MC only); p^{true}_{T}";
              			mOutList->Add(new TH1D(hname, hdesc, 600, 0, ptmaxbin));
              			
              			hname = Form("hRCmatchedpT_pTl%i_R0%.0lf_centbin%d",pTl,fR[r]*10, centbin);
              			hdesc = "RC matched jets (p_{T,lead} on MC only); p^{true}_{T}";
              			mOutList->Add(new TH1D(hname, hdesc, 600, 0, ptmaxbin));
            		}
               	}
							
        } // if (!isMcMode()) {
	} // for(unsigned int r = 0; r < fR.size(); r++) 
		

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
	vector<PseudoJet> neutraljetTracks; //from bemc towers only
	vector<PseudoJet> fullTracks;
	vector<PseudoJet> MCjetTracks;
	vector<PseudoJet> MCjetTowers;	

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
	double vr = TMath::Sqrt(mPrimVtx.x()*mPrimVtx.x() + mPrimVtx.y()*mPrimVtx.y());
	if (vr > 2) cout << "Vr: " << vr << endl; 
	mRefmultCorrUtil->setEvent(fRunNumber, refMult, mPicoDst->event()->ZDCx(), vz);
	int centrality = mRefmultCorrUtil->centrality9(); //0 = 0-5 %,..., 8 = 70-80 %
	float Weight = 1.0;
	Weight = mRefmultCorrUtil->weight();
	float weight = Weight*fWeight; //centrality weight * cross section weight
	static_cast<TH1D*>(mOutList->FindObject("hweight"))->Fill(weight);
	static_cast<TH1D*>(mOutList->FindObject("hcent"))->Fill(centrality, weight);
		
	//if (centrality > 1) return kStOk; //REMEMBER NOW ONLY CENTRAL
	
	//THIS FUNCTION WILL NOT WORK ON EMBEDDING, UNLESS THE EMBEDDING IS INTO HT EVENTS				
	//if (!FindTriggerTowers(2)) return kStOk; //2 = HT2, don't continue if there is no HT2-trigger tower with sufficient energy
	
	
	//MC tracks
	int noMCtracks = mPicoDst->numberOfMcTracks();
	for (int i = 0; i < noMCtracks; i++){
		StPicoMcTrack *mctrk = (StPicoMcTrack*)mPicoDst->mcTrack(i);
		if(mctrk->idVtxStart() > 1) continue; //only primary tracks
		int geantId = mctrk->geantId();
		double mcpt = mctrk->pt();
		double mceta = mctrk->eta();
		if ((geantId > 2 && geantId < 7) || fabs(mceta) > 1.0 || mcpt < 0.2) continue;
		double mcp = mctrk->ptot(); 
		TVector3 mcmom = mctrk->p();
   		double mcphi =mcmom.Phi();
  		if(mcphi<0.0) mcphi += 2.0*TMath::Pi(); 
    		if(mcphi>2.0*TMath::Pi()) mcphi -= 2.0*TMath::Pi();

		static_cast<TH2D*>(mOutList->FindObject("hMcEtaPhi"))->Fill(mceta, mcphi, weight); 

		double mcpx,mcpy,mcpz;
		mcpx = mcmom.x();
		mcpy = mcmom.y();
		mcpz = mcmom.z();

		double mcE = mctrk->energy();
		
		PseudoJet inputMcParticle(mcpx, mcpy, mcpz, mcE);
		//PseudoJet inputNeutralMcParticle(mcpx, mcpy, mcpz, mcp); //assume m = 0
		//cout << inputNeutralMcParticle.perp() << endl; 
		if (!mctrk->charge()) {inputMcParticle.set_user_index(geantId);MCjetTowers.push_back(inputMcParticle);} 
		else {MCjetTracks.push_back(inputMcParticle);}  	
	}
	


	//RC part
	
	
	GetCaloTrackMomentum(mPicoDst,mPrimVtx); //fill array Sump with momenta of tracks which are matched to BEMC

	int TOWE = 0;
	for (int iTow = 0; iTow < 4800; iTow++){ //get btow info
		StPicoBTowHit *towHit = mPicoDst->btowHit(iTow);
		if (!towHit || towHit->isBad()) continue; //if the tower is marked as bad or missing info
		int realtowID = towHit->numericIndex2SoftId(iTow);
		if (BadTowerMap[realtowID]) continue; //exclude bad towers (map in JetInfo.h)

		double towE = GetTowerCalibEnergy(iTow+1); //get tower energy
		TOWE=towE; //just keep track of the original energy for trigger approximation
		towE-= fHadronCorr*Sump[iTow]; //subtract hadronic energy deposition
		if (towE < 0) towE = 0;
						
		StEmcGeom* mEmcGeom;
		mEmcGeom = StEmcGeom::getEmcGeom("bemc");
		float Toweta_tmp = 0, Towphi = 0;
		mEmcGeom->getEtaPhi(realtowID,Toweta_tmp,Towphi);
		float Toweta = vertexCorrectedEta(Toweta_tmp, vz); //max eta 1.05258 max difference: ET = 0.124452 for E = 0.2, if we cut on |Vz| < 30 cm
		static_cast<TH2D*>(mOutList->FindObject("heta_phi_tow"))->Fill(Toweta, Towphi+TMath::Pi(), weight);
		double ET = towE/cosh(Toweta);
		if (ET > 30) {/*cout << towE << endl;*/ return kStOK;} //discard events with E > 30 GeV towers 
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
		//THIS LINE WILL NOT WORK
		//if (find(Triggers.begin(), Triggers.end(), realtowID)!=Triggers.end()) inputTower.set_user_index(2); //mark trigger towers with user_index 2
		if (TOWE > fTrgthresh) inputTower.set_user_index(2); //mark trigger towers with user_index 2	
		neutraljetTracks.push_back(inputTower);}
	} //end get btow info

	//loop over primary tracks
	for (unsigned int i = 0; i < mIdxPicoParticles.size(); i++) {
        	StPicoTrack *trk = mPicoDst->track(mIdxPicoParticles[i]);
		double pT = trk->pMom().Perp(); //using primary tracks
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

	//match tracks
	MatchTracks(MCjetTracks,jetTracks);

	fullTracks = neutraljetTracks;
	fullTracks.insert(fullTracks.end(), jetTracks.begin(), jetTracks.end()); //commenting this line will cause only neutral jets, MAX NEUTRAL FRACTION HAS TO BE TURNED OFF

	//make charged jets
	//background estimation - charged jets
	/*JetDefinition jet_def_bkgd(kt_algorithm, fRBg);
	AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(fGhostMaxrap, 1, 0.01));
	if (centrality == 0 || centrality == 1) nJetsRemove = 2;//remove two hardest jets in central collisions, one in others
	Selector selector = SelectorAbsEtaMax(1.0) * (!SelectorNHardest(nJetsRemove)) * SelectorPtMin(0.01);
	JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
	bkgd_estimator.set_particles(jetTracks);

	float rho = bkgd_estimator.rho();		
	float rho_sigma = bkgd_estimator.sigma();
	static_cast<TH1D*>(mOutList->FindObject("hrho"))->Fill(rho, weight);
	static_cast<TH2D*>(mOutList->FindObject("hrho_mult"))->Fill(refMult, rho, weight);
	*/
	// full jets
	JetDefinition fjet_def_bkgd(kt_algorithm, fRBg);
	AreaDefinition farea_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(fGhostMaxrap, 1, 0.01));
	if (centrality == 0 || centrality == 1) nJetsRemove = 2;//remove two hardest jets in central collisions, one in others
	Selector fselector = SelectorAbsEtaMax(1.0) * (!SelectorNHardest(nJetsRemove)) * SelectorPtMin(0.01);
	JetMedianBackgroundEstimator fbkgd_estimator(fselector, fjet_def_bkgd, farea_def_bkgd);
	fbkgd_estimator.set_particles(fullTracks);

	float frho   = fbkgd_estimator.rho();
	//float frho_sigma = fbkgd_estimator.sigma();
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
	        float maxRapJet = 1 - fR[i];
		//MC jets first
		
		MCjetTracks.insert(MCjetTracks.end(),MCjetTowers.begin(),MCjetTowers.end());
		//setup fastjet
		JetDefinition Mcjet_def(antikt_algorithm, fR[i]);
        	AreaDefinition Mcarea_def(active_area_explicit_ghosts, GhostedAreaSpec(fGhostMaxrap, 1, 0.01));
            //run jet reconstruction
		ClusterSequenceArea Mcclust_seq_hard(MCjetTracks, Mcjet_def, Mcarea_def);
		vector<PseudoJet> Mcjets_all = sorted_by_pt(Mcclust_seq_hard.inclusive_jets(fJetPtMin));
          	Selector McFiducial_cut_selector = SelectorAbsEtaMax(maxRapJet)* SelectorPtMin(0.01)* SelectorPtMax(3*fpThatmin); //throw out jets with pT larger than 3*pThat to eliminate high-weight fluctuations; // Fiducial cut for jets
            	vector<PseudoJet> McJets = McFiducial_cut_selector(Mcjets_all);
		vector<double> McPtLeads;

		for(unsigned int pjet = 0; pjet < McJets.size(); pjet++) {
                	float phi_jet = McJets[pjet].phi();
                	float eta_jet = McJets[pjet].eta();
               		float pT_jet = McJets[pjet].perp();
                	float area_jet = McJets[pjet].area();
               	 	vector<PseudoJet> constituents = sorted_by_pt(McJets[pjet].constituents());							
                	float pTlead = constituents[0].perp();
			McPtLeads.push_back(pTlead);
                	for(Int_t pTl = 0; pTl < npTlead; pTl++) {
                        	if(pTl < pTlead) static_cast<TH1D*>(mOutList->FindObject(Form("hMCpT_pTl%i_R0%.0lf_centbin%d",pTl,fR[i]*10, centrality)))->Fill(pT_jet, weight);
                    	}
  
            	} // for(unsigned int pjet = 0; pjet < McJets.size(); pjet++)
        	
        	
        	//RC jets
		//setup fastjet
		//CHARGED JETS
	/*	JetDefinition jet_def(antikt_algorithm, fR[i]);
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
             
           	 } // for(unsigned int pjet = 0; pjet < jets.size(); pjet++)

		*/
		//full jet reconstruction
		JetDefinition fjet_def(antikt_algorithm, fR[i]);
		AreaDefinition farea_def(active_area_explicit_ghosts, GhostedAreaSpec(fGhostMaxrap, 1, 0.01));
        	//run full jet reconstruction
		ClusterSequenceArea fclust_seq_hard(fullTracks, fjet_def, farea_def);
            	vector<PseudoJet> fjets_all = sorted_by_pt(fclust_seq_hard.inclusive_jets(fJetPtMin));
            	Selector fFiducial_cut_selector = SelectorAbsEtaMax(maxRapJet); // Fiducial cut for jets
            	vector<PseudoJet> fjets = fFiducial_cut_selector(fjets_all);
		int naccJets = 0;
		vector<double> RcPtLeads;
		vector<PseudoJet> RcJets;
		for(unsigned int pjet = 0; pjet < fjets.size(); pjet++) {
		        vector<PseudoJet> constituents = sorted_by_pt(fjets[pjet].constituents());
			//look only for jets with associated trigger tower
			bool istriggerjet = false;
			for(unsigned int ic = 0; ic < constituents.size(); ++ic) {
				if (constituents[ic].user_index() == 2) {istriggerjet = true; break;}
			}
			if (!istriggerjet) continue; 
                	float phi_jet = fjets[pjet].phi();
                	float eta_jet = fjets[pjet].eta();
                	float pT_jet = fjets[pjet].perp();
                	float area_jet = fjets[pjet].area();
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
			
			float neutralpT = 0;
			//is the leading particle charged or neutral?
			for(unsigned int ic = 0; ic < constituents.size(); ++ic) {
				if (constituents[ic].user_index() >= 0 && constituents[ic].user_index()< 3) neutralpT+=constituents[ic].pt(); //select towers
				float ceta = constituents[ic].eta();
				float cphi = constituents[ic].phi();
				static_cast<TH1D*>(mOutList->FindObject(Form("hfceta_R0%.0lf",fR[i]*10)))->Fill(ceta, weight);
				static_cast<TH1D*>(mOutList->FindObject(Form("hfcphi_R0%.0lf",fR[i]*10)))->Fill(cphi, weight);
			}		
			float nfraction = neutralpT/pT_jet;
			//cout << "neutral fraction " << nfraction << endl;
			static_cast<TH1D*>(mOutList->FindObject(Form("hNF_R0%.0lf",fR[i]*10)))->Fill(nfraction, weight);
			if (nfraction > maxneutralfrac) continue; //keep only jets with reasonable neutral energy fraction of jet pT (default 0.95)	
			naccJets++;		
							
	                float pTlead = constituents[0].perp();
	                RcPtLeads.push_back(pTlead);
	                RcJets.push_back(fjets[pjet]);
			float leadindex = constituents[0].user_index();
			//cout << leadindex << endl;
			if (leadindex >= 0) { static_cast<TH1D*>(mOutList->FindObject(Form("hfjetpTleadNeutral_R0%.0lf_centbin%d",fR[i]*10, centrality)))->Fill(pTlead, weight);} 
			else { static_cast<TH1D*>(mOutList->FindObject(Form("hfjetpTleadCharged_R0%.0lf_centbin%d",fR[i]*10, centrality)))->Fill(pTlead, weight);}
			static_cast<TH1D*>(mOutList->FindObject(Form("hfjetpT_R0%.0lf_centbin%d",fR[i]*10, centrality)))->Fill(pT_jet, weight);
                	static_cast<TH1D*>(mOutList->FindObject(Form("hfjetpTlead_R0%.0lf_centbin%d",fR[i]*10, centrality)))->Fill(pTlead, weight);
			static_cast<TH1D*>(mOutList->FindObject(Form("hfpT_R0%.0lf",fR[i]*10)))->Fill(pT_jet, weight);
			//static_cast<TH2D*>(mOutList->FindObject(Form("heta_phi_R0%.0lf", fR[i]*10)))->Fill(eta_jet, phi_jet);
			static_cast<TH1D*>(mOutList->FindObject(Form("hfeta_R0%.0lf", fR[i]*10)))->Fill(eta_jet, weight);
			static_cast<TH1D*>(mOutList->FindObject(Form("hfphi_R0%.0lf", fR[i]*10)))->Fill(phi_jet, weight);
	                static_cast<TH2D*>(mOutList->FindObject(Form("hfnparticlesinjet_R0%.0lf",fR[i]*10)))->Fill(nparticles, pTlead);

	                for(Int_t pTl = 0; pTl < npTlead; pTl++) {
                     	   if(pTl < pTlead) {static_cast<TH1D*>(mOutList->FindObject(Form("hfpT_pTl%i_R0%.0lf_centbin%d",pTl,fR[i]*10, centrality)))->Fill(pTcorr_jet, weight);}
                    	}


		} // for(unsigned int pjet = 0; pjet < jets.size(); pjet++)
            	static_cast<TH1D*>(mOutList->FindObject(Form("hfNjets_R0%.0lf_centbin%d", fR[i]*10,centrality)))->Fill(naccJets, weight);
           	//cout << "total jet energy in this event: " << totaljetE << endl;
           	
           	//JET MATCHING
           	//cout << "JET MATCHING" << endl;
           	//cout << "Number of MC jets: " << McJets.size() << " number of RC jets: " << RcJets.size() << endl;
		if (RcJets.size() == 0) continue;
		vector<pair<PseudoJet, PseudoJet>> Matched;
		vector<pair<double, double>> MatchedpTleads;
		vector<pair<double, double>> MatchedNeutralFraction;
		//vector<pair<int, int>> MatchedNNeutral, MatchedNCharged, MatchedNTot;
		MatchJets(McJets, RcJets, McPtLeads, RcPtLeads, &Matched, &MatchedpTleads, &MatchedNeutralFraction, /*&MatchedNNeutral, &MatchedNCharged, &MatchedNTot, */fR[i]);
		//cout << deltaR << " " << deltapT << " " << pTtrue << endl;

		for (unsigned int j = 0; j < Matched.size(); j++) {
			double pT_det = Matched[j].second.perp();
			double pT_true = Matched[j].first.perp();
			double pT_corr_det = pT_det - Matched[j].second.area()*frho;
			double MCmatchedeta = Matched[j].first.eta();
			double MCmatchedphi = Matched[j].first.phi();
			double RCmatchedeta = Matched[j].second.eta();
			double RCmatchedphi = Matched[j].second.phi();				

			static_cast<TH2D*>(mOutList->FindObject(Form("hetaphi_MCmatched_R0%.0lf_centbin%i",fR[i]*10, centrality)))->Fill(MCmatchedeta,MCmatchedphi, weight);
			static_cast<TH2D*>(mOutList->FindObject(Form("hetaphi_RCmatched_R0%.0lf_centbin%i",fR[i]*10, centrality)))->Fill(RCmatchedeta,RCmatchedphi, weight);	

			//double pTlead = MatchedpTleads[j].second; //pTlead cut on detector level only
			double MatchedMCpTlead = MatchedpTleads[j].first;
			double MatchedRCpTlead = MatchedpTleads[j].second;			
			double pTlead = min(MatchedMCpTlead,MatchedRCpTlead); //pTlead cut on both levels 
			double matchedNF = MatchedNeutralFraction[j].first;
			if (matchedNF < 0.01) continue; //throw out track-only MC jets
			static_cast<TH2D*>(mOutList->FindObject(Form("hpTleads_R0%.0lf_centbin%i",fR[i]*10, centrality)))->Fill(MatchedRCpTlead,MatchedMCpTlead,weight);
			for(Int_t pTl = 0; pTl < npTlead; pTl++) {
				if(pTl < pTlead) {
				static_cast<TH2D*>(mOutList->FindObject(Form("hResponseMatrix_pTl%i_R0%.0lf_centbin%i",pTl,fR[i]*10,centrality)))->Fill(pT_corr_det, pT_true, weight);
				} 
			 	if (pTl < MatchedpTleads[j].first){
			 	static_cast<TH1D*>(mOutList->FindObject(Form("hMCmatchedpT_pTl%i_R0%.0lf_centbin%d",pTl,fR[i]*10, centrality)))->Fill(pT_true, weight);
			 	}
			}
		}
		Matched.clear();
		MatchedpTleads.clear();
		MatchedNeutralFraction.clear();
       		} //if (!isMcMode())
	} //for (unsigned int i = 0; i < fR.size(); i++)

	//clear vectors/arrays
	Sump.fill(0);
	Triggers.clear();
	
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

//----------------------------------------------------------------------------- 
//Fill array with ID of towers which are marked as HTlevel and are above the desired trigger threshold 
// IN DIJET EMBEDDING INTO MB, THERE IS NO HT TRIGGERS
//----------------------------------------------------------------------------- 
 Int_t StPicoHFJetMaker::FindTriggerTowers(Int_t level = 2) {
 
 	if (level < 1 || level > 3) {cout << "Wrong trigger level, this function cannot be used" << endl; return -1;}
 	
 	UInt_t ntrg = mPicoDst->numberOfEmcTriggers();
	for (int i = 0; i < ntrg; i++){ //get btow info
		StPicoEmcTrigger *trg = mPicoDst->emcTrigger(i);
		if (!trg) {/*cout << "no trigger info" << endl;*/ continue;} 
		if ((level == 2 && !trg->isHT2()) || (level == 1 && !trg->isHT1()) || (level == 3 && !trg->isHT3())) {/*cout << "not HT2 trigger" << endl;*/ continue;}
		int towid = trg->id();
		if (BadTowerMap[towid]) continue; 
		float energy = GetTowerCalibEnergy(towid);
		if (energy > fTrgthresh) Triggers.push_back(towid);
	} 	
	
	return Triggers.size();
 }