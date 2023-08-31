#define StMiniMcTree_cxx
#include "StMiniMcTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
/*
  #include "StRoot/StarRoot/TRVector.h"
  #include "StRoot/StarRoot/TRArray.h"
  #include "StRoot/StarRoot/TRMatrix.h"
  #include "StRoot/StarRoot/TRSymMatrix.h"
  #include "StRoot/StarRoot/TRDiagMatrix.h"
  #include "StRoot/StarRoot/TRSymMatrix.h"
  #include "TRVector.h"
  #include "TRMatrix.h"
  #include "TRSymMatrix.h"
*/
#include "TROOT.h"
#include "Riostream.h"
#include <stdio.h>
#include "TSystem.h"
#include "TMath.h"
#include "cuts.h"

void StMiniMcTree::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L StMiniMcTree.C
//      Root > StMiniMcTree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


	int doCentral=atoi(gSystem->Getenv("CENTRAL")); //central or peripheral collisions
	int nfit=atoi(gSystem->Getenv("NFIT")); //number of fit points (used only for output file description)
	TString particle=gSystem->Getenv("PARTICLE"); //pi,K,p
	bool doGlobal=atoi(gSystem->Getenv("GLOBAL")); //global or primary tracks
	int geantID1=-1, geantID2=-1; //geantID of the embedded particle (and antiparticle)

	cout << "doing   " << doCentral << endl;
	if(particle=="P")
	{
		geantID1=14;
		geantID2=15;
	}
	else if(particle=="Pi")
	{
		geantID1=8;
		geantID2=9;
	}
	else if(particle=="K")
	{
		geantID1=11;
		geantID2=12;
	}

   if (fChain == 0) return;
   
   //TFile *myfile = new TFile("results_minimctree_all.root","RECREATE");
   TFile *myfile = new TFile(Form("MC_%s_cent%i_nfit%i.root",particle.Data(),doCentral,nfit),"RECREATE");
   myfile->cd();
   
   int npTbins=200;
   double pTmax=20;
   TH2D* hcent_mult=new TH2D("hcent_mult", "centrality bin vs refMult", 11, -0.5 , 10.5, 150,0, 600); 
   TH2D* hpTRec_pTMc=new TH2D("hpTRec_pTMc", "pT reco vs pT MC", npTbins, 0, pTmax, npTbins, 0,pTmax);
   TH1D* hpTMC_MC=new TH1D("hpTMC_MC", "pTMC of MC tracks", npTbins, 0, pTmax); //for efficiency calculation
   TH1D* hpTMC_Rec=new TH1D("hpTMC_Rec", "pTMC of matched tracks", npTbins, 0, pTmax); //for efficiency calculation
   TH1I* hgeant=new TH1I("hgeant","geant id",1000,0.5,1000.5);
   TH1D* hMP_CommonHit1 = new TH1D("hcommonhit1","distribution of common hits",100,-0.5,99.5); //Jana added
   TH1D* hMP_CommonHit2 = new TH1D("hcommonhit2","distribution of common hits",100,-0.5,99.5); //Jana added      

   //Jana added the sumw2 for the error propagation
   hcent_mult->Sumw2();
   hpTRec_pTMc->Sumw2();
   hpTMC_MC->Sumw2();
   hpTMC_Rec->Sumw2();
   hgeant->Sumw2();

   fChain->SetBranchStatus("*",0);  // disable all branches
   fChain->SetBranchStatus("mVertexZ",1);  // activate branchname
   fChain->SetBranchStatus("mCentrality",1);  // activate branchname
   fChain->SetBranchStatus("mCentralMult",1);  // activate branchname
   fChain->SetBranchStatus("mNMatchedPair",1);  // activate branchname
   fChain->SetBranchStatus("mMatchedPairs.*",1);  // activate branchname
   fChain->SetBranchStatus("mNMcTrack",1);  // activate branchname
   fChain->SetBranchStatus("mMcTracks.*",1);  // activate branchname
   fChain->SetBranchStatus("mNMatGlobPair",1);  // activate branchname
   fChain->SetBranchStatus("mMatGlobPairs.*",1);  // activate branchname

   Long64_t nentries = fChain->GetEntries();
   cout << "nentries = " << nentries << endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
		int enm=fChain->GetEntryNumber(jentry);
		//if(enm%100==0)cout<<"entry #:"<<enm<<endl;
		nb=fChain->GetEntry(jentry);  
		nbytes += nb;

		if(mVertexZ>zcut)continue;
   
		//		cout<<"centrality:"<<mCentrality<<" multiplicity:"<<mCentralMult<<endl;

		if(mCentralMult<=refMultCutMin[doCentral] || mCentralMult>refMultCutMax[doCentral]) continue;
		//		cout<<"cenrality passed " <<endl;  
		hcent_mult->Fill(mCentrality,mCentralMult);
		
		Double_t PMC_trk1=0,PtMC_trk1=0,PR_trk1=0,PtPR_trk1=0,EtaPr_trk1=0,PhiPr_trk1=0,dcaXYGl_trk1=0,dcaGl_trk1=0,dcaXYGlMc_trk1=0,Chi2Pr_trk1=0;
		Int_t TPC_trk1=0,geant_trk1=0,parentGeant_trk1=0,parentKey_trk1,Q_trk1=0, key_trk1=0, Recokey_trk1=0;
		
		//		if(!doGlobal && mNMatchedPair<1) continue; //Jana commented out
		//	if(doGlobal && mNMatGlobPair<1) continue; //Jana commented out
		
		if(jentry%1000 == 0){cout << " jentry = " << jentry << " eventId = " <<  mEventId << "# of mMatGlobPairs = " << mNMatGlobPair << endl;}
		//	cout << "mctracks " << endl;
		for(Int_t trk =0; trk < mNMcTrack; trk++) //loop over MC tracks
		  {
		    //track cuts
		    if(TMath::Abs(mMcTracks_mEtaMc[trk])>EtaMax) continue;
		    //		    if(mMcTracks_mNHitMc[trk] <TpcMcCut) continue; //Jana: not necessary to check this
		    int geant_mctrk  = mMcTracks_mGeantId[trk];
		    if((geantID1>0 && geant_mctrk!=geantID1) && (geantID2>0 && geant_mctrk!=geantID2)) continue; //geant id doesn't match the embedded particle/antiparticle geant id
		    if(mMcTracks_mPtMc[trk]<pTmin || mMcTracks_mPtMc[trk]>pTmax) continue;
		    
		    //fill histograms
		    hpTMC_MC->Fill(mMcTracks_mPtMc[trk]);
		    //  cout << trk << "   " << mMcTracks_mPtMc[trk] << endl;
		    cout <<  "mc  " << mMcTracks_mPtMc[trk] << "   " <<  mMcTracks_mIsValid[trk] << "  " <<  mMcTracks_mIsPrimary[trk] << endl;
		  }
		
		if(!doGlobal) //primary tracks
		  {
		    //cout << "analyzing primary:  "<< mNMatchedPair << endl;
		    //cout << EtaMax << "  " << pTmin << "  " << pTmax << endl;
		    //cout << TpcCut << "  " << RatioHit << "   " << FlagMax << endl;
		    //cout << DcaCutMax << "  " << Chi2Cut << endl;
		    // cout << "primary" << endl;
		    for(Int_t trk =0; trk < mNMatchedPair; trk++)
		      { 
			//track cuts
			//if(TMath::Abs(mMatchedPairs_mEtaMc[trk]) > EtaMax) continue;
			if(TMath::Abs(mMatchedPairs_mEtaPr[trk]) > EtaMax) continue;
			if(mMatchedPairs_mPtPr[trk]<pTmin || mMatchedPairs_mPtPr[trk]>pTmax) continue;
			hMP_CommonHit1->Fill(mMatchedPairs_mNCommonHit[trk]);
			if(mMatchedPairs_mFitPts[trk] < TpcCut) continue;
			if(((float)mMatchedPairs_mFitPts[trk]/(float)mMatchedPairs_mNPossible[trk]) <= RatioHit) continue; //Jana was "<"
			if(mMatchedPairs_mFlag[trk]>FlagMax || mMatchedPairs_mFlag[trk]<=0) continue;

			//			cout <<  mMatchedPairs_mNCommonHit[trk] << endl;
			hMP_CommonHit2->Fill(mMatchedPairs_mNCommonHit[trk]);
			if(mMatchedPairs_mNCommonHit[trk]<10) continue; //Jana added after discussion with Rongrong
			cout << "matched  " << mMatchedPairs_mPtPr[trk] <<  "  " <<  mMatchedPairs_mIsPrimary[trk] << endl;
			//			if(mMatchedPairs_mIsPrimary[trk]!=1) continue; //Jana added after discussio with Rongrong this is doing really weird things!!!
			dcaGl_trk1      = mMatchedPairs_mDcaGl[trk];
			Chi2Pr_trk1 = mMatchedPairs_mChi2Pr[trk];
			if(dcaGl_trk1>=DcaCutMax) continue; // Jana was ">"
			if(Chi2Pr_trk1>Chi2Cut) continue;
			
			geant_trk1       = mMatchedPairs_mGeantId[trk];
			if((geantID1>0 && geant_trk1!=geantID1) && (geantID2>0 && geant_trk1!=geantID2)) continue; //geant id doesn't match the embedded particle/antiparticle geant id
			
			PtPR_trk1        = mMatchedPairs_mPtPr[trk];
			PtMC_trk1        = mMatchedPairs_mPtMc[trk];
			
			//fill histograms
			hpTRec_pTMc->Fill(PtPR_trk1,PtMC_trk1); 
			hpTMC_Rec->Fill(PtMC_trk1);
			hgeant->Fill(geant_trk1);
			//		cout << trk << "  " << PtMC_trk1 << endl;
		      }//track loop
		  }//primary tracks
		else //global tracks
		  {
		    
		//*******************************************
		    
		    for(Int_t trk =0; trk < mNMatGlobPair; trk++)
		      { 
			
			if(mMatGlobPairs_mPtGl[trk]<pTmin || mMatGlobPairs_mPtGl[trk]>pTmax) continue;
			if(mMatGlobPairs_mFitPts[trk] < TpcCut) continue;
			if(((float)mMatGlobPairs_mFitPts[trk]/(float)mMatGlobPairs_mNPossible[trk]) < RatioHit) continue;
			if(TMath::Abs(mMatGlobPairs_mEtaPr[trk]) > EtaMax) continue;
			//if(mMatGlobPairs_mFlag[trk]>FlagMax || mMatGlobPairs_mFlag[trk]<=0) continue;
			
			dcaGl_trk1      = mMatGlobPairs_mDcaGl[trk];
			//Chi2Pr_trk1 = mMatGlobPairs_mChi2Pr[trk];
			//if(dcaGl_trk1>DcaCutMax) continue;
			//if(Chi2Pr_trk1>Chi2Cut) continue;
			
			geant_trk1       = mMatGlobPairs_mGeantId[trk];
			if((geantID1>0 && geant_trk1!=geantID1) && (geantID2>0 && geant_trk1!=geantID2)) continue; //geant id doesn't match the embedded particle/antiparticle geant id
			
			PtPR_trk1        = mMatGlobPairs_mPtGl[trk];
			PtMC_trk1        = mMatGlobPairs_mPtMc[trk];
			
	  		//fill histograms
	  		hpTRec_pTMc->Fill(PtPR_trk1,PtMC_trk1); 
			hpTMC_Rec->Fill(PtMC_trk1);
	  		hgeant->Fill(geant_trk1);
		      }//track loop
		  }//global tracks
   }//event loop
   myfile->Write();
   hpTMC_Rec->Divide(hpTMC_MC);
   hpTMC_Rec->Write("heffi");
   myfile->Close();  
}
