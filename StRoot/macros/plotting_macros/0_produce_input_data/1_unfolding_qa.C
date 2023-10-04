#include "../utils/common.h"
#include "../utils/statistic_tests.C"
#include "../utils/utils.C"

void unfolding_qa(C_systems system=cent, TString uname="GPC2",TString inputDir="./unfolded_data",/*int binning=0,*/bool doToymodel=0, bool verbose=0)
{
   bool ppHT=1; //0: use only MB data for pp | 1: use combined MB+HT data for pp
	const int nunf=2; //number of unfoldings (usually 2: SVD+Bayes)
	//TString C_unfoldings_name[]={"Bayes","SVD"}; //moved to utils/common.h
	const int nprior=16; //16
	TString prior_type[]={"","flat","pythia","powlaw3","powlaw45","powlaw5","powlaw55","tsalis_1","tsalis_2","tsalis_3","tsalis_4","tsalis_5","tsalis_6","tsalis_7","tsalis_8","tsalis_9"};
	int prior_start=2;
	TString nbins="VAR";
	//const int C_nR=3; //moved to utils/common.h
	//const float C_Rs[]={0.2,0.3,0.4,0.5}; //moved to utils/common.h
	int nbinnings=3;
	if(doToymodel)nbinnings=1;
	const int binnings[]={1,4,0,2,3}; 
	int ptl_start=5;
   int ptl_end=7; 
   if(system==peri)
   {
      ptl_start=4;
      ptl_end=7;
   }
	else if(system==pp) //pp data
	{
      ptl_start=3;
      ptl_end=6;
	}
	//const float C_pTls[]={0,1,3,4,5,6,7}; moved to utils/common.h
	int maxiter[]={8,8}; //number of unfolding iterations
	
	float RAAtoy=(system==cent) ? 0.5 : 0.7; //toymodel RAA of hard jet spectrum
	int nsuf=11; //11
	const TString suff[]={"normal","pythiaUG","RRho02","RRho04","nrem-1", "pp","v2","m5","p5","u","g"}; //,"nfit13","nfit14","nfit18","nfit20"};
  	
	float mincutoff[]={10,10,10}; //low-pT cutoff for backfolding test
	float mincutof_curv[]={0,7,7}; //mincutoff for curvature test for central,  peripheral and pp collisons
	float maxcutoff[]={30,22,25}; //high-pT cutoff
	
	TString dtype="data";
	if(doToymodel)dtype="toy";
	TString toyraa="";
	if(doToymodel) toyraa=Form("_RAA%.1lf",RAAtoy);
	
	
		Color_t colorList[]={kBlack,kBlue,kRed,kGreen,kMagenta,kOrange};
	  /*TH1 *frame = new TH1I("frame", "", 1000, 0, +100);
	  frame->GetXaxis()->SetRangeUser(0, 30);
  frame->GetYaxis()->SetRangeUser(1E-6, 1E6);
   TCanvas *cspectra = new TCanvas("cspectra","cspectra",10,10,900,600);
  cspectra->cd();
   cspectra->SetLogy();
	frame->DrawCopy("");*/
	
	for(int bn=0;bn<nbinnings;bn++)
	{
		int binning=binnings[bn];
		cout<<"binning: "<<binnings[bn]<<"("<<bn+1<<"/"<<nbinnings<<")"<<endl;
	for(int sf=0; sf<nsuf;sf++)
	{
		//if(doToymodel && sf>0)continue;
	cout<<"type: "<<suff[sf].Data()<<"("<<sf+1<<"/"<<nsuf<<")"<<endl;
		
	TString resDir="inclusive";
	if(doToymodel) resDir="toymodel";
	if(system==peri) resDir+="_peripheral";
	else if(system==cent) resDir+="_central";
	else if(system==pp) 
    {
        if(ppHT) resDir+="HT_pp";
        else resDir+="_pp";
    }

		if(suff[sf]=="normal" || suff[sf]=="pp" || suff[sf]=="v2" || suff[sf]=="m5" || suff[sf]=="p5" || suff[sf]=="u" || suff[sf]=="g")  //these unfoldings are put in the inclusive_X_main directory
			resDir+="_main";
		else
			resDir+=Form("_%s",suff[sf].Data());
	for(int r=0; r<C_nR; r++){
		//if(r==0)continue;
		//if(r==2)continue;
		cout<<"R="<<C_Rs[r]<<endl;
      for(int ptl=ptl_start; ptl<=ptl_end; ptl++){
			cout<<"   pTlead>"<<C_pTls[ptl]<<endl;
            //if(C_pTls[ptl]==2) continue;
			for(int unf=0; unf<nunf; unf++){
				cout<<"      "<<C_unfoldings_name[unf].Data()<<endl;
				for(int prior=prior_start; prior<nprior; prior++){
               if (prior==3)continue;
					//cout<<"         "<<prior_type[prior].Data()<<endl;

	TString cortype=(system==pp) ? "dete" : "BGD";
	TString wrkdir = Form("%s/%s/Unfolded_R%.1lf_%s_%sbins_bining%i_%s%s_%s_%s",inputDir.Data(),resDir.Data(),C_Rs[r],C_unfoldings_name[unf].Data(),nbins.Data(),binning,cortype.Data(),toyraa.Data(),uname.Data(),suff[sf].Data());
	//if(!doToymodel || (doToymodel && nsuf>1))wrkdir+=Form("_%s",suff[sf].Data());
	TString filein = Form("%s/%s/unfolded_SignalSpectrum_R%.1lf_pTthresh%.1lf.root", wrkdir.Data(),prior_type[prior].Data(), C_Rs[r], C_pTls[ptl]);
	if(verbose) cout<<"opening file:"<<filein.Data()<<endl;
	TFile* funfolding = new TFile(filein.Data(), "OPEN");
	TString fileout = Form("%s/%s/unfoldingQA_R%.1lf_pTthresh%.1lf.root", wrkdir.Data(),prior_type[prior].Data(), C_Rs[r], C_pTls[ptl]);
	TFile* fqa = new TFile(fileout.Data(), "RECREATE");
	TFile* ftrue;
	TFile* fjets;

	const int vsize=maxiter[unf]; //size of the output vectors
	//chi2 test
	TVectorD chi2change(vsize); //succeessive iterations
	TVectorD chi2backf(vsize);
	TVectorD chi2backf_trunk(vsize);// backfolded vs measured

	//Kolmogorov Smirnov test
   TVectorD KSchange(vsize); //succeessive iterations
	TVectorD KSbackf(vsize);
	TVectorD KSbackf_trunk(vsize);// backfolded vs measured
	
	//Camberra test
	TVectorD Cambchange(vsize); //succeessive iterations
	TVectorD Cambbackf(vsize);
	TVectorD Cambbackf_trunk(vsize);// backfolded vs measured
	
	//Curvature test
	TVectorD Curvature(vsize);
	
	//toymodel closure test
	TVectorD Toyclosure_Camb(vsize);
	TVectorD Toyclosure_KS(vsize);
	TVectorD Toyclosure_chi2(vsize);
	
	TH1D* htrue;
	int nevtrue;
	int nevunf;
	if(doToymodel)
	{
		//unfolded vs true (TOYMODEl only)

		fjets=new TFile(Form("../../plotting_out/root/%s/%s/root_RAA%.1lf_%s/histos_jets_R%.1lf_pTcut0.2.root",resDir.Data(),uname.Data(),RAAtoy,uname.Data(),C_Rs[r]),"OPEN");
		TH1I* hevents=(TH1I*)fjets->Get("hevents");
		nevunf=hevents->GetEntries();
		cout<<"number of reconstructed events:"<<nevunf<<endl;
		delete hevents;
		fjets->Close();

		ftrue=new TFile(Form("../../plotting_out/root/%s/%s/root_RAA%.1lf_%s/jetonly_R%.1lf_pTcut0.2.root",resDir.Data(),uname.Data(),RAAtoy,uname.Data(),C_Rs[r]),"OPEN");
		TH1I* hevt=(TH1I*)ftrue->Get("hevents");
		nevtrue=hevt->GetEntries();
		cout<<"number of generated events:"<<nevtrue<<endl;
		delete hevt;
		
		funfolding->cd();
		TDirectoryFile* diter_1 = (TDirectoryFile*)funfolding->Get("iter1");
		TH1D* hunf_1=(TH1D*)diter_1->Get("hunfolded");
		
		ftrue->cd();
		TH2D *hPtRecpTleading = (TH2D*)ftrue->Get("fhPtRecpTleading");
		Int_t firstbin = hPtRecpTleading->GetXaxis()->FindBin(C_pTls[ptl]);
		Int_t lastbin = hPtRecpTleading->GetNbinsX();
		TH1D *htemp = (TH1D*) hPtRecpTleading->ProjectionY("htemp", firstbin, lastbin);
		float test_cont=htemp->GetBinContent(htemp->FindBin(10));
		cout<<"test content:"<<test_cont<<endl;
		htemp->Scale((float) nevunf/nevtrue/*,"width"*/);
		htrue = (TH1D*) rebinhisto(htemp,hunf_1,"hjettruth",0);
		/*htrue->SetLineColor(colorList[r+3]);
		htrue->SetLineStyle(2);
		htrue->DrawClone("same");*/
		cout<<"n events true: "<<nevtrue<<"n events unf: "<<nevunf<<endl;
		delete htemp;
		delete hPtRecpTleading;


	}
	for(int iter=1; iter<maxiter[unf]; iter++){ //first iteration: iter=0, we start from 2nd iteration (iter=1), since the first iteration is the same as the prior function (in SVD case)
		if(verbose) cout<<"iter:"<<iter<<endl;
						TDirectoryFile* dinput = (TDirectoryFile*)funfolding->Get("input");
						TDirectoryFile* diter_n = (TDirectoryFile*)funfolding->Get(Form("iter%i",iter));
						TDirectoryFile* diter_m = (TDirectoryFile*)funfolding->Get(Form("iter%i",iter-1));
						
						
						//sucessive iterations
						TH1D* hunf_n=(TH1D*)diter_n->Get("hunfolded");
						TH1D* hunf_m=(TH1D*)diter_m->Get("hunfolded");
						//cout<<"chi2 test"<<endl;
						chi2change[iter]=chi2_test(hunf_n,hunf_m,maxcutoff[system]);
						//cout<<"Cambera test"<<endl;
						Cambchange[iter]=Canberra_test(hunf_n,hunf_m,maxcutoff[system]);
						//cout<<" test"<<endl;
						KSchange[iter]=KS_test(hunf_n,hunf_m,maxcutoff[system]);
						
						//backfolded vs measured
						TH1D* hbf=(TH1D*) diter_n->Get("hbackfolded");
						TH1D* hmeasured=(TH1D*) dinput->Get("hmeasured");
		
						chi2backf[iter]=chi2_test(hbf,hmeasured,maxcutoff[system]);
						Cambbackf[iter]=Canberra_test(hbf,hmeasured,maxcutoff[system]);
						KSbackf[iter]=KS_test(hbf,hmeasured,maxcutoff[system]);
						
						//test of trunkated distributions
						chi2backf_trunk[iter]=chi2_test(hbf,hmeasured,maxcutoff[system],mincutoff[system]);
						Cambbackf_trunk[iter]=Canberra_test(hbf,hmeasured,maxcutoff[system],mincutoff[system]);
						KSbackf_trunk[iter]=KS_test(hbf,hmeasured,maxcutoff[system],mincutoff[system]);
						
						//curvature test
						Curvature[iter]=Curv_test(hunf_n,maxcutoff[system],mincutof_curv[system]);
						//cout<<iter<<" curvature:"<<Curvature[iter]<<endl;
						/*if(iter==4){
						hunf_n->SetLineColor(colorList[r]);
						hunf_n->DrawClone("same");}*/
						if(doToymodel) //toymodel closure test
						{
							//hunf_n->Scale(1.0,"width");
							float tc=Canberra_test(hunf_n,htrue,maxcutoff[system],5);
							float tk=KS_test(hunf_n,htrue,maxcutoff[system],5);
							float tx=chi2_test(hunf_n,htrue,maxcutoff[system],5);
							Toyclosure_Camb[iter]=tc;
							Toyclosure_KS[iter]=tk;
							Toyclosure_chi2[iter]=tx;
							//cout<<iter+1<<" Toy test: "<<tc<<endl;
						}

					}//iteration

               //write the output
               fqa->cd();
					
					chi2change.Write("chi2change");
					Cambchange.Write("Cambchange");
					KSchange.Write("KSchange");
					
					chi2backf.Write("chi2backf");
					Cambbackf.Write("Cambbackf");
					KSbackf.Write("KSbackf");
					
					chi2backf_trunk.Write("chi2backf_trunk");
					Cambbackf_trunk.Write("Cambbackf_trunk");
					KSbackf_trunk.Write("KSbackf_trunk");
					
					Curvature.Write("Curvature");
					
					if(doToymodel)
					{
						Toyclosure_Camb.Write("Toyclosure_Camb");
						Toyclosure_KS.Write("Toyclosure_KS");
						Toyclosure_chi2.Write("Toyclosure_chi2");
					}
						fqa->Close();
					funfolding->Close();
				if(doToymodel)ftrue->Close();
				}//prior
			}//unfolding
		}//pTlead
	}//R
	}//suffix
	}//binning

	
}
