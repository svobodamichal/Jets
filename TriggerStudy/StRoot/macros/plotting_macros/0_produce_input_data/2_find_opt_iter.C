#include <iostream>
#include <fstream>
#include "../utils/common.h"
#include "../utils/qa_cuts.h"
#include "../utils/utils.C"

using namespace std;

void find_opt_iter(C_systems system=cent, int doToymodel=0, TString suffix="GPC2",C_tests test=camb, float toytest_limit=0.10)
{

	//const int ntests=2;
	//const int first_test=1;
	//TString test_name[]={"chi2","Camb","KS"};
	//const int C_nunfoldings=2; //moved to utils/common.h
	//TString C_unfoldings_name[]={"Bayes","SVD"}; //moved to utils/common.h
	const int nprior=16;
	TString prior_type[]={"","flat","pythia","powlaw3","powlaw45","powlaw5",
		"powlaw55","tsalis_1","tsalis_2","tsalis_3","tsalis_4","tsalis_5","tsalis_6","tsalis_7","tsalis_8","tsalis_9"};
	int prior_start=2;
	int maxiter=8;
	TString nbins="VAR";
    bool ppHT=1; //0: use only MB data for pp | 1: use combined MB+HT data for pp

   const int nR=3; 
	const float Rs[]={0.2,0.3,0.4,0.5}; 
	const int npTlead=7;
   const float pTls[]={0,1,2,3,4,5,6,7}; 
    
	int nbinnings=3;
	if(doToymodel)nbinnings=1;
	const int binnings[]={1,4,0,2,3};

	int nsuf=17;  //13
	if(doToymodel)nsuf=1;
	const TString suff[]={"normal","pythiaU","pythiaG","pythiaUG","nfit13","nfit14","nfit18","nfit20","RRho02","RRho04","nrem-1", "pp","v2","m5","p5","u","g"};
	
	TString cent_name[]={"_central","_peripheral","_pp"};
    if(ppHT) cent_name[2]="HT_pp";
	int ptl_start=5;
	int ptl_end=7;
   if(system==peri)
   {
      ptl_start=4;
      ptl_end=7;
   }
	else if(system==pp)
	{
      ptl_start=3;
      ptl_end=6;
	}
	float RAAtoy=(system==cent) ? 0.5 : 0.7; //toymodel RAA of hard jet spectrum
   
	TString dtype="data";
	if(doToymodel>0)dtype="toy";
	
		
	for(int bn=0;bn<nbinnings;bn++)
	{
		int binning=binnings[bn];
		cout<<"binning:"<<binning<<endl;
		
	//for(int test=first_test; test<ntests; test++){	
	cout<<"test: "<<C_test_name[test].Data()<<endl;
	
	//control histograms
	TString fnamectrl=Form("./out_optIter/ctrl_histos%s_binning%i_%s.root",cent_name[system].Data(),binnings[bn],C_test_name[test].Data());
	if(doToymodel>0)fnamectrl=Form("./out_optIter/ctrl_histos%s_binning%i_%s_toy%i.root",cent_name[system].Data(),binnings[bn],C_test_name[test].Data(),doToymodel);
	TFile *fchistos=new TFile(fnamectrl,"RECREATE");
	TH1D* hoptiter[nR][npTlead][C_nunfoldings];
	TH1D* htest_bf[nR][npTlead][C_nunfoldings];
	TH1D* htest_suc[nR][npTlead][C_nunfoldings];
	TH1D* htest_cur[nR][npTlead][C_nunfoldings];
	
	TH2D* htest_bf2d[nR][npTlead][C_nunfoldings];
	TH2D* htest_suc2d[nR][npTlead][C_nunfoldings];
	TH2D* htest_cur2d[nR][npTlead][C_nunfoldings];

	bool doHistos=true;
	for(int sf=0; sf<nsuf;sf++)
	{
	if(sf>0)doHistos=false;
		//if(doToymodel>0 && sf>0)continue;
		cout<<sf+1<<"/"<<nsuf<<" - "<<suff[sf].Data()<<endl;
	ofstream fout;
	TString fname=Form("../utils/optIter/optIter_%s%s_%s_%s_bin%i.txt",dtype.Data(),cent_name[system].Data(),suff[sf].Data(),C_test_name[test].Data(),binning);
   fout.open(fname);
   fout << "//prior|R|pTlead|unfolding|centrality|binning|test"<< "\n";
   
	TString trigg_cent="MB";
	if(doToymodel>0) trigg_cent="toymodel";
	trigg_cent+=cent_name[system];
	TString toyraa="";
	if(doToymodel>0) toyraa=Form("_RAA%.1lf",RAAtoy);
	
	float xmin[/*ntests*/]={0,0,0};
	float xmax[]={100,0.5,0.05};
	int ntestbins[]={50,50,50};
	float ymax_scl[]={10,2,4};
	/*
	TString hname1[C_nunfoldings];
	TString hname2[C_nunfoldings];
	TString hname3[C_nunfoldings];
	
	for(int unf=0; unf<C_nunfoldings; unf++){
	hname1[unf]=Form("htest_bf_unf%i",unf);
	htest_bf[unf]=new TH1D(hname1[unf],"backfoled-measured test",50,xmin[test],xmax[test]);
	hname2[unf]=Form("htest_suc_unf%i",unf);
	htest_suc[unf]=new TH1D(hname2[unf],"succeessive iterations test",50,xmin[test],xmax[test]);
	hname3[unf]=Form("htest_cur_unf%i",unf);
	htest_cur[unf]=new TH1D(hname3[unf],"curvature test",20,0,10);
	}*/

	for(int r=0; r<nR; r++){
		if(r==0)continue;
		if(r==2)continue;
      fout << "//R=" << Rs[r]<< "\n";
      for(int ptl=ptl_start; ptl<=ptl_end; ptl++){
         fout << "//pTlead=" << pTls[ptl]<< "\n";
			for(int unf=0; unf<C_nunfoldings; unf++){
				
				//create histograms
				if(doHistos)
				{
				TString hname=Form("hiter_R0%.0lf_pTl%.0lf_unf%i",Rs[r]*10,pTls[ptl],unf);
				hoptiter[r][ptl][unf]=new TH1D(hname,"optimal iterations",10,-1.5,8.5);
				
				TString hname_bf=Form("htest_bf_R0%.0lf_pTl%.0lf_unf%i",Rs[r]*10,pTls[ptl],unf);
				htest_bf[r][ptl][unf]=new TH1D(hname_bf,"backfolded-measured test",ntestbins[test],xmin[test],xmax[test]);
				TString hname_suc=Form("htest_suc_R0%.0lf_pTl%.0lf_unf%i",Rs[r]*10,pTls[ptl],unf);
				htest_suc[r][ptl][unf]=new TH1D(hname_suc,"succeessive iter. test",ntestbins[test],xmin[test],xmax[test]);
				TString hname_cur=Form("htest_cur_R0%.0lf_pTl%.0lf_unf%i",Rs[r]*10,pTls[ptl],unf);
				htest_cur[r][ptl][unf]=new TH1D(hname_cur,"curvature test",100,0,100);
				
				TString hname_bf2d=Form("htest_bf_2D_R0%.0lf_pTl%.0lf_unf%i",Rs[r]*10,pTls[ptl],unf);
				htest_bf2d[r][ptl][unf]=new TH2D(hname_bf2d,"backfolded-measured test",ntestbins[test],xmin[test],xmax[test],ntestbins[test],xmin[test],ymax_scl[test]*xmax[test]);
				TString hname_suc2d=Form("htest_suc_2D_R0%.0lf_pTl%.0lf_unf%i",Rs[r]*10,pTls[ptl],unf);
				htest_suc2d[r][ptl][unf]=new TH2D(hname_suc2d,"succeessive iter. test",ntestbins[test],xmin[test],xmax[test],ntestbins[test],xmin[test],ymax_scl[test]*xmax[test]);
				TString hname_cur2d=Form("htest_cur_2D_R0%.0lf_pTl%.0lf_unf%i",Rs[r]*10,pTls[ptl],unf);
				htest_cur2d[r][ptl][unf]=new TH2D(hname_cur2d,"curvature test",100,0,100,ntestbins[1],xmin[test],ymax_scl[test]*xmax[test]);
				}
            fout << "//" << C_unfoldings_name[unf].Data()<< "\n";
				for(int prior=prior_start; prior<nprior; prior++){
               if (prior==3)continue;
	TString unfType=(system==pp) ? "dete" : "BGD";
	TString wrkdir = Form("../../plotting_out/root/%s/%s/Unfolded_R%.1lf_%s_%sbins_bining%i_%s%s_%s",trigg_cent.Data(),suffix.Data(), Rs[r],C_unfoldings_name[unf].Data(),nbins.Data(),binning,unfType.Data(),toyraa.Data(),suffix.Data());
	//if(!doToymodel || (doToymodel && nsuf>1))wrkdir+=Form("_%s",suff[sf].Data());
	TString filen = Form("%s/%s/unfoldingQA_R%.1lf_pTthresh%.1lf.root", wrkdir.Data(),prior_type[prior].Data(), Rs[r], pTls[ptl]);
	TFile* fqa = new TFile(filen.Data(), "OPEN");
	TVectorD *vec_iterch = (TVectorD*)fqa->Get(Form("%schange",C_test_name[test].Data())); //ratio of succeessive iterations
	TVectorD *vec_backf=(TVectorD*)fqa->Get(Form("%sbackf_trunk",C_test_name[test].Data())); //ratio backfolded/measured
	TVectorD *vec_curv=(TVectorD*)fqa->Get("Curvature"); //curvature of the unfolded spectrum
	TVectorD *vec_toy_Camb;
	TVectorD *vec_toy_KS;
	TVectorD *vec_toy_chi2;
	
	if(doToymodel>0) 
	{
		vec_toy_Camb=(TVectorD*)fqa->Get("Toyclosure_Camb"); //toymodel closure test value - distance of the unfolded solution from true
		vec_toy_KS=(TVectorD*)fqa->Get("Toyclosure_KS"); //toymodel closure test value - distance of the unfolded solution from true
		vec_toy_chi2=(TVectorD*)fqa->Get("Toyclosure_chi2"); //toymodel closure test value - distance of the unfolded solution from true
	}

	int optIter=-1;
	float chi2SVD[10];
	for(int iter=1; iter<maxiter; iter++){ //first iteration: iter=0, we start from 2nd iteration (iter=1), since the first iteration is the same as the prior function (in SVD case)
                  if(optIter>0) continue; // we already found a good iteration
                  float iter_change=vec_iterch(iter);
                  float backf_vs_meas=vec_backf(iter);
						float curvature=vec_curv(iter);
						float toytest[]={0,0,0};
							if(doToymodel>0) 
							{
								toytest[0]=vec_toy_chi2(iter);
								toytest[1]=vec_toy_Camb(iter);
								toytest[2]=vec_toy_KS(iter);
							}
						//apply quality cuts
						if(qa_check(iter_change,backf_vs_meas,curvature, system, unf, binning, test, r, ptl))optIter=iter+1;
               
					if(doHistos){
						if(doToymodel>0){			
						htest_bf2d[r][ptl][unf]->Fill(backf_vs_meas,toytest[test]);
						//cout<<"bf vs measured:"<<backf_vs_meas<<"toytest:"<<toytest[test]<<endl;
						htest_suc2d[r][ptl][unf]->Fill(iter_change,toytest[test]);
						htest_cur2d[r][ptl][unf]->Fill(curvature,toytest[test]);
						}
						if(doToymodel==1 && toytest[1]>toytest_limit) continue; //save only good toymodel solutions
						else if (doToymodel==2 && toytest[1]<=toytest_limit) continue; //save only bad toymodel solutions
						
						htest_bf[r][ptl][unf]->Fill(backf_vs_meas);
						htest_suc[r][ptl][unf]->Fill(iter_change);
						htest_cur[r][ptl][unf]->Fill(curvature);
						
					}//fill histograms	
						
					}//iteration

               TString line;
					if(optIter>0) 
					{
						line=Form("diter_%s%s_%s[%i][%i][%i][%i][%i][%i][%i]=%i; //backf_vs_meas=%.1lf, iter_change=%.1lf, curvature=%.0lf",dtype.Data(),cent_name[system].Data(),suff[sf].Data(),prior,r,ptl,unf,system,binning,test,optIter,vec_backf(optIter-1),vec_iterch(optIter-1),vec_curv(optIter-1));
					}
					else line=Form("diter_%s%s_%s[%i][%i][%i][%i][%i][%i][%i]=%i;",dtype.Data(),cent_name[system].Data(),suff[sf].Data(),prior,r,ptl,unf,system,binning,test,optIter);
					//write out the code
					fout << line.Data() << "\n";
					fqa->Close();
					
					if(doHistos) hoptiter[r][ptl][unf]->Fill(optIter);
				}//prior

			if(doHistos){			
			fchistos->cd();
			hoptiter[r][ptl][unf]->Write(hname);
			htest_bf[r][ptl][unf]->Write(hname_bf);
			htest_suc[r][ptl][unf]->Write(hname_suc);
			htest_cur[r][ptl][unf]->Write(hname_cur);
			
			if(doToymodel>0)
			{
			htest_bf2d[r][ptl][unf]->Write(hname_bf2d);
			htest_suc2d[r][ptl][unf]->Write(hname_suc2d);
			htest_cur2d[r][ptl][unf]->Write(hname_cur2d);
			}
			}//save histograms			
			}//unfolding
		}//pTlead
		}//R
		}//suffix
		fchistos->Close();
	//}//test
	}//binning
	fout.close();

}
