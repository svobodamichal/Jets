//This macro does:
//1) plots number of good unfolded solutions for different Rs and pTleadings
//2) using toymodel simulation input, it determines and prints out recommended values of QA cuts used to find optimal iterations (these cut values are stored in ../utils/qa_cuts.h
//PREREQUIREMENTS: 
//1) you need to run ../unfold_qa.C on unfolded data, in order to calculate test statistics (e.g. relative distance of successive iterations)
//2) you need to run ../find_optiter.C, which produces control histograms used as input for this macro (find_optiter.C also calculates optimal iterations, but for this you need first to set the cut values in ../utils/qa_cuts.C)


#include "../../utils/common.h"
#include "../../utils/qa_cuts.h"


void show_iter(C_tests test=camb,int bin=1,C_systems system=cent,bool toymodel=0,TString grsuff="gif"/*,int r_show=1, int pTl_show=2*/) 
{
	float distanceCut[C_nR]={0.25,0.25,0.35};
	//const int ntests=3;
	//const TString C_test_name[C_ntests]={"chi2","Camb","KS"};
	const TString test_name2[C_ntests]={"#chi^{2} test","relat. distance","Kol.-Smir."};
	//const int nunfoldings=2; //number of unfoldings to use //moved to utils/common.h
	//const TString C_unfoldings_name[C_nunfoldings]={"Bayes","SVD"};//moved to utils/common.h
	//const int nR=3;
	//const float Rs[]={0.2,0.3,0.4,0.5};
	TString sys_name[3]={"","_peripheral","_pp"};
	int ptl_start=5;
   int ptl_end=7;
   if(system==peri)
   {
      ptl_start=4;
      ptl_end=6;
   }
   else if(system==pp)
   {
      ptl_start=3;
      ptl_end=5;
   }
   int nptl=ptl_end-ptl_start+1;
	//const float pTls[]={0,1,3,4,5,6,7};
	
	TString fname;
	TFile *fin;
	ofstream fout;
	fout.open("new_qa_cuts.h");

	//input file
	fname=Form("input/ctrl_histos%s_binning%i_%s.root",sys_name[system].Data(),bin,C_test_name[test].Data());
	if(toymodel) fname=Form("input/ctrl_histos%s_binning%i_%s_toy1.root",sys_name[system].Data(),bin,C_test_name[test].Data());
	cout<<"opening file:"<<fname.Data()<<endl;
	fin=new TFile(fname,"OPEN");

	TH1D* hcontrol[C_nR][C_npTlead][C_nunfoldings];
	Color_t hcolor[C_nunfoldings]={kRed,kBlue};
	Color_t test_color[]={kBlue,kRed};
	TString drawcmd[C_nunfoldings]={"","same"};
	
	TH1D* htest_bf[C_nR][C_npTlead][C_nunfoldings];
	TH1D* htest_suc[C_nR][C_npTlead][C_nunfoldings];
	TH1D* htest_cur[C_nR][C_npTlead][C_nunfoldings];
	
	TH2D* htest_bf2d[C_nR][C_npTlead][C_nunfoldings];
	TH2D* htest_suc2d[C_nR][C_npTlead][C_nunfoldings];
	TH2D* htest_cur2d[C_nR][C_npTlead][C_nunfoldings];
	
	float critical_value_bf[C_ntests][C_nsystems][C_nbinnings][C_nunfoldings][C_nR][C_npTlead];
	float critical_value_suc[C_ntests][C_nsystems][C_nbinnings][C_nunfoldings][C_nR][C_npTlead];
	
	TLatex *latex = new TLatex();
	latex->SetNDC();
	latex->SetTextSize(0.04);

	const int can_x=400;
	const int can_y=400;
	TCanvas* c1=new TCanvas("c1","iterations",10,10,C_nR*can_x,nptl*can_y);
	c1->Divide(C_nR,nptl);
	int counter=1;
	for(int ptl=ptl_start; ptl<=ptl_end; ptl++)
	{
		for(int r=0; r<C_nR; r++)
		{
			c1->cd(counter);
			counter++;
			for(int unf=0; unf<C_nunfoldings; unf++)
			{
				TString hname=Form("hiter_R0%.0lf_pTl%.0lf_unf%i",C_Rs[r]*10,C_pTls[ptl],unf);
				hcontrol[r][ptl][unf]=(TH1D*) fin->Get(hname);
				hcontrol[r][ptl][unf]->SetLineColor(hcolor[unf]);
				hcontrol[r][ptl][unf]->Draw(drawcmd[unf]);
				latex->DrawLatex(0.35, 0.8,Form("R=%.1lf, pTlead>%.0lf",C_Rs[r],C_pTls[ptl]));
			}//unf
		}//r
	}//pTlead	
	
	TString ytitle[C_ntests]={"#chi^{2}(true,unfolded)", "R(true,unfolded)", "#Delta_{KS}(true,unfolded)"};
		

	TCanvas* c_bf[C_nunfoldings];
	TCanvas* c_suc[C_nunfoldings];
	TCanvas* c_cur[C_nunfoldings];
	
	for(int unf=0; unf<C_nunfoldings; unf++)
	{
	if(toymodel)
	{
	c_bf[unf]=new TCanvas(Form("c_bf_%i",unf),Form("backfolded_measured_%s",C_unfoldings_name[unf].Data()),10,10,C_nR*can_x,nptl*can_y);
	
	c_bf[unf]->Divide(C_nR,nptl);
	int counter=1;
	for(int ptl=ptl_start; ptl<=ptl_end; ptl++)
	{
		for(int r=0; r<C_nR; r++)
		{
			c_bf[unf]->cd(counter);
			gPad->SetLeftMargin(0.13);
			counter++;

				TString hname=Form("htest_bf_2D_R0%.0lf_pTl%.0lf_unf%i",C_Rs[r]*10,C_pTls[ptl],unf);
				htest_bf2d[r][ptl][unf]=(TH2D*) fin->Get(hname);
				//htest_bf2d[r][ptl][unf]->SetLineColor(hcolor[unf]);
				htest_bf2d[r][ptl][unf]->Draw("COLZ");
				htest_bf2d[r][ptl][unf]->GetYaxis()->SetTitle(ytitle[test]);
				htest_bf2d[r][ptl][unf]->GetYaxis()->SetTitleSize(0.04);
				htest_bf2d[r][ptl][unf]->GetYaxis()->SetTitleOffset(1.5);
				//latex->DrawLatex(0.35, 0.8,"backfolded/measured");
				latex->DrawLatex(0.35, 0.8,Form("%s, %s",C_unfoldings_name[unf].Data(),test_name2[test].Data()));
				latex->DrawLatex(0.35, 0.7,Form("R=%.1lf, pTlead>%.0lf",C_Rs[r],C_pTls[ptl]));
				
				//***********
				//FIT
				//***********
				float a=0;
				float b=0;
				TF1* fitfnc=fit_histo(htest_bf2d[r][ptl][unf]);
				a=fitfnc->GetParameter(0);
				b=fitfnc->GetParameter(1);
				fitfnc->DrawClone("same");
				float optCut=(a>0) ? ((distanceCut[r]-b)/a) : 0;
				latex->DrawLatex(0.35, 0.3,Form("opt. cut:%.2lf",optCut));
				delete fitfnc;
				
				critical_value_bf[test][system][bin][unf][r][ptl]=optCut;
				
		}//r
	}//pTlead	
	c_bf[unf]->SaveAs(Form("./obr/toy/test_bf_%s_unf%i_bin%i.%s",C_test_name[test].Data(),unf,bin,grsuff.Data()));
	
	c_suc[unf]=new TCanvas(Form("c_suc_%i",unf),Form("successive_iterations_%s",C_unfoldings_name[unf].Data()),10,10,C_nR*can_x,nptl*can_y);
	c_suc[unf]->Divide(C_nR,nptl);
	gPad->SetLeftMargin(0.15);
	int counter=1;
	for(int ptl=ptl_start; ptl<=ptl_end; ptl++)
	{
		for(int r=0; r<C_nR; r++)
		{
			c_suc[unf]->cd(counter);
			counter++;
		
				TString hname=Form("htest_suc_2D_R0%.0lf_pTl%.0lf_unf%i",C_Rs[r]*10,C_pTls[ptl],unf);
				htest_suc2d[r][ptl][unf]=(TH2D*) fin->Get(hname);
				//htest_suc2d[r][ptl][unf]->SetLineColor(hcolor[unf]);
				htest_suc2d[r][ptl][unf]->GetYaxis()->SetTitle(ytitle[test]);
				htest_suc2d[r][ptl][unf]->GetYaxis()->SetTitleSize(0.04);
				htest_suc2d[r][ptl][unf]->GetYaxis()->SetTitleOffset(1.2);
				htest_suc2d[r][ptl][unf]->Draw("COLZ");
				latex->DrawLatex(0.35, 0.8,Form("%s, %s",C_unfoldings_name[unf].Data(),test_name2[test].Data()));
				latex->DrawLatex(0.35, 0.7,Form("R=%.1lf, pTlead>%.0lf",C_Rs[r],C_pTls[ptl]));
				
				//***********
				//FIT
				//***********
				float a=0;
				float b=0;
				TF1* fitfnc=fit_histo(htest_suc2d[r][ptl][unf]);
				a=fitfnc->GetParameter(0);
				b=fitfnc->GetParameter(1);
				fitfnc->DrawClone("same");
				float optCut=(a>0) ? ((distanceCut[r]-b)/a) : 0;
				latex->DrawLatex(0.35, 0.3,Form("opt. cut:%.2lf",optCut));
				delete fitfnc;
				
                if(unf==SVD)optCut=1;
				critical_value_suc[test][system][bin][unf][r][ptl]=optCut;
				
		}//r
	}//pTlead	
	c_suc[unf]->SaveAs(Form("./obr/toy/test_suc_%s_unf%i_bin%i.%s",C_test_name[test].Data(),unf,bin,grsuff.Data()));
	
	c_cur[unf]=new TCanvas(Form("c_cur_%i",unf),Form("Curvature test, %s",C_unfoldings_name[unf].Data()),10,10,C_nR*can_x,nptl*can_y);
	c_cur[unf]->Divide(C_nR,nptl);
	int counter=1;
	for(int ptl=ptl_start; ptl<=ptl_end; ptl++)
	{
		for(int r=0; r<C_nR; r++)
		{
			c_cur[unf]->cd(counter);
			gPad->SetLeftMargin(0.13);
			counter++;

				TString hname=Form("htest_cur_2D_R0%.0lf_pTl%.0lf_unf%i",C_Rs[r]*10,C_pTls[ptl],unf);
				htest_cur2d[r][ptl][unf]=(TH2D*) fin->Get(hname);
				//htest_cur2d[r][ptl][unf]->SetLineColor(hcolor[unf]);
				htest_cur2d[r][ptl][unf]->GetYaxis()->SetTitle(ytitle[test]);
				htest_cur2d[r][ptl][unf]->GetYaxis()->SetTitleSize(0.04);
				htest_cur2d[r][ptl][unf]->GetYaxis()->SetTitleOffset(1.5);
				htest_cur2d[r][ptl][unf]->Draw("COLZ");
				latex->DrawLatex(0.35, 0.8,Form("%s, curvature",C_unfoldings_name[unf].Data()));
				latex->DrawLatex(0.35, 0.7,Form("R=%.1lf, pTlead>%.0lf",C_Rs[r],C_pTls[ptl]));
		}//r
	}//pTlead
	c_cur[unf]->SaveAs(Form("./obr/toy/test_cur_%s_unf%i_bin%i.%s",C_test_name[test].Data(),unf,bin,grsuff.Data()));

	}//do toymodel
	else //data
	{
	
	

	c_bf[unf]=new TCanvas(Form("c_bf_%i",unf),Form("backfolded_measured_%s",C_unfoldings_name[unf].Data()),10,10,C_nR*can_x,nptl*can_y);
	c_bf[unf]->Divide(C_nR,nptl);
	int counter=1;
	for(int ptl=ptl_start; ptl<=ptl_end; ptl++)
	{
		for(int r=0; r<C_nR; r++)
		{
			c_bf[unf]->cd(counter);
			gPad->SetLeftMargin(0.13);
			counter++;

				TString hname=Form("htest_bf_R0%.0lf_pTl%.0lf_unf%i",C_Rs[r]*10,C_pTls[ptl],unf);
				htest_bf[r][ptl][unf]=(TH1D*) fin->Get(hname);
				//htest_bf2d[r][ptl][unf]->SetLineColor(hcolor[unf]);
				htest_bf[r][ptl][unf]->Draw();
				htest_bf[r][ptl][unf]->GetYaxis()->SetTitle("");
				htest_bf[r][ptl][unf]->GetYaxis()->SetTitleSize(0.04);
				htest_bf[r][ptl][unf]->GetYaxis()->SetTitleOffset(1.5);
				//latex->DrawLatex(0.35, 0.8,"backfolded/measured");
				latex->DrawLatex(0.35, 0.8,Form("%s, %s",C_unfoldings_name[unf].Data(),test_name2[test].Data()));
				latex->DrawLatex(0.35, 0.7,Form("R=%.1lf, pTlead>%.0lf",C_Rs[r],C_pTls[ptl]));
				
               	float xcutbf=backfold_cut[test][system][bin][unf][r][ptl];
                //fout<<"xcutbf:"<<xcutbf<<endl;
				float max=htest_bf[r][ptl][unf]->GetMaximum();
				TLine *one = new TLine(xcutbf, 0, xcutbf, max);
				one->SetLineWidth(2);
				one->SetLineStyle(2);
				one->SetLineColor(kBlack);
				one->DrawClone("same");
		}//r
	}//pTlead	
	c_bf[unf]->SaveAs(Form("./obr/data/test_bf_%s_unf%i_bin%i.%s",C_test_name[test].Data(),unf,bin,grsuff.Data()));
	
	c_suc[unf]=new TCanvas(Form("c_suc_%i",unf),Form("successive_iterations_%s",C_unfoldings_name[unf].Data()),10,10,C_nR*can_x,nptl*can_y);
	c_suc[unf]->Divide(C_nR,nptl);
	gPad->SetLeftMargin(0.15);
	int counter=1;
	for(int ptl=ptl_start; ptl<=ptl_end; ptl++)
	{
		for(int r=0; r<C_nR; r++)
		{
			c_suc[unf]->cd(counter);
			counter++;
		
				TString hname=Form("htest_suc_R0%.0lf_pTl%.0lf_unf%i",C_Rs[r]*10,C_pTls[ptl],unf);
				htest_suc[r][ptl][unf]=(TH1D*) fin->Get(hname);
				//htest_suc2d[r][ptl][unf]->SetLineColor(hcolor[unf]);
				htest_suc[r][ptl][unf]->GetYaxis()->SetTitle("");
				htest_suc[r][ptl][unf]->GetYaxis()->SetTitleSize(0.04);
				htest_suc[r][ptl][unf]->GetYaxis()->SetTitleOffset(1.2);
				htest_suc[r][ptl][unf]->Draw("");
				latex->DrawLatex(0.35, 0.8,Form("%s, %s",C_unfoldings_name[unf].Data(),test_name2[test].Data()));
				latex->DrawLatex(0.35, 0.7,Form("R=%.1lf, pTlead>%.0lf",C_Rs[r],C_pTls[ptl]));
				
               	float xcutsuc=change_cut[test][system][bin][unf][r][ptl];
				float max=htest_suc[r][ptl][unf]->GetMaximum();
				TLine *one = new TLine(xcutsuc, 0, xcutsuc, max);
				one->SetLineWidth(2);
				one->SetLineStyle(2);
				one->SetLineColor(kBlack);
				one->DrawClone("same");
		}//r
	}//pTlead	
	c_suc[unf]->SaveAs(Form("./obr/data/test_suc_%s_unf%i_bin%i.%s",C_test_name[test].Data(),unf,bin,grsuff.Data()));
	
	c_cur[unf]=new TCanvas(Form("c_cur_%i",unf),Form("Curvature test, %s",C_unfoldings_name[unf].Data()),10,10,C_nR*can_x,nptl*can_y);
	c_cur[unf]->Divide(C_nR,nptl);
	int counter=1;
	for(int ptl=ptl_start; ptl<=ptl_end; ptl++)
	{
		for(int r=0; r<C_nR; r++)
		{
			c_cur[unf]->cd(counter);
			gPad->SetLeftMargin(0.13);
			counter++;

				TString hname=Form("htest_cur_R0%.0lf_pTl%.0lf_unf%i",C_Rs[r]*10,C_pTls[ptl],unf);
				htest_cur[r][ptl][unf]=(TH1D*) fin->Get(hname);
				//htest_cur2d[r][ptl][unf]->SetLineColor(hcolor[unf]);
				htest_cur[r][ptl][unf]->GetYaxis()->SetTitle("");
				htest_cur[r][ptl][unf]->GetYaxis()->SetTitleSize(0.04);
				htest_cur[r][ptl][unf]->GetYaxis()->SetTitleOffset(1.5);
				htest_cur[r][ptl][unf]->Draw("");
				latex->DrawLatex(0.35, 0.8,Form("%s, curvature",C_unfoldings_name[unf].Data()));
				latex->DrawLatex(0.35, 0.7,Form("R=%.1lf, pTlead>%.0lf",C_Rs[r],C_pTls[ptl]));
				
               	float xcutcur=curv_cut_default[unf];
				float max=htest_cur[r][ptl][unf]->GetMaximum();
				TLine *one = new TLine(xcutcur, 0, xcutcur, max);
				one->SetLineWidth(2);
				one->SetLineStyle(2);
				one->SetLineColor(kBlack);
				one->DrawClone("same");
		}//r
	}//pTlead
	c_cur[unf]->SaveAs(Form("./obr/data/test_cur_%s_unf%i_bin%i.%s",C_test_name[test].Data(),unf,bin,grsuff.Data()));

	}//data
	
	//print out recommended critical values of QA tests
	if(toymodel){
	if(unf==0){
	fout<<"//test:"<<C_test_name[test].Data()<<endl;
	fout<<"//system:"<<C_system_name[system].Data()<<endl;
	fout<<"//bin:"<<bin<<endl;}
	fout<<"//unfolding:"<<C_unfoldings_name[unf].Data()<<endl;
	fout<<"//backfolded vs measured"<<endl;
	
	
	for(int r=0; r<C_nR; r++){
		fout<<"//R="<<C_Rs[r]<<", distance cut="<<distanceCut[r]<<endl;
	for(int ptl=ptl_start; ptl<=ptl_end; ptl++){
		fout<<"backfold_cut["<<test<<"]["<<system<<"]["<<bin<<"]["<<unf<<"]["<<r<<"]["<<ptl<<"]="<<critical_value_bf[test][system][bin][unf][r][ptl]<<";"<<endl;
	}//pTlead
	}//r
	
	fout<<"//successive iterations"<<endl;
	for(int r=0; r<C_nR; r++){
		fout<<"//R="<<C_Rs[r]<<", distance cut="<<distanceCut[r]<<endl;
	for(int ptl=ptl_start; ptl<=ptl_end; ptl++){
		fout<<"change_cut["<<test<<"]["<<system<<"]["<<bin<<"]["<<unf<<"]["<<r<<"]["<<ptl<<"]="<<critical_value_suc[test][system][bin][unf][r][ptl]<<";"<<endl;
	}//pTlead
	}//r
    }//toymodel
	
	}//unfolding
	
	fout.close();
	/*
	TCanvas* c2=new TCanvas("c2","tests",10,10,nunfoldings*can_x,3*can_y);
	c2->Divide(nunfoldings,3);
	
	for(int unf=0; unf<C_nunfoldings; unf++)
	{
		for(int pl=0; pl<nplots; pl++)
		{
				TString hname=Form("htest_bf_R0%.0lf_pTl%.0lf_unf%i",C_Rs[r_show]*10,C_pTls[pTl_show],unf);
				htest_bf[unf][pl]=(TH1D*) fin[pl]->Get(hname);
				htest_bf[unf][pl]->Rebin(2);
				hname=Form("htest_suc_R0%.0lf_pTl%.0lf_unf%i",C_Rs[r_show]*10,C_pTls[pTl_show],unf);
				htest_suc[unf][pl]=(TH1D*) fin[pl]->Get(hname);
				htest_suc[unf][pl]->Rebin(2);
				hname=Form("htest_cur_R0%.0lf_pTl%.0lf_unf%i",C_Rs[r_show]*10,C_pTls[pTl_show],unf);
				htest_cur[unf][pl]=(TH1D*) fin[pl]->Get(hname);
				
 
				c2->cd(1+unf);
				htest_bf[unf][pl]->SetLineColor(test_color[pl]);
				htest_bf[unf][pl]->Draw(drawcmd[pl]);
				if(unf==0)latex->DrawLatex(0.9, 0.95,Form("%s",C_test_name[test].Data()));
				latex->DrawLatex(0.45, 0.75,Form("%s",C_unfoldings_name[unf].Data()));
				//latex->DrawLatex(0.35, 0.7,"backfolded/measured");
				
				float xcut=backfold_cut[system][test][unf];
				float max=htest_bf[unf][pl]->GetMaximum();
				  TLine *one = new TLine(xcut, 0, xcut, max);
					one->SetLineWidth(2);
					one->SetLineStyle(2);
					one->SetLineColor(kBlack);
					one->DrawClone("same");
				
				c2->cd(3+unf);
				htest_suc[unf][pl]->SetLineColor(test_color[pl]);
				htest_suc[unf][pl]->Draw(drawcmd[pl]);
				//latex->DrawLatex(0.35, 0.7,"successive iterations");
				
				xcut=change_cut[system][test][unf];
				max=htest_suc[unf][pl]->GetMaximum();
				  TLine *one = new TLine(xcut, 0, xcut, max);
					one->SetLineWidth(2);
					one->SetLineStyle(2);
					one->SetLineColor(kBlack);
					one->DrawClone("same");
				
				c2->cd(5+unf);
				htest_cur[unf][pl]->SetLineColor(test_color[pl]);
				htest_cur[unf][pl]->Draw(drawcmd[pl]);
				//latex->DrawLatex(0.35, 0.7,"curvature");
				
				xcut=curv_cut[unf];
				max=htest_cur[unf][pl]->GetMaximum();
				  TLine *one = new TLine(xcut, 0, xcut, max);
					one->SetLineWidth(2);
					one->SetLineStyle(2);
					one->SetLineColor(kBlack);
					one->DrawClone("same");
				
			}//toymodel type [good | bad]
	}//unfolding
		c2->SaveAs(Form("./tests%i_toy%i.%s",test,toymodel,grsuff.Data()));*/
}//This is the end, my friend...

//function for fitting 2D histograms
TF1* fit_histo(TH2D* histo)
{
	int nbins_tmp=histo->GetNbinsX()*histo->GetNbinsY();
	//fout<<"nbins:"<<nbins_tmp<<endl;
	const int nbinsxy=2500; //nbins_tmp; 
	//fout<<"nbins:"<<nbinsxy<<endl;
	float x_val[nbinsxy];
	float y_val[nbinsxy];
	float xmax=0;
	int k=0;
	for(int i=1; i<=histo->GetNbinsX(); i++){
		for(int j=1; j<=histo->GetNbinsY(); j++){
			float xc=histo->GetXaxis()->GetBinCenter(i);
			float yc=histo->GetYaxis()->GetBinCenter(j);
			float val=histo->GetBinContent(i,j);
			if(val>0)
			{
				if(xc>xmax) xmax=xc;
				x_val[k]=xc;
				y_val[k]=yc;
				k++;
			}
		}//j
	}//i
			
	TGraph* graph=new TGraph(k, x_val, y_val);
	TF1* ffit=new TF1("ffit","[0]*x+[1]",0,xmax);
	ffit->SetParameters(1,0);
	graph->Fit("ffit");
	delete graph;
	return ffit;
}
			
			
			
			
			
			
			
			
			
			
			
		
