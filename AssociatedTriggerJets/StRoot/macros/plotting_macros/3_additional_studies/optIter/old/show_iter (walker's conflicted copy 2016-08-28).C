#include "../utils/qa_cuts.h"

void show_iter(int test=0,int bin=1,bool peripheral=0,bool toymodel=0,TString grsuff="pdf",int unf_show=0,int r_show=1, int pTl_show=2) 
{
	const int ntests=3;
	const TString test_name[ntests]={"chi2","Camb","KS"};
	const TString test_name2[ntests]={"#chi^{2} test","relat. distance","Kol.-Smir."};
	const int nunf=2;
	const TString unfname[nunf]={"Bayes","SVD"};
	const int nr=3;
	const float R[]={0.2,0.3,0.4,0.5};
	TString cent_name[2]={"","_peripheral"};
	int ptl_start=2;
   int ptl_end=4;
   if(peripheral)
   {
      ptl_start=1;
      ptl_end=3;
   }
   int nptl=ptl_end-ptl_start+1;
	const float pTthresh[]={3,4,5,6,7};
	

	TString fname;
	TFile *fin;

	//input file
	fname=Form("ctrl_histos%s_binning%i_%s.root",cent_name[peripheral].Data(),bin,test_name[test].Data());
	if(toymodel) fname=Form("ctrl_histos%s_binning%i_%s_toy1.root",cent_name[peripheral].Data(),bin,test_name[test].Data());
	fin=new TFile(fname,"OPEN");

	TH1D* hcontrol[nr][5][2];
	Color_t hcolor[nunf]={kRed,kBlue};
	Color_t test_color[]={kBlue,kRed};
	TString drawcmd[nunf]={"","same"};
	
	TH1D* htest_bf[nr][5][nunf];
	TH1D* htest_suc[nr][5][nunf];
	TH1D* htest_cur[nr][5][nunf];
	
	TH2D* htest_bf2d[nr][5][nunf];
	TH2D* htest_suc2d[nr][5][nunf];
	TH2D* htest_cur2d[nr][5][nunf];
	
	TLatex *latex = new TLatex();
	latex->SetNDC();
	latex->SetTextSize(0.04);

	const int can_x=400;
	const int can_y=400;
	TCanvas* c1=new TCanvas("c1","iterations",10,10,nr*can_x,nptl*can_y);
	c1->Divide(nr,nptl);
	int counter=1;
	for(int ptl=ptl_start; ptl<=ptl_end; ptl++)
	{
		for(int r=0; r<nr; r++)
		{
			c1->cd(counter);
			counter++;
			for(int unf=0; unf<nunf; unf++)
			{
				TString hname=Form("hiter_R0%.0lf_pTl%.0lf_unf%i",R[r]*10,pTthresh[ptl],unf);
				hcontrol[r][ptl][unf]=(TH1D*) fin->Get(hname);
				hcontrol[r][ptl][unf]->SetLineColor(hcolor[unf]);
				hcontrol[r][ptl][unf]->Draw(drawcmd[unf]);
				latex->DrawLatex(0.35, 0.8,Form("R=%.1lf, pTlead>%.0lf",R[r],pTthresh[ptl]));
			}//unf
		}//r
	}//pTlead	
	
	TString ytitle[ntests]={"#chi^{2}(true,unfolded)", "R(true,unfolded)", "#Delta_{KS}(true,unfolded)"};
		
	if(toymodel)
	{
		
	TCanvas* c_bf=new TCanvas("c_bf","backfolded_measured",10,10,nr*can_x,nptl*can_y);
	c_bf->Divide(nr,nptl);
	int counter=1;
	for(int ptl=ptl_start; ptl<=ptl_end; ptl++)
	{
		for(int r=0; r<nr; r++)
		{
			c_bf->cd(counter);
			gPad->SetLeftMargin(0.13);
			counter++;

				TString hname=Form("htest_bf_2D_R0%.0lf_pTl%.0lf_unf%i",R[r]*10,pTthresh[ptl],unf_show);
				htest_bf2d[r][ptl][unf_show]=(TH2D*) fin->Get(hname);
				//htest_bf2d[r][ptl][unf_show]->SetLineColor(hcolor[unf]);
				htest_bf2d[r][ptl][unf_show]->Draw("COLZ");
				htest_bf2d[r][ptl][unf_show]->GetYaxis()->SetTitle(ytitle[test]);
				htest_bf2d[r][ptl][unf_show]->GetYaxis()->SetTitleSize(0.04);
				htest_bf2d[r][ptl][unf_show]->GetYaxis()->SetTitleOffset(1.5);
				//latex->DrawLatex(0.35, 0.8,"backfolded/measured");
				latex->DrawLatex(0.35, 0.8,Form("%s, %s",unfname[unf_show].Data(),test_name2[test].Data()));
				latex->DrawLatex(0.35, 0.7,Form("R=%.1lf, pTlead>%.0lf",R[r],pTthresh[ptl]));
		}//r
	}//pTlead	
	c_bf->SaveAs(Form("./obr/toy/test_bf_%s_unf%i_bin%i.%s",test_name[test].Data(),unf_show,bin,grsuff.Data()));
	
	TCanvas* c_suc=new TCanvas("c_suc","successive_iterations",10,10,nr*can_x,nptl*can_y);
	c_suc->Divide(nr,nptl);
	gPad->SetLeftMargin(0.15);
	int counter=1;
	for(int ptl=ptl_start; ptl<=ptl_end; ptl++)
	{
		for(int r=0; r<nr; r++)
		{
			c_suc->cd(counter);
			counter++;
		
				TString hname=Form("htest_suc_2D_R0%.0lf_pTl%.0lf_unf%i",R[r]*10,pTthresh[ptl],unf_show);
				htest_suc2d[r][ptl][unf_show]=(TH2D*) fin->Get(hname);
				//htest_suc2d[r][ptl][unf_show]->SetLineColor(hcolor[unf]);
				htest_suc2d[r][ptl][unf_show]->GetYaxis()->SetTitle(ytitle[test]);
				htest_suc2d[r][ptl][unf_show]->GetYaxis()->SetTitleSize(0.04);
				htest_suc2d[r][ptl][unf_show]->GetYaxis()->SetTitleOffset(1.2);
				htest_suc2d[r][ptl][unf_show]->Draw("COLZ");
				latex->DrawLatex(0.35, 0.8,Form("%s, %s",unfname[unf_show].Data(),test_name2[test].Data()));
				latex->DrawLatex(0.35, 0.7,Form("R=%.1lf, pTlead>%.0lf",R[r],pTthresh[ptl]));
		}//r
	}//pTlead	
	c_suc->SaveAs(Form("./obr/toy/test_suc_%s_unf%i_bin%i.%s",test_name[test].Data(),unf_show,bin,grsuff.Data()));
	
	TCanvas* c_cur=new TCanvas("c_cur","Curvature test",10,10,nr*can_x,nptl*can_y);
	c_cur->Divide(nr,nptl);
	int counter=1;
	for(int ptl=ptl_start; ptl<=ptl_end; ptl++)
	{
		for(int r=0; r<nr; r++)
		{
			c_cur->cd(counter);
			gPad->SetLeftMargin(0.13);
			counter++;

				TString hname=Form("htest_cur_2D_R0%.0lf_pTl%.0lf_unf%i",R[r]*10,pTthresh[ptl],unf_show);
				htest_cur2d[r][ptl][unf_show]=(TH2D*) fin->Get(hname);
				//htest_cur2d[r][ptl][unf_show]->SetLineColor(hcolor[unf]);
				htest_cur2d[r][ptl][unf_show]->GetYaxis()->SetTitle(ytitle[test]);
				htest_cur2d[r][ptl][unf_show]->GetYaxis()->SetTitleSize(0.04);
				htest_cur2d[r][ptl][unf_show]->GetYaxis()->SetTitleOffset(1.5);
				htest_cur2d[r][ptl][unf_show]->Draw("COLZ");
				latex->DrawLatex(0.35, 0.8,Form("%s, curvature",unfname[unf_show].Data()));
				latex->DrawLatex(0.35, 0.7,Form("R=%.1lf, pTlead>%.0lf",R[r],pTthresh[ptl]));
		}//r
	}//pTlead
	c_cur->SaveAs(Form("./obr/toy/test_cur_%s_unf%i_bin%i.%s",test_name[test].Data(),unf_show,bin,grsuff.Data()));

	
	}//do toymodel
	else //data
	{
	
	float xcutbf=backfold_cut[peripheral][test][unf_show];
	float xcutsuc=change_cut[peripheral][test][unf_show];
	float xcutcur=curv_cut[bin][unf_show];
	

	TCanvas* c_bf=new TCanvas("c_bf","backfolded_measured",10,10,nr*can_x,nptl*can_y);
	c_bf->Divide(nr,nptl);
	int counter=1;
	for(int ptl=ptl_start; ptl<=ptl_end; ptl++)
	{
		for(int r=0; r<nr; r++)
		{
			c_bf->cd(counter);
			gPad->SetLeftMargin(0.13);
			counter++;

				TString hname=Form("htest_bf_R0%.0lf_pTl%.0lf_unf%i",R[r]*10,pTthresh[ptl],unf_show);
				htest_bf[r][ptl][unf_show]=(TH1D*) fin->Get(hname);
				//htest_bf2d[r][ptl][unf_show]->SetLineColor(hcolor[unf]);
				htest_bf[r][ptl][unf_show]->Draw();
				htest_bf[r][ptl][unf_show]->GetYaxis()->SetTitle("");
				htest_bf[r][ptl][unf_show]->GetYaxis()->SetTitleSize(0.04);
				htest_bf[r][ptl][unf_show]->GetYaxis()->SetTitleOffset(1.5);
				//latex->DrawLatex(0.35, 0.8,"backfolded/measured");
				latex->DrawLatex(0.35, 0.8,Form("%s, %s",unfname[unf_show].Data(),test_name2[test].Data()));
				latex->DrawLatex(0.35, 0.7,Form("R=%.1lf, pTlead>%.0lf",R[r],pTthresh[ptl]));
				
				float max=htest_bf[r][ptl][unf_show]->GetMaximum();
				TLine *one = new TLine(xcutbf, 0, xcutbf, max);
				one->SetLineWidth(2);
				one->SetLineStyle(2);
				one->SetLineColor(kBlack);
				one->DrawClone("same");
		}//r
	}//pTlead	
	c_bf->SaveAs(Form("./obr/data/test_bf_%s_unf%i_bin%i.%s",test_name[test].Data(),unf_show,bin,grsuff.Data()));
	
	TCanvas* c_suc=new TCanvas("c_suc","successive_iterations",10,10,nr*can_x,nptl*can_y);
	c_suc->Divide(nr,nptl);
	gPad->SetLeftMargin(0.15);
	int counter=1;
	for(int ptl=ptl_start; ptl<=ptl_end; ptl++)
	{
		for(int r=0; r<nr; r++)
		{
			c_suc->cd(counter);
			counter++;
		
				TString hname=Form("htest_suc_R0%.0lf_pTl%.0lf_unf%i",R[r]*10,pTthresh[ptl],unf_show);
				htest_suc[r][ptl][unf_show]=(TH1D*) fin->Get(hname);
				//htest_suc2d[r][ptl][unf_show]->SetLineColor(hcolor[unf]);
				htest_suc[r][ptl][unf_show]->GetYaxis()->SetTitle("");
				htest_suc[r][ptl][unf_show]->GetYaxis()->SetTitleSize(0.04);
				htest_suc[r][ptl][unf_show]->GetYaxis()->SetTitleOffset(1.2);
				htest_suc[r][ptl][unf_show]->Draw("");
				latex->DrawLatex(0.35, 0.8,Form("%s, %s",unfname[unf_show].Data(),test_name2[test].Data()));
				latex->DrawLatex(0.35, 0.7,Form("R=%.1lf, pTlead>%.0lf",R[r],pTthresh[ptl]));
				
				float max=htest_suc[r][ptl][unf_show]->GetMaximum();
				TLine *one = new TLine(xcutsuc, 0, xcutsuc, max);
				one->SetLineWidth(2);
				one->SetLineStyle(2);
				one->SetLineColor(kBlack);
				one->DrawClone("same");
		}//r
	}//pTlead	
	c_suc->SaveAs(Form("./obr/data/test_suc_%s_unf%i_bin%i.%s",test_name[test].Data(),unf_show,bin,grsuff.Data()));
	
		TCanvas* c_cur=new TCanvas("c_cur","Curvature test",10,10,nr*can_x,nptl*can_y);
	c_cur->Divide(nr,nptl);
	int counter=1;
	for(int ptl=ptl_start; ptl<=ptl_end; ptl++)
	{
		for(int r=0; r<nr; r++)
		{
			c_cur->cd(counter);
			gPad->SetLeftMargin(0.13);
			counter++;

				TString hname=Form("htest_cur_R0%.0lf_pTl%.0lf_unf%i",R[r]*10,pTthresh[ptl],unf_show);
				htest_cur[r][ptl][unf_show]=(TH1D*) fin->Get(hname);
				//htest_cur2d[r][ptl][unf_show]->SetLineColor(hcolor[unf]);
				htest_cur[r][ptl][unf_show]->GetYaxis()->SetTitle("");
				htest_cur[r][ptl][unf_show]->GetYaxis()->SetTitleSize(0.04);
				htest_cur[r][ptl][unf_show]->GetYaxis()->SetTitleOffset(1.5);
				htest_cur[r][ptl][unf_show]->Draw("");
				latex->DrawLatex(0.35, 0.8,Form("%s, curvature",unfname[unf_show].Data()));
				latex->DrawLatex(0.35, 0.7,Form("R=%.1lf, pTlead>%.0lf",R[r],pTthresh[ptl]));
				
				float max=htest_cur[r][ptl][unf_show]->GetMaximum();
				TLine *one = new TLine(xcutcur, 0, xcutcur, max);
				one->SetLineWidth(2);
				one->SetLineStyle(2);
				one->SetLineColor(kBlack);
				one->DrawClone("same");
		}//r
	}//pTlead
	c_cur->SaveAs(Form("./obr/data/test_cur_%s_unf%i_bin%i.%s",test_name[test].Data(),unf_show,bin,grsuff.Data()));

	}
	
	
	/*
	TCanvas* c2=new TCanvas("c2","tests",10,10,nunf*can_x,3*can_y);
	c2->Divide(nunf,3);
	
	for(int unf=0; unf<nunf; unf++)
	{
		for(int pl=0; pl<nplots; pl++)
		{
				TString hname=Form("htest_bf_R0%.0lf_pTl%.0lf_unf%i",R[r_show]*10,pTthresh[pTl_show],unf);
				htest_bf[unf][pl]=(TH1D*) fin[pl]->Get(hname);
				htest_bf[unf][pl]->Rebin(2);
				hname=Form("htest_suc_R0%.0lf_pTl%.0lf_unf%i",R[r_show]*10,pTthresh[pTl_show],unf);
				htest_suc[unf][pl]=(TH1D*) fin[pl]->Get(hname);
				htest_suc[unf][pl]->Rebin(2);
				hname=Form("htest_cur_R0%.0lf_pTl%.0lf_unf%i",R[r_show]*10,pTthresh[pTl_show],unf);
				htest_cur[unf][pl]=(TH1D*) fin[pl]->Get(hname);
				
 
				c2->cd(1+unf);
				htest_bf[unf][pl]->SetLineColor(test_color[pl]);
				htest_bf[unf][pl]->Draw(drawcmd[pl]);
				if(unf==0)latex->DrawLatex(0.9, 0.95,Form("%s",test_name[test].Data()));
				latex->DrawLatex(0.45, 0.75,Form("%s",unfname[unf].Data()));
				//latex->DrawLatex(0.35, 0.7,"backfolded/measured");
				
				float xcut=backfold_cut[peripheral][test][unf];
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
				
				xcut=change_cut[peripheral][test][unf];
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
}
