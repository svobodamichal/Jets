#include <iostream>
#include <cmath>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"

using namespace std;

void plotpriors(){
	TFile *fin = TFile::Open("priors_default.root"/*,"update"*/);
	TCanvas *can = new TCanvas("can","",1600,1600);
	//TF1 *ftsalis[9];
	can->SetTopMargin(0.01);
	can->SetRightMargin(0.03);	
	gPad->SetLogy();
	TLegend *leg = new TLegend(0.7,0.6,0.9,0.9);
	
	TFile *fout = new TFile("priors_mtsalis.root", "recreate");
	
	TF1* sigmoid = new TF1("sigmoid","1./(1+exp([0]*(-x+[1])))",0,100);
	TF1* ftxs = new TF1("txs","1e-9/(1+exp(-x+25))+x*(1+x/([0]*[1]))^-[1]",0,100);
	ftxs->SetParameters(1.4,10);
	double x0=20,x1 = 40;	
	TF1* ftxht = new TF1("txht",Form("(1+([2]*(1+tanh((x-%.1f)/[3]))/2-[2]*(1-tanh((x-%.1f)/[3]))/2))*x*(1+x/([0]*[1]))^-[1]",x0,x1),0,100);
	ftxht->SetParameters(0.9,12,0.1,1);
	ftxht->SetLineColor(kRed);
	ftxht->SetLineWidth(4);

	TH1* txs[9];
	sigmoid->SetParameters(1,25);
	
	Double_t *params;
	
	/*for (int i = 0; i<9; i++)
	{
		ftsalis[i]=(TF1*)fin->Get(Form("tsalis_%i",i+1));
		//cout << ftsalis[i]->GetNpar() << " " <<  ftsalis[i]->GetParName(1)<< endl;
		params = ftsalis[i]->GetParameters();
		//cout << "function: Tsalis_"<<i+1 << endl;
		//cout << "parameters: " << params[0]<< " " << params[1] << " " << params[2] << endl << endl;  
		ftsalis[i]->SetLineColor(i+1);
		if (!i) ftsalis[i]->Draw(""); 
		else ftsalis[i]->Draw("lsame");
		//txs[i] = (TH1D*)sigmoid->GetHistogram();
		//txs[i]->Multiply(ftsalis[i]);
		//cout << txs[i]->GetMinimum() << endl;
		//txs[i]->Draw("same");
	}
	ftxs->SetLineColor(1);
	cout << ftxs->Eval(25) << endl;
	ftxs->Draw("lsame");
	ftxht->Draw("lsame");
	ftsalis[0]->GetHistogram()->SetYTitle("a.u.");
	ftsalis[0]->GetHistogram()->GetYaxis()->SetRangeUser(1e-10,1);
	sigmoid->Draw("same");

	can->SaveAs("ftsalis.png");
	can->Clear();
	*/
	
	TF1* fpowlaw = new TF1("powlaw","1/(x^[0])",1,100);
	fpowlaw->SetParameter(0,5);
	//TF1 *fgauspol = new TF1("gauspol","gaus(0)*(1/(1+x^[3]))",0,100);
	//fgauspol->SetParameters(1,0,12,3);
	
		Color_t colors[8] = {kMagenta,kGreen+2,kCyan,kRed+2,kOrange,kViolet,kSpring,kYellow+2};
	
	TF1 *fexpol = new TF1("expol","expo(0)+(1/(1+x^[2]))",0,100); 
	fexpol->SetParameters(0,-0.5,5);
	TF1* fmtsalis[6];
	Double_t pararray[6][4] = {{0.6,12,0.1,1},{0.9,12,0.1,1},{1.2,12,0.1,1},{0.6,16,0.1,1},{0.9,16,0.1,1},{1.2,16,0.1,1}};
	for (int f = 0; f < 6; f++){
		fmtsalis[f]=(TF1*)ftxht->Clone(Form("mtsalis_%i",f+1));
		fmtsalis[f]->SetParameters(pararray[f]);
		fmtsalis[f]->SetLineColor(colors[f]+2);
		fmtsalis[f]->SetLineWidth(3);
		fmtsalis[f]->SetNpx(10000);
		fmtsalis[f]->GetYaxis()->SetRangeUser(1e-10,10);
		if (!f) fmtsalis[f]->Draw(""); 
		else fmtsalis[f]->Draw("lsame");
		fmtsalis[f]->Write();
		leg->AddEntry(fmtsalis[f],Form("mtsalis_%i",f+1));
		//fin->cd();
		fmtsalis[f]->Write();	
	}
	//fpowlaw->Draw("same");
	//fgauspol->Draw("same");
	//fgauspol->Write();
	fexpol->SetLineColor(kMagenta+1);
	//fexpol->Draw("");
	fexpol->Write();
	//leg->AddEntry(fexpol,"expol");
	//cout << fexpol->Eval(0) << endl;
	
	TF1* fmpowlaw = new TF1("mpowlaw_template", "((x < 5) ? gaus(1) : [4]/(1+x^[0]))",0,100);
	TF1* mpowlaw[5];
	TString pname[] = {"powlaw4","powlaw45","powlaw5","powlaw55","powlaw6"};
	Double_t powers[] = {4,4.5,5,5.5,6};
	Double_t amplitude[] = {5e2,1e3,2e3,5e3,1e4};	
	x0=8;
	TF1* gammapolorig = new TF1("gammapol", Form("TMath::GammaDist(x,[0],[1],[2])*(1-tanh((x-%.1f)/[4]))/2+1/(1+x^[3])*(1+tanh((x-%.1f)/[4]))/2",x0,x0),0,100); //
	TF1* origamma = new TF1("origamma", "TMath::GammaDist(x,[0],[1],[2])",0,100);
	TF1* gamma[5];
	Double_t pararraygamma[5][3] = {{5,0,2.5},{5,0,2},{5,0,3},{4,0,2.5},{6,0,2.5}};
	TF1* gammapol[5];
	Double_t pargammapol[5][5] = {{5,0,2.5,4,1},{8,0,2.5,4.5,1},{5,0,2.5,5,1},{5,0,1.5,4,1},{5,0,2.5,3.5,1}};
	
	
	for (int f = 0;f<5;f++){
		mpowlaw[f] = (TF1*)fmpowlaw->Clone(Form("m%s",pname[f].Data()));
		mpowlaw[f]->SetParameters(powers[f],1,6,1.5,amplitude[f]);
		//cout << "{" << mpowlaw[f]->GetParameter(0) <<"," << mpowlaw[f]->GetParameter(1)<< ","<< mpowlaw[f]->GetParameter(2)<< ","<< mpowlaw[f]->GetParameter(3) << "," <<mpowlaw[f]->GetParameter(4) << "}" << endl;
		mpowlaw[f]->SetLineColor(colors[f]);
		mpowlaw[f]->SetNpx(10000);
		mpowlaw[f]->GetYaxis()->SetRangeUser(1e-10,10);
		/*if (!f) mpowlaw[f]->Draw(""); else*/ mpowlaw[f]->Draw("lsame");
		mpowlaw[f]->Write();
		leg->AddEntry(mpowlaw[f],mpowlaw[f]->GetName());
		
		gamma[f] = (TF1*)origamma->Clone(Form("gamma_%i",f+1));
		gamma[f]->SetParameters(pararraygamma[f]);
		gamma[f]->SetLineColor(f+1);
		gamma[f]->SetNpx(10000);	
		//gamma[f]->Draw("lsame");
		//leg->AddEntry(gamma[f],Form("gamma_%i",f+1));
		//gamma[f]->Write();
		
		gammapol[f]=(TF1*)gammapolorig->Clone(Form("gammapol_%i",f+1));
		gammapol[f]->SetParameters(pargammapol[f]);
		gammapol[f]->SetLineColor(colors[f]);
		gammapol[f]->SetNpx(10000);
		//gammapol[f]->Draw("lsame");
		gammapol[f]->Write();
		//leg->AddEntry(gammapol[f],Form("gammapol_%i",f+1));
		
	}
	gammapolorig->SetParameters(5,0,2.5,4,1);
	//cout << gammapol->Eval(10) << " " << gammapol->Eval(2) << endl;
	gammapolorig->SetLineColor(kMagenta+1);
	origamma->SetParameters(5,0,2.5);
	origamma->SetLineColor(kGreen+1);
	//origamma->Draw("lsame");
	//gammapol->Draw("lsame");
	
	x0=15;	
	TF1* gausxpol = new TF1("gausxpol",Form("gaus(3)*[2]*(1-tanh((x-%.1f)/[1]))/2 + (200/(15+x^[0]))*[2]*(1+tanh((x-%.1f)/[1]))/2",x0,x0),0,100);
	gausxpol->SetParameters(4.5,2,1,1,10,1);
	gausxpol->SetLineColor(kMagenta+1);
	gausxpol->SetNpx(10000);
//	gausxpol->Draw("lsame");

	TF1 *gauspol[5];
	Double_t pargauspol[5][6] = {{4.5,2,1,1,10,2},{5.5,2,1,1,10,2},{5.,2,1,1,10,2},{5,2,1,1,9,3},{5,2,1,1,8,2}};
	for (int f = 0;f<5;f++){
		gauspol[f] = (TF1*)gausxpol->Clone(Form("gauspol_%i",f+1));
		gauspol[f]->SetParameters(pargauspol[f]);
		gauspol[f]->SetLineColor(f+1);
		gauspol[f]->Draw("lsame");
		gauspol[f]->Write();
		leg->AddEntry(gauspol[f],gauspol[f]->GetName());
	}
	
	leg->Draw("same");

	//can->SaveAs("mfunctions_new.png");
	can->SaveAs("mfunctions_modified.png");	
	
	
}
