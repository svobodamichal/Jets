#include <iostream>
#include <cmath>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"

using namespace std;

void plotpriors_orig(){
	TFile *fin = TFile::Open("priors_default.root"/*,"update"*/);
	TCanvas *can = new TCanvas("can","",1600,1600);
	TF1 *ftsalis[9];
	TF1 *fpowlaw[3];
	can->SetTopMargin(0.01);
	can->SetRightMargin(0.03);	
	gPad->SetLogy();
	TLegend *leg = new TLegend(0.7,0.6,0.9,0.9);
	Color_t colors[8] = {kMagenta,kGreen+2,kCyan,kRed+2,kOrange,kViolet,kSpring,kYellow+2};
	
	Double_t *params;
	
	for (int i = 0; i<9; i++)
	{
		ftsalis[i]=(TF1*)fin->Get(Form("tsalis_%i",i+1));
		//cout << ftsalis[i]->GetNpar() << " " <<  ftsalis[i]->GetParName(1)<< endl;
		params = ftsalis[i]->GetParameters(); 
		ftsalis[i]->SetLineColor(i+2);
		if (!i) ftsalis[i]->Draw(""); 
		else ftsalis[i]->Draw("lsame");
	}
	
	ftsalis[0]->GetHistogram()->SetYTitle("a.u.");
	ftsalis[0]->GetHistogram()->GetYaxis()->SetRangeUser(1e-10,1);
	Double_t power[3] = {45.,5.,55.};
	
	for (int j = 0; j<3; j++)
	{
		fpowlaw[j]=(TF1*)fin->Get(Form("powlaw%.0f",power[j]));
		fpowlaw[j]->SetLineColor(colors[j]);
	}
	
	TFile *fin2 = TFile::Open("pythia6_all.root");	
	TH1D *hPythia6 = fin2->Get("hMcmatchedpT_pTl5_R02_cent1");
	hPythia6->Scale(1e5);
	hPythia6->SetMarkerStyle(20);
	hPythia6->Draw("same");
	
	can->SaveAs("priors.png");
	can->Clear();

	
	
}
