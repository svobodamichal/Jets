#include "util.h"


//will take input histograms from .root file, normalize by the number of events and plot spectra for multiple pTlead cuts. Output: .root file, .pdf and .png images. 

using namespace std;


void ploteffi(string prod = "MC")
{
 
	TFile *f1 = TFile::Open(Form("pythia6_%s.root", prod.c_str()));
	TFile *fout1 = new TFile("pythia6_effi.root","recreate");

	int k=0;
	double R;
	//int pTlead = 0;
	array<double, 3> Rarr = {0.2, 0.3, 0.4};
	array<TString, 2> centrality = {"central","peripheral"};
	array<TString, 2> centbin = {"0-10","60-80"};
	
	TH1D *hMcpT_central[10][3];
	TH1D *hMcMatchedpT_central[10][3];
	TH1D *heffi_central[10][3];
	
	TH1D *hMcpT_peripheral[10][3];
	TH1D *hMcMatchedpT_peripheral[10][3];
	TH1D *heffi_peripheral[10][3];



	TCanvas* can = new TCanvas("can", "can", 1600, 1400);
	TLegend* leg = new TLegend(0.39, 0.34, 0.70, 0.69);
	leg->SetBorderSize(0);
  	leg->SetFillColor(0);


	TLatex* latex = new TLatex();
	latex->SetTextSize(0.05);
	latex->SetNDC();
	double lleft = 0.15; 
	double lbottom = 0.10;
	double lstep = 0.06;

	for(double &R : Rarr)
	{
		gPad->SetLogz();
		can->SetRightMargin(0.15);
		can->SetTopMargin(0.015);
		can->SetBottomMargin(0.11);
		can->SetLeftMargin(0.11);

		gStyle->SetOptStat(0);
		gStyle->SetTextFont(62);
		gStyle->SetTextSize(0.03);
		gStyle->SetLegendFont(62);
		gStyle->SetLegendTextSize(0.05);
		TGaxis::SetMaxDigits(3);
		
		
		//cout << "first" << endl;
		for(int pTlead =0; pTlead < 10; pTlead++)
		{
			hMcpT_central[pTlead][k]=(TH1D*)f1->Get(Form("hMcpT_pTl%i_R0%.0f_cent1", pTlead, R*10));
			hMcpT_central[pTlead][k]->Rebin(5);
			hMcMatchedpT_central[pTlead][k]=(TH1D*)f1->Get(Form("hMcmatchedpT_pTl%i_R0%.0f_cent1", pTlead, R*10));	
			hMcMatchedpT_central[pTlead][k]->Rebin(5);
			heffi_central[pTlead][k]=(TH1D*)hMcMatchedpT_central[pTlead][k]->Clone(Form("heffi_pTl%i_R0%.0f_central", pTlead, R*10));
			//heffi_central[pTlead][k]->Divide(hMcpT_central[pTlead][k]);	
			heffi_central[pTlead][k]->Divide(hMcMatchedpT_central[pTlead][k],hMcpT_central[pTlead][k],1,1,"B");	
					
			hMcpT_peripheral[pTlead][k]=(TH1D*)f1->Get(Form("hMcpT_pTl%i_R0%.0f_cent8", pTlead, R*10));
			hMcpT_peripheral[pTlead][k]->Rebin(5);
			hMcMatchedpT_peripheral[pTlead][k]=(TH1D*)f1->Get(Form("hMcmatchedpT_pTl%i_R0%.0f_cent8", pTlead, R*10));	
			hMcMatchedpT_peripheral[pTlead][k]->Rebin(5);
			heffi_peripheral[pTlead][k]=(TH1D*)hMcMatchedpT_peripheral[pTlead][k]->Clone(Form("heffi_pTl%i_R0%.0f_peripheral", pTlead, R*10));
			//heffi_peripheral[pTlead][k]->Divide(hMcpT_peripheral[pTlead][k]);			
			heffi_peripheral[pTlead][k]->Divide(hMcMatchedpT_peripheral[pTlead][k],hMcpT_peripheral[pTlead][k],1,1,"B");	
		
			heffi_central[pTlead][k]->SetTitle("");
			heffi_central[pTlead][k]->SetMarkerStyle(20);
			heffi_central[pTlead][k]->SetMarkerColor(kBlack);
			heffi_central[pTlead][k]->GetXaxis()->SetRangeUser(0.,65.);
			heffi_central[pTlead][k]->GetXaxis()->SetTitle("p^{MC}_{T,jet}");			
			heffi_central[pTlead][k]->GetYaxis()->SetTitle("matching efficiency");			
			heffi_central[pTlead][k]->Draw();
			//latex->DrawLatex(lleft + 0*lstep, lbottom+1*lstep, Form("PYTHIA6 p+p #otimes Au+Au %s%%",centbin[0].Data()));	
			latex->DrawLatex(lleft + 0*lstep, lbottom+1*lstep, "PYTHIA6 p+p #otimes Au+Au MB");				
			latex->DrawLatex(lleft + 0*lstep, lbottom+2*lstep, Form("Anti-k_{T}, R = %.1f, #it{p}_{T}^{lead} > %d GeV/#it{c}", R, pTlead));
	
			//can->SaveAs(Form("effi_R0%.0f_pTlead%i_central.png", R*10, pTlead));
			//can->Clear();

			heffi_peripheral[pTlead][k]->SetTitle("");
			heffi_peripheral[pTlead][k]->SetMarkerStyle(21);
			heffi_peripheral[pTlead][k]->SetMarkerColor(kBlue);
			heffi_peripheral[pTlead][k]->GetXaxis()->SetRangeUser(0.,65.);
			heffi_peripheral[pTlead][k]->Draw("same");
			//latex->DrawLatex(lleft + 0*lstep, lbottom+1*lstep, Form("PYTHIA6 p+p #otimes Au+Au %s%%",centbin[1].Data()));	
			//latex->DrawLatex(lleft + 0*lstep, lbottom+2*lstep, Form("Anti-k_{T}, R = %.1f, #it{p}_{T}^{lead} > %d GeV/#it{c}", R, pTlead));
			//can->SaveAs(Form("effi_R0%.0f_pTlead%i_peripheral.png", R*10, pTlead));
			
			if (pTlead ==0){leg->AddEntry(heffi_central[pTlead][k],"central 0-10%");
			leg->AddEntry(heffi_peripheral[pTlead][k],"peripheral 60-80%");	}		
			leg->Draw("same");
			
			can->SaveAs(Form("effi_R0%.0f_pTlead%i_cent.png", R*10, pTlead));			
			can->Clear();	
		}
		//can->Clear();
		leg->Clear();
		
		heffi_central[0][k]->Draw();
		heffi_central[5][k]->Draw("same"); 	
		heffi_central[5][k]->SetMarkerColor(kBlue);	
		heffi_central[7][k]->Draw("same"); 	
		heffi_central[7][k]->SetMarkerColor(kRed);		
		latex->DrawLatex(lleft + 0*lstep, lbottom+1*lstep, "PYTHIA6 p+p #otimes Au+Au 0-10%");				
		latex->DrawLatex(lleft + 0*lstep, lbottom+2*lstep, Form("Anti-k_{T}, R = %.1f", R));	
		leg->AddEntry(heffi_central[0][k],"#it{p}_{T}^{lead} > 0 GeV/#it{c}");	
		leg->AddEntry(heffi_central[5][k],"#it{p}_{T}^{lead} > 5 GeV/#it{c}");			
		leg->AddEntry(heffi_central[7][k],"#it{p}_{T}^{lead} > 7 GeV/#it{c}");	
		leg->Draw("same");			
		can->SaveAs(Form("effi_R0%.0f_central_pTlead.png", R*10));		
		can->Clear();
		
		heffi_peripheral[0][k]->Draw();
		heffi_peripheral[0][k]->SetMarkerColor(kBlack);	
		heffi_peripheral[5][k]->Draw("same"); 	
		heffi_peripheral[5][k]->SetMarkerColor(kBlue);	
		heffi_peripheral[7][k]->Draw("same"); 	
		heffi_peripheral[7][k]->SetMarkerColor(kRed);		
		latex->DrawLatex(lleft + 0*lstep, lbottom+1*lstep, "PYTHIA6 p+p #otimes Au+Au 60-80%");				
		latex->DrawLatex(lleft + 0*lstep, lbottom+2*lstep, Form("Anti-k_{T}, R = %.1f", R));	
		leg->Draw("same");			
		can->SaveAs(Form("effi_R0%.0f_peripheral_pTlead.png", R*10));	
		leg->Clear();
		k++;
	}
	fout1->cd();
	fout1->Write();	
	fout1->Close();	

}
