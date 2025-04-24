#include "util.h"

using namespace std;
//will plot ratios of two consecutive unfolding iterations
void unfoldedratios()
{
	TString utype = gSystem->Getenv("UTYPE");
	TString centrality = gSystem->Getenv("CENTRALITY");
	TString cent = "60-80 %";
	if (centrality == "central") cent = "0-10 %";
	TString prior = gSystem->Getenv("PRIOR");
	double R = atof(gSystem->Getenv("R")), pTlead  = atof(gSystem->Getenv("PTLEAD"));
	TString infile_path = Form("./%s/%s/%s/unfolded_SignalSpectrum_R%.1f_pTthresh%.1f.root", utype.Data(), centrality.Data(), prior.Data(), R, pTlead);
	TFile *fin = TFile::Open(infile_path.Data());
	if (fin->IsOpen()) cout << "file " << infile_path.Data() << " is open" << endl;
	TDirectoryFile *inputdir = (TDirectoryFile*)fin->Get("input");

	TFile *fout = new TFile("spectra_unfolded.root","update");

	TCanvas* can = new TCanvas("can", "can", 1600, 1600);	
	TLegend* leg = new TLegend(0.78, 0.65, 0.94, 0.95);

	TLatex* latex = new TLatex();
	latex->SetTextSize(0.035);
	latex->SetTextFont(42);
	latex->SetNDC();
	double lleft = 0.14; 
	double lbottom = 0.52;
	double lstep = 0.06;

	//Color_t colors[8] = {kBlack,kRed,kBlue,kYellow,kCyan,kOrange,kMagenta,kGreen};
	Color_t colors[8] = {kMagenta,kGreen,kYellow,kCyan,kRed,kOrange,kBlue,kBlack};


	gPad->SetLogz();
	gStyle->SetOptStat(0);
	can->SetRightMargin(0.13);
	TH2D *hresponseM = (TH2D*)inputdir->Get("hresponse");

	//normalize each row to 1
	TH2D *hresponse = (TH2D*)hresponseM->Clone("hresponse");
	double integral = -1;
		for (int y = 1; y < hresponseM->GetNbinsY()+1; y++) {
		integral = hresponseM->Integral(hresponseM->FindFirstBinAbove(0), hresponseM->FindLastBinAbove(0), y, y);
		if (integral == 0) integral = 1;
		for (int x = 1; x < hresponseM->GetNbinsX()+1; x++) {
		hresponse->SetBinContent(x,y, hresponseM->GetBinContent(x,y)/integral);
				}							
		}
	
	//hresponse->GetZaxis()->SetRangeUser(1e-5,1e3);
	hresponse->GetZaxis()->SetLabelFont(42);
	hresponse->GetXaxis()->SetTitleFont(42);
	hresponse->GetXaxis()->SetLabelFont(42);
	hresponse->GetYaxis()->SetTitleFont(42);
	hresponse->GetYaxis()->SetLabelFont(42);
	hresponse->GetXaxis()->SetTitleSize(0.04);
	hresponse->GetXaxis()->SetLabelSize(0.04);
	hresponse->GetYaxis()->SetTitleSize(0.04);
	hresponse->GetYaxis()->SetLabelSize(0.04);
	hresponse->GetYaxis()->SetTitleOffset(1.15);
	hresponse->GetYaxis()->SetTitle("#it{p}_{T}^{true} [GeV/#it{c}]");
	hresponse->GetXaxis()->SetTitle("#it{p}_{T}^{meas} [GeV/#it{c}]");
	hresponse->SetTitle("");
	hresponse->Draw("colz0");

// Set x-axis to start from 0
    hresponse->GetXaxis()->SetRangeUser(0, hresponse->GetXaxis()->GetXmax());

// Move TLatex text to the right edge of the plot
    double lright = 0.85; // Adjust this value as needed
    latex->SetTextSize(0.025); // Decrease text size
    latex->SetTextAlign(31); // Align text to the right
    latex->DrawLatex(lright, lbottom-2.5*lstep, "#it{Response matrix}");
    latex->DrawLatex(lright, lbottom-3.5*lstep, Form("anti-#it{k}_{T}, #it{R} = %.1f", R));
    latex->DrawLatex(lright, lbottom-4.5*lstep, Form("#it{p}_{T,min}^{lead} = %.0f GeV/#it{c}", pTlead));
    latex->DrawLatex(lright, lbottom-5.5*lstep, "#sqrt{#it{s}_{NN}} = 200 GeV");
    latex->DrawLatex(lright, lbottom-6.5*lstep, Form("PYTHIA6 #it{p}+#it{p} #otimes STAR Au+Au %s", cent.Data()));

	can->SaveAs(Form("responseM_R0%.0f_%s_%s_pTlead%.0f.png",R*10,centrality.Data(), prior.Data(), pTlead));
	can->SaveAs(Form("responseM_R0%.0f_%s_%s_pTlead%.0f.pdf",R*10,centrality.Data(), prior.Data(), pTlead));
	can->Clear();

	gPad->SetLogy();
	TH1D *hprior = (TH1D*)inputdir->Get("hprior");
	hprior->Draw("");
	can->SaveAs(Form("prior_R0%.0f_%s_pTlead%.0f.png",R*10,prior.Data(), pTlead));
	can->SaveAs(Form("prior_R0%.0f_%s_pTlead%.0f.pdf",R*10,prior.Data(), pTlead));
	can->Clear();

	TH1D *hmeasured = (TH1D*)inputdir->Get("hmeasured");
	hmeasured->Draw("");
	can->SaveAs(Form("measured_R0%.0f_%s_pTlead%.0f.png",R*10,centrality.Data(), pTlead));
	can->SaveAs(Form("measured_R0%.0f_%s_pTlead%.0f.pdf",R*10,centrality.Data(), pTlead));
	can->Clear();

	TCanvas* can1 = new TCanvas("can1", "can1", 1600, 1600);
	can1->Divide(1,2,1e-5,1e-5);

	can1->cd();
	can1->SetRightMargin(0.05);
	gStyle->SetOptStat(0);
	gPad->SetLogy();
	gPad->SetTopMargin(0.02);
	gPad->SetBottomMargin(0.30);
	TDirectoryFile *iterationdir[8];
	TH1D *hunfolded[8]; 
	TH1D *hratio[7];
	int j=0; //will increment directory no 
	for (int i = 0; i < 8; i++)
	{
	 //iterationdir[i] = (TDirectoryFile*)fin->Get(Form("iter%d", i));
	 iterationdir[i] = (TDirectoryFile*)fin->Get(Form("iter%d", j));
	 hunfolded[i] = (TH1D*)iterationdir[i]->Get("hunfolded");
	 hunfolded[i]->Scale(1.0,"width");
	 hunfolded[i]->SetYTitle("dN/d#it{p}_{T} [a. u.]");
 	 hunfolded[i]->GetYaxis()->SetTitleOffset(1.2);
 	 hunfolded[i]->GetYaxis()->SetLabelFont(42);
 	 hunfolded[i]->GetYaxis()->SetTitleFont(42);
	 hunfolded[i]->SetMarkerColor(colors[i]);
	 hunfolded[i]->SetMarkerStyle(27-i);
	 hunfolded[i]->SetMarkerSize(2);
	 hunfolded[i]->SetTitle("");
	 hunfolded[i]->GetXaxis()->SetRangeUser(0,60);

	 //leg->AddEntry(hunfolded[i], Form("iteration: %d",i+1));
	 leg->AddEntry(hunfolded[i], Form("iteration: %d",j+1));
	//j+=5;	//increment by x each step
		j++;	//increment by 1 each step

	 if (!i) continue;
	 hratio[i-1]=(TH1D*)hunfolded[i-1]->Clone(Form("hratio%d",i));
	 hratio[i-1]->Divide(hunfolded[i]);
	 hratio[i-1]->GetXaxis()->SetLabelSize(0.13);
	 hratio[i-1]->SetYTitle("ratio iter i/i+1");
	 //hratio[i-1]->SetYTitle("ratio iter i/i+2");
	 hratio[i-1]->SetXTitle("#it{p}_{T,jet} [GeV/#it{c}]");
	 hratio[i-1]->GetYaxis()->SetTitleSize(0.13);
	 hratio[i-1]->GetXaxis()->SetTitleSize(0.13);
	 hratio[i-1]->GetXaxis()->SetTitleFont(42);
	 hratio[i-1]->GetXaxis()->SetLabelFont(42);
	 hratio[i-1]->GetYaxis()->SetLabelSize(0.1);
 	 hratio[i-1]->GetYaxis()->SetTitleOffset(0.35);
 	 hratio[i-1]->GetXaxis()->SetTitleOffset(1);
	 hratio[i-1]->GetYaxis()->SetRangeUser(0.8,1.2);
	}	

	hunfolded[0]->Draw("ep");
	for (int i = 1; i < 8; i++)
	{
		hunfolded[i]->Draw("ep same");
	}

	leg->SetBorderSize(0);
  leg->SetFillColor(0);
	leg->Draw("same");


	latex->DrawLatex(lleft+5*lstep, lbottom+7*lstep, "#it{Raw data}");
	latex->DrawLatex(lleft+5*lstep, lbottom+6*lstep, Form("STAR Au+Au %s", cent.Data()));
	latex->DrawLatex(lleft+6*lstep, lbottom+5*lstep, "#sqrt{#it{s}_{NN}} = 200 GeV");
	latex->DrawLatex(lleft, lbottom-0*lstep, Form("anti-#it{k}_{T}, #it{R} = %.1f", R));
	latex->DrawLatex(lleft, lbottom-lstep, Form("#it{p}_{T,min}^{lead} = %.0f GeV/#it{c}", pTlead));
	latex->DrawLatex(lleft, lbottom-2*lstep, Form("unfolding: %s", utype.Data()));
	latex->DrawLatex(lleft, lbottom-3*lstep, Form("prior: %s ", prior.Data()));
	
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
	pad2->SetTopMargin(0);
	pad2->SetLogy(0);
	pad2->SetBottomMargin(0.30);
	pad2->SetRightMargin(0.05);
  // pad2[cent]->SetGridx(); // vertical grid
	pad2->Draw();
	pad2->cd();       // pad2 becomes the current pad
	hratio[0]->Draw("ep");
	for (int i = 1; i < 7; i++)
	{
		hratio[i]->Draw("ep same");
	}

	TLine *line = new TLine();
	line->DrawLine(0,1,60,1);

	can1->SaveAs(Form("unfolded_%s_R0%.0f_%s_pTlead%.0f_%s.png",utype.Data(),R*10,centrality.Data(),pTlead,prior.Data()));
	can1->SaveAs(Form("unfolded_%s_R0%.0f_%s_pTlead%.0f_%s.pdf",utype.Data(),R*10,centrality.Data(),pTlead,prior.Data()));
	can1->Clear();
	
	TH1D* hbf[8];
	TH1D* hbm[8];
	for (int i = 0;i<8;i++){
		hbf[i]=(TH1D*)iterationdir[i]->Get("hbackfolded");
		hbf[i]->Scale(1.0,"width");
		hbf[i]->SetMarkerStyle(21);
		hbf[i]->SetMarkerColor(kRed);
		hbf[i]->SetTitle("");
		hbm[i]=(TH1D*)iterationdir[i]->Get("hbfmratio");
		hbm[i]->SetMarkerStyle(21);
		hbm[i]->SetTitle("");
	}
	
	
	hunfolded[7]->GetXaxis()->SetRangeUser(-20,60);
	hunfolded[7]->Draw("ep");
	hbf[7]->Draw("epsame");
	hmeasured->Scale(1.0,"width");
	cout << hbf[7]->Integral(0,-1) << " hmeasured " << hmeasured->Integral(0,-1) << " unfolded " << hunfolded[7]->Integral(0,-1) << endl;
	cout << hbf[7]->Integral(0,-1,"width") << " hmeasured " << hmeasured->Integral(0,-1,"width") << " unfolded " << hunfolded[7]->Integral(0,-1,"width") << endl;	
	for (int x = 0; x<=hbf[7]->GetNbinsX();x++)
	{
		cout << hunfolded[7]->GetBinContent(x)*hbf[7]->GetBinWidth(x) << endl;
	}
	hmeasured->SetMarkerStyle(21);
	hmeasured->SetMarkerColor(kBlue);
	hmeasured->Draw("epsame");
	leg->Clear();
	leg->AddEntry(hunfolded[7],"unfolded");
	leg->AddEntry(hbf[7],"backfolded");
	leg->AddEntry(hmeasured,"measured");	
	leg->Draw("same");
	
	latex->DrawLatex(lleft+5*lstep, lbottom+7*lstep, "#it{Raw data}");
	latex->DrawLatex(lleft+5*lstep, lbottom+6*lstep, Form("STAR Au+Au %s", cent.Data()));
	latex->DrawLatex(lleft+6*lstep, lbottom+5*lstep, "#sqrt{#it{s}_{NN}} = 200 GeV");
	latex->DrawLatex(lleft, lbottom-0*lstep, Form("anti-#it{k}_{T}, #it{R} = %.1f", R));
	latex->DrawLatex(lleft, lbottom-lstep, Form("#it{p}_{T,min}^{lead} = %.0f GeV/#it{c}", pTlead));
	latex->DrawLatex(lleft, lbottom-2*lstep, Form("unfolding: %s", utype.Data()));
	latex->DrawLatex(lleft, lbottom-3*lstep, Form("prior: %s ", prior.Data()));
	//can1->Update();	
	//pad2->Draw();
	pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
	pad2->SetTopMargin(0);
	pad2->SetLogy(0);
	pad2->SetBottomMargin(0.30);
	pad2->SetRightMargin(0.05);
  // pad2[cent]->SetGridx(); // vertical grid
	pad2->Draw();
	pad2->cd();       // pad2 becomes the current pad
	
	hbm[7]->GetYaxis()->SetRangeUser(0.3,3);
    hbm[7]->GetXaxis()->SetRangeUser(0,60);

    hbm[7]->GetXaxis()->SetLabelSize(0.13);
	 hbm[7]->SetYTitle("ratio bf/meas");
	 //hbm[7]->SetYTitle("ratio iter i/i+2");
	 hbm[7]->SetXTitle("#it{p}_{T,jet} [GeV/#it{c}]");
	 hbm[7]->GetYaxis()->SetTitleSize(0.13);
	 hbm[7]->GetXaxis()->SetTitleSize(0.13);
	 hbm[7]->GetXaxis()->SetTitleFont(42);
	 hbm[7]->GetXaxis()->SetLabelFont(42);
	 hbm[7]->GetYaxis()->SetLabelSize(0.1);
 	 hbm[7]->GetYaxis()->SetTitleOffset(0.35);
 	 hbm[7]->GetXaxis()->SetTitleOffset(1);
	hbm[7]->Draw("ep");
	line->DrawLine(0,1,60,1);
	
	can1->SaveAs(Form("backfolded_%s_R0%.0f_%s_pTlead%.0f_%s.png",utype.Data(),R*10,centrality.Data(),pTlead,prior.Data()));
	can1->SaveAs(Form("backfolded_%s_R0%.0f_%s_pTlead%.0f_%s.pdf",utype.Data(),R*10,centrality.Data(),pTlead,prior.Data()));	

	fout->cd();
	hunfolded[3]->Write(Form("unfolded_%s_R0%.0f_%s_pTlead%.0f_%s_it4",utype.Data(),R*10,centrality.Data(),pTlead,prior.Data()));
    hunfolded[4]->Write(Form("unfolded_%s_R0%.0f_%s_pTlead%.0f_%s_it5",utype.Data(),R*10,centrality.Data(),pTlead,prior.Data()));
    hunfolded[5]->Write(Form("unfolded_%s_R0%.0f_%s_pTlead%.0f_%s_it6",utype.Data(),R*10,centrality.Data(),pTlead,prior.Data()));
	hunfolded[7]->Write(Form("unfolded_%s_R0%.0f_%s_pTlead%.0f_%s_it8",utype.Data(),R*10,centrality.Data(),pTlead,prior.Data()));	


}
