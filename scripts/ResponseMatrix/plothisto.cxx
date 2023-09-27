#include "util.h"


//will take input histograms from .root file, normalize by the number of events and plot spectra for multiple pTlead cuts. Output: .root file, .pdf and .png images. 

using namespace std;


void plothisto(string prod = "combined_response")
{
 
	TFile *f1 = TFile::Open(Form("Pythia6_%s.root", prod.c_str()));
	TFile *fout1 = new TFile(Form("pythia6_normalized_%s.root", prod.c_str()), "recreate");



	int k=0;
	TString c[10][3][9];
	//TString cf[9][2];
 	

	//TDirectoryFile* list1 = (TDirectoryFile *)f1->Get("TowerMCtask"); 
	TList* list1 = (TList *)f1->Get("stPicoHFJetMaker"); 	

	//TH1D* hevents = (TH1D*)list1->Get("hEventStat");
	TH1D* hevents = (TH1D*)list1->FindObject("hEventStat1");	
	//TH1D* hcent = (TH1D*)list1->FindObject("hcent");
	//double Nevents = hevents->GetBinContent(3);
	double Nevents = hevents->GetBinContent(1); //normalize by number of THROWN Pythia events 
	cout << "Nevents: " << Nevents << endl;
	double Neventsbin;

	double R;
	//int pTlead = 0;
	array<double, 3> Rarr = {0.2, 0.3, 0.4};
	array<TString, 9> centrality = {"","central", "", "", "", "","","","peripheral"};
	array<TString, 9> centbin = {"","0-10", "", "", "", "","","","60-80"};
	
	TH2D* hResponseM_tmp[10][3][9];
	TH2D* hResponseM[10][3][9];
	TH2D* hResponseMrebin[10][3][9]; //extend axis range

	TH2D *hptleads[3][9];
	TH2D *hptleads_tmp[3][9];

	TH1D *hMcpT[10][3][9];
	TH1D *hRcpT[10][3][9];	
	TH1D *hMcMatchedpT[10][3][9];
	TH1D *hRcMatchedpT[10][3][9];	
	TH1D *heffi[10][3][9];


	TCanvas* can = new TCanvas("can", "can", 1600, 1400);
	TLegend* leg = new TLegend(0.49, 0.54, 0.70, 0.89);



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
			for(int cent = 0; cent < 9; cent++)
			{
				//cout << "first" << endl;
				for(int pTlead =0; pTlead < 10; pTlead++)
				{

			//	hMcpT[pTlead][k][cent]=(TH1D*)list1->FindObject(Form("hMcpT_pTl%i_R0%.0f_cent%i", pTlead, R*10,cent));
				hMcpT[pTlead][k][cent]=(TH1D*)list1->FindObject(Form("hMCpT_pTl%i_R0%.0f_centbin%i", pTlead, R*10,cent));	
				hMcpT[pTlead][k][cent]->Scale(1./Nevents);
				
				//hRcpT[pTlead][k][cent]=(TH1D*)list1->FindObject(Form("hfpT_pTl%i_R0%.0f_centbin%i", pTlead, R*10,cent));	
				//hRcpT[pTlead][k][cent]->Scale(1./Nevents);
				
			//	hMcMatchedpT[pTlead][k][cent]=(TH1D*)list1->FindObject(Form("hMcmatchedpT_pTl%i_R0%.0f_cent%i", pTlead, R*10,cent));	
				hMcMatchedpT[pTlead][k][cent]=(TH1D*)list1->FindObject(Form("hMCmatchedpT_pTl%i_R0%.0f_centbin%i", pTlead, R*10,cent));		
				hMcMatchedpT[pTlead][k][cent]->Scale(1./Nevents);
//				heffi[pTlead][k][cent]=(TH1D*)hMcMatchedpT[pTlead][k][cent]->Clone(Form("heffi_pTl%i_R0%.0f_cent%i", pTlead, R*10,cent));
//				heffi[pTlead][k][cent]->Divide(hMcpT[pTlead][k][cent]);		

				hRcMatchedpT[pTlead][k][cent]=(TH1D*)list1->FindObject(Form("hRCmatchedpT_pTl%i_R0%.0f_centbin%i", pTlead, R*10,cent));
				hRcMatchedpT[pTlead][k][cent]->Scale(1./Nevents);
			
			

				c[pTlead][k][cent] = Form("hResponseMatrix_pTl%i_R0%.0f_centbin%i", pTlead, R*10,cent);
				//cf[j][k] = Form("hfjetpTlead_R0%.0f_centbin%i", R*10, centbin);
				//} else {c[j][k] = Form("hjetpT_R0%.0f_centbin%i_corr", R*10, centbin);cf[j][k] = Form("hfjetpT_R0%.0f_centbin%i_corr", R*10, centbin);}
		

				cout << "name: " << c[pTlead][k][cent] << endl;
				hResponseM_tmp[pTlead][k][cent] = (TH2D*)list1->FindObject(c[pTlead][k][cent]);	
				if (hResponseM_tmp[pTlead][k][cent]->GetEntries() == 0) {cout << "empty histogram: " << hResponseM_tmp[pTlead][k][cent]->GetName() << ", skipping! " << endl; continue;}
				hResponseM_tmp[pTlead][k][cent]->SetTitle("");
				//hResponseM_tmp[pTlead][k][cent]->RebinY(10);
			//hrescharged[j][k]->Rebin(5);
				hResponseM_tmp[pTlead][k][cent]->Scale(1./Nevents);
					
				//normalize each column in pTtrue to 1
				hResponseM[pTlead][k][cent] = (TH2D*)hResponseM_tmp[pTlead][k][cent]->Clone(Form("hResponseMatrix_pTl%i_R0%.0f_centbin%i", pTlead, R*10,cent));
				double integral = -1;
				/*for (int y = 1; y < hResponseM_tmp[pTlead][k][cent]->GetNbinsY()+1; y++) {
					integral = hResponseM_tmp[pTlead][k][cent]->Integral(hResponseM_tmp[pTlead][k][cent]->FindFirstBinAbove(0), hResponseM_tmp[pTlead][k][cent]->FindLastBinAbove(0), y, y);
					if (integral == 0) integral = 1;
						for (int x = 1; x < hResponseM_tmp[pTlead][k][cent]->GetNbinsX()+1; x++) {
						//if (hResponseM_tmp[pTlead][k][cent]->GetBinContent(x,y) > 0) cout << hResponseM_tmp[pTlead][k][cent]->GetBinError(x,y)/hResponseM_tmp[pTlead][k][cent]->GetBinContent(x,y) << endl;
						hResponseM[pTlead][k][cent]->SetBinContent(x,y, hResponseM_tmp[pTlead][k][cent]->GetBinContent(x,y)/integral);
						}							
					}*/	
		
				//hResponseM[pTlead][k]->GetZaxis()->SetRangeUser(1e-6,1);
				hResponseM[pTlead][k][cent]->SetXTitle("p_{T}^{det} [GeV/#it{c}]");

				if (!(!strcmp(centrality[cent].Data(),"peripheral") || !strcmp(centrality[cent].Data(),"central"))) continue;
			//if (centbin == 7|| centbin == 8) hrescharged[j][k]->Scale(2);
				if (hResponseM_tmp[pTlead][k][cent-1]->GetEntries()!=0) {
				cout << "adding "<< hResponseM[pTlead][k][cent]->GetName() << " and " << hResponseM[pTlead][k][cent-1]->GetName()<<endl;
				//adding
				hResponseM[pTlead][k][cent]->Add(hResponseM[pTlead][k][cent-1]);
				//hResponseM[pTlead][k][cent]->Write();
				
				//heffi[pTlead][k][cent]->Add(heffi[pTlead][k][cent-1]);
				//heffi[pTlead][k][cent]->Write();
				hMcpT[pTlead][k][cent]->Add(hMcpT[pTlead][k][cent-1]);
				hMcpT[pTlead][k][cent]->Write();
				hMcMatchedpT[pTlead][k][cent]->Add(hMcMatchedpT[pTlead][k][cent-1]);
				hMcMatchedpT[pTlead][k][cent]->Write();
				hRcMatchedpT[pTlead][k][cent]->Add(hRcMatchedpT[pTlead][k][cent-1]);
				hRcMatchedpT[pTlead][k][cent]->Write();
				}
				hResponseM[pTlead][k][cent]->Draw("colz");
				hResponseM[pTlead][k][cent]->Write();


		latex->DrawLatex(lleft + 0*lstep, lbottom+1*lstep, Form("PYTHIA6 p+p #otimes Au+Au %s%%",centbin[cent].Data()));	
		latex->DrawLatex(lleft + 0*lstep, lbottom+2*lstep, Form("Anti-k_{T}, R = %.1f, #it{p}_{T}^{lead} > %d GeV/#it{c}", R, pTlead));
		//latex->DrawLatex(lleft, lbottom+2*lstep, Form("Char. jets, anti-k_{T}, R = %.1f", R));				
	
						can->SaveAs(Form("ResponseMatrix_R0%.0f_pTlead%i_%s.png", R*10, pTlead, centbin[cent].Data()));
						can->SaveAs(Form("ResponseMatrix_R0%.0f_pTlead%i_%s.pdf", R*10, pTlead, centbin[cent].Data()));
						can->Clear();


					//extend bin range
						/*hResponseMrebin[pTlead][k][cent] = new TH2D(Form("hResponseMatrix_pTl%i_R0%.0lf_%s",pTlead,R*10,centrality[cent].Data()), Form("Response Matrix for p_{T}lead>%i ; p_{T}^{det} (GeV/c); p_{T}^{true} (GeV/c)",pTlead), 800, -100, 100, 800, -100, 100);
						hResponseMrebin[pTlead][k][cent]->Sumw2();
						for (int y = 1; y < hResponseM[pTlead][k][cent]->GetNbinsY()+1; y++) {
							for (int x = 1; x < hResponseM[pTlead][k][cent]->GetNbinsX()+1; x++) {
							cout << x << " " << y << " " <<  hResponseM[pTlead][k][cent]->GetXaxis()->GetBinCenter(x) << " " << hResponseMrebin[pTlead][k][cent]->GetXaxis()->GetBinCenter(x+440) << " " << hResponseM[pTlead][k][cent]->GetYaxis()->GetBinCenter(y) << " " << hResponseMrebin[pTlead][k][cent]->GetYaxis()->GetBinCenter(y+440) << endl;
							hResponseMrebin[pTlead][k][cent]->SetBinContent(x+440,y+440, hResponseM[pTlead][k][cent]->GetBinContent(x,y));
							hResponseMrebin[pTlead][k][cent]->SetBinError(x+440,y+440, hResponseM[pTlead][k][cent]->GetBinError(x,y));
							}							
						}	
						hResponseMrebin[pTlead][k][cent]->Write();
						delete hResponseMrebin[pTlead][k][cent];
*/
					}

			/*hptleads_tmp[k][cent] = (TH2D*)list1->FindObject(Form("hpTleads_R0%.0f_cent%i", R*10,cent));
			hptleads[k][cent] = (TH2D*)hptleads_tmp[k][cent]->Clone(Form("hpTleads_R0%.0f", R*10));
			double intgrl = -1;
						for (int y = 1; y < hptleads_tmp[k][cent]->GetNbinsY()+1; y++) {
							intgrl = hptleads_tmp[k][cent]->Integral(hptleads_tmp[k][cent]->FindFirstBinAbove(0), hptleads_tmp[k][cent]->FindLastBinAbove(0), y, y);
							if (intgrl == 0) intgrl = 1;
							for (int x = 1; x < hptleads_tmp[k][cent]->GetNbinsX()+1; x++) {
							hptleads[k][cent]->SetBinContent(x,y, hptleads_tmp[k][cent]->GetBinContent(x,y)/intgrl);
							}							
						}	
			hptleads[k][cent]->SetTitle("");
			//hptleads[k]->GetZaxis()->SetRangeUser(1e-8, 1e-4);
			hptleads[k][cent]->Draw("colz");

			latex->DrawLatex(lleft + 4*lstep, lbottom+6*lstep, Form("PYTHIA6 p+p #otimes Au+Au %s%%",centbin[cent].Data()));	
			latex->DrawLatex(lleft + 4*lstep, lbottom+7*lstep, Form("Anti-k_{T}, R = %.1f", R));
			
			can->SaveAs(Form("pTleads_R0%.0f.png", R*10));
			can->SaveAs(Form("pTleads_R0%.0f.pdf", R*10));
			can->Clear();
			*/
				}
				k++;
			}
		can->Clear();
		//leg->Clear();

	//fout1->Close();		

}
