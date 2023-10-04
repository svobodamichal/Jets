#include "util.h"


//will take input histograms from .root file and add together histos for central and peripheral collisions. Output: .root file. 

using namespace std;


//void plotspectra(string prod = "P18ih")
void preparehistosforunfolding(string prod = "trig")
{
 
//vpd mb 30+5
	/*TFile *f1 = TFile::Open(Form("result_%s_low.root", prod.c_str()));
	TFile *fout1 = new TFile(Form("spec_%s_low.root", prod.c_str()), "recreate");*/
//	TFile *f1 = TFile::Open(Form("out_%s_nocent.root", prod.c_str()));
//	TFile *fout1 = new TFile(Form("spec_%s_nocent.root", prod.c_str()), "recreate");
	//TFile *f1 = TFile::Open(Form("%s_lowmid.root", prod.c_str()));
	//TFile *fout1 = new TFile(Form("spec_%s_lowmid.root", prod.c_str()), "recreate");
	/*TFile *f1 = TFile::Open(Form("%s_nolumi.root", prod.c_str()));
	TFile *fout1 = new TFile(Form("spec_%s_nolumi.root", prod.c_str()), "recreate");*/
	TFile *f1 = TFile::Open(Form("out_HT2_%s.root", prod.c_str()));
	TFile *fout1 = new TFile(Form("histos_HT2_%s.root", prod.c_str()), "recreate");

	int j=0;
	int k=0;
	int l=0;

	TString cf[7][9];
 	

	TList* list1 = (TList *)f1->Get("stPicoHFJetMaker"); 


	double R;
	double deltaeta;
	const double twopi = 2*TMath::Pi();
	int pTlead = 0;
	array<int,7> pTarr = {0,3,4,5,6,7,9};
	array<double, 3> Rarr = {0.2, 0.3, 0.4};
	array<int, 9> centbins = {0,1,2,3,4,5,6,7,8};
	array<string, 9> centrality={"0-5%", "5-10%", "10-20%", "20-30%","30-40%", "40-50%", "50-60%", "60-70%", "70-80%"};
	//array<string, 9> centrality={"70-80 %", "60-70 %", "50-60 %", "40-50 %", "30-40 %", "20-30 %", "10-20 %", "5-10 %", "0-5 %"};
	array<string, 9> centrality1={"0-10%", "0 %", "0 %", "0 %","0 %", "0 %", "40-60%", "0 %", "60-80%"};


	TH1D* hfull[7][9];
	TH1D* hresfull[7][9];


	for(double &R : Rarr)
				{

					k=0;
				for(int &centbin : centbins)
				{
						j=0;	
					for(int &pTlead : pTarr)
					{
									//cout << "third" << endl;
					//if (pTlead > 0) {			
						cf[j][k] = Form("hfpT_pTl%i_R0%.0f_centbin%i", pTlead, R*10, centbin);
						//} else {c[j][k] = Form("hjetpT_R0%.0f_centbin%i_corr", R*10, centbin);cf[j][k] = Form("hfjetpT_R0%.0f_centbin%i_corr", R*10, centbin);}
		

	//cout << "name: " << c[j][k] << endl;
		
			hfull[j][k] = (TH1D*)list1->FindObject(cf[j][k]);	
			if (hfull[j][k]->GetEntries() == 0) {cout << "empty histogram: " << hfull[j][k]->GetName() << ", skipping! " << endl; continue;}
			hresfull[j][k] = (TH1D*)hfull[j][k]->Clone(hfull[j][k]->GetName());
			hresfull[j][k]->SetTitle("");
			hresfull[j][k]->SetXTitle("p_{T,jet}^{reco} [GeV/#it{c}]");
			hresfull[j][k]->Rebin(10);
			//if (centbin == 7|| centbin == 8) hresfull[j][k]->Scale(2);
			j++;
			}

		k++;
		}
	//cout << "adding" << endl;
		for (int i = 0; i < 7; i++){
		hresfull[i][8]->Add(hresfull[i][7]);

		hresfull[i][6]->Add(hresfull[i][5]);

		hresfull[i][0]->Add(hresfull[i][1]);

		}
		
		for (int i = 0; i<7; i++) {
		TString namecent = Form("hfpT_pTl%d_R0%.0f_peripheral", pTarr[i], R*10);
		hresfull[i][8]->Write(namecent);
		namecent = Form("hfpT_pTl%d_R0%.0f_central", pTarr[i], R*10);
		hresfull[i][0]->Write(namecent);
		
		}

		
		l++;

		//leg->Clear();	
	}	

	fout1->Close();		

}
