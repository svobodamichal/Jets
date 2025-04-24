#include "util.h"
#include "utils.C"
#include "binning.h"



//will take input histograms from .root file, normalize by the number of events and plot spectra for multiple pTlead cuts. Output: .root file, .pdf and .png images. 

using namespace std;


void plothisto(string prod = "combined_response")
{
 
	TFile *f1 = TFile::Open(Form("Pythia6_%s.root", prod.c_str()));
	TFile *fout1 = new TFile(Form("pythia6_normalized_%s.root", prod.c_str()), "recreate");



	int k=0;
	TString c[10][3][7];
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
    array<TString, 8> centrality = {"","central", "", "", "", "","","peripheral"};
    array<TString, 8> centbin = {"","0-10", "", "", "", "","","60-80"};

	TH2D* hResponseM_tmp[10][3][7];
	TH2D* hResponseM[10][3][7];
	TH2D* hResponseMrebin[10][3][7]; //extend axis range
    TH2D* hDeltaR[3][7];
    TH2D* hDeltaRw[3][7];
    TH2D* hNconst[3][7];
    TH2D* hNconstw[3][7];
    TH2D* hMCNconst[3][7];
    TH2D* hMCNconstw[3][7];
    TH2D* hEtaPhiMc_Rc[3];
    TH2D* hEtaPhiMc_Rcw[3];

	TH2D *hptleads[3][7];
	TH2D *hptleads_tmp[3][7];

	TH1D *hMcpT[10][3][7];
	TH1D *hRcpT[10][3][7];
	TH1D *hMcMatchedpT[10][3][7];
	TH1D *hRcMatchedpT[10][3][7];
	TH1D *heffi[10][3][7];
    TH1D *hEtaMc_Rc[3];
    TH1D *hPhiMc_Rc[3];
    TH1D *hEtaMc_Rcw[3];
    TH1D *hPhiMc_Rcw[3];
    TH1D *hMcMatchedMCpT[10][3][7];

    TH1D *hET_tow;


    TCanvas* can = new TCanvas("can", "can", 1600, 1400);
	TLegend* leg = new TLegend(0.49, 0.54, 0.70, 0.89);



	TLatex* latex = new TLatex();
	latex->SetTextSize(0.05);
	latex->SetNDC();
	double lleft = 0.15; 
	double lbottom = 0.10;
	double lstep = 0.06;

    hET_tow=(TH1D*)list1->FindObject(Form("hET_tow"));
    hET_tow->Scale(1./Nevents);
    hET_tow->Write();


    //VARIABLE BINNING - use arrays defined in "binning.h"
    Int_t binChoice = 6; //choose binning
    float pTrange=100; //pTrange in input histograms (has to be same in all histos!)
    int tmp_nbins;
    int tmp_nbins2;
    float* tmp_array=NULL;
    float* tmp_array2=NULL;
    if(binChoice==0)
    {
        tmp_nbins=nbinsarr0;
        tmp_nbins2=nbinsarr0;
        tmp_array=binarr0;
        tmp_array2=binarr0;
    }
    else if(binChoice==1)
    {
        tmp_nbins=nbinsarr1a;
        tmp_nbins2=nbinsarr1b;
        tmp_array=binarr1a;
        tmp_array2=binarr1b;
    }
    else if(binChoice==2)
    {
        tmp_nbins=nbinsarr2a;
        tmp_nbins2=nbinsarr2b;
        tmp_array=binarr2a;
        tmp_array2=binarr2b;
    }
    else if(binChoice==3)
    {
        tmp_nbins=nbinsarr3a;
        tmp_nbins2=nbinsarr3b;
        tmp_array=binarr3a;
        tmp_array2=binarr3b;
    }

    else if(binChoice==4)
    {
        tmp_nbins=nbinsarr4a;
        tmp_nbins2=nbinsarr4b;
        tmp_array=binarr4a;
        tmp_array2=binarr4b;
    }

    else if(binChoice==5)
    {
        tmp_nbins=nbinsarr5a;
        tmp_nbins2=nbinsarr5b;
        tmp_array=binarr5a;
        tmp_array2=binarr5b;
    }

    else if(binChoice==6)
    {
        tmp_nbins=nbinsarr6a;
        tmp_nbins2=nbinsarr6b;
        tmp_array=binarr6a;
        tmp_array2=binarr6b;
    }

    else if(binChoice==7)
    {
        tmp_nbins=nbinsarr7a;
        tmp_nbins2=nbinsarr7b;
        tmp_array=binarr7a;
        tmp_array2=binarr7b;
    }

    else if(binChoice==8)
    {
        tmp_nbins=nbinsarr8a;
        tmp_nbins2=nbinsarr8b;
        tmp_array=binarr8a;
        tmp_array2=binarr8b;
    }

    else if(binChoice==9)
    {
        tmp_nbins=nbinsarr9a;
        tmp_nbins2=nbinsarr9b;
        tmp_array=binarr9a;
        tmp_array2=binarr9b;
    }

    const Int_t newbins=tmp_nbins;
    float *pTbinArray=tmp_array;

    const Int_t newbins2=tmp_nbins2;
    float *pTbinArray2=tmp_array2;

    cout<<"binning choice:"<<binChoice<<endl;
    for(int i=0; i<=newbins; i++){
        cout<<pTbinArray[i]<<",";
    }cout<<endl;
    for(int i=0; i<=newbins2; i++){
        cout<<pTbinArray2[i]<<",";
    }
    cout<<endl;





	for(double &R : Rarr)
	{

        hEtaPhiMc_Rc[k]=(TH2D*)list1->FindObject(Form("hEtaPhi_MC-RC_R0%.0f", R*10));
        hEtaPhiMc_Rc[k]->Scale(1./Nevents);

        hEtaMc_Rc[k]=(TH1D*)list1->FindObject(Form("heta_MCRC_R0%.0f", R*10));
        hEtaMc_Rc[k]->Scale(1./Nevents);

        hPhiMc_Rc[k]=(TH1D*)list1->FindObject(Form("hphi_MCRC_R0%.0f", R*10));
        hPhiMc_Rc[k]->Scale(1./Nevents);

        hEtaPhiMc_Rc[k]->Write();

        hEtaMc_Rc[k]->Write();

        hPhiMc_Rc[k]->Write();


        hEtaPhiMc_Rcw[k]=(TH2D*)list1->FindObject(Form("hEtaPhi_MC-RCw_R0%.0f", R*10));
        hEtaPhiMc_Rcw[k]->Scale(1./Nevents);

        hEtaMc_Rcw[k]=(TH1D*)list1->FindObject(Form("heta_MCRCw_R0%.0f", R*10));
        hEtaMc_Rcw[k]->Scale(1./Nevents);

        hPhiMc_Rcw[k]=(TH1D*)list1->FindObject(Form("hphi_MCRCw_R0%.0f", R*10));
        hPhiMc_Rcw[k]->Scale(1./Nevents);

        hEtaPhiMc_Rcw[k]->Write();

        hEtaMc_Rcw[k]->Write();

        hPhiMc_Rcw[k]->Write();


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
			for(int cent = 1; cent < 8; cent++)
			{

                /*   hDeltaR[k][cent]=(TH2D*)list1->FindObject(Form("hDeltaR_R0%.0f_centbin%i", R*10,cent));
                   hDeltaR[k][cent]->Scale(1./Nevents);

                   hDeltaRw[k][cent]=(TH2D*)list1->FindObject(Form("hDeltaRw_R0%.0f_centbin%i", R*10,cent));
                   hDeltaRw[k][cent]->Scale(1./Nevents);

                   hNconst[k][cent]=(TH2D*)list1->FindObject(Form("hNconst_R0%.0f_centbin%i", R*10,cent));
                   hNconst[k][cent]->Scale(1./Nevents);

                   hNconstw[k][cent]=(TH2D*)list1->FindObject(Form("hNconstw_R0%.0f_centbin%i", R*10,cent));
                   hNconstw[k][cent]->Scale(1./Nevents);

                   hMCNconst[k][cent]=(TH2D*)list1->FindObject(Form("hMCNconst_R0%.0f_centbin%i", R*10,cent));
                   hMCNconst[k][cent]->Scale(1./Nevents);

                   hMCNconstw[k][cent]=(TH2D*)list1->FindObject(Form("hMCNconstw_R0%.0f_centbin%i", R*10,cent));
                   hMCNconstw[k][cent]->Scale(1./Nevents);

                   hDeltaR[k][cent]->Write();

                   hDeltaRw[k][cent]->Write();

                   hNconst[k][cent]->Write();

                   hNconstw[k][cent]->Write();

                   hMCNconst[k][cent]->Write();

                   hMCNconstw[k][cent]->Write();*/

				//cout << "first" << endl;
				for(int pTlead =0; pTlead < 10; pTlead++) {


                    //	hMcpT[pTlead][k][cent]=(TH1D*)list1->FindObject(Form("hMcpT_pTl%i_R0%.0f_cent%i", pTlead, R*10,cent));
                    hMcpT[pTlead][k][cent] = (TH1D *) list1->FindObject(
                            Form("hMCpT_pTl%i_R0%.0f_centbin%i", pTlead, R * 10, cent));
                    hMcpT[pTlead][k][cent]->Scale(1. / Nevents);

                    //hRcpT[pTlead][k][cent]=(TH1D*)list1->FindObject(Form("hfpT_pTl%i_R0%.0f_centbin%i", pTlead, R*10,cent));
                    //hRcpT[pTlead][k][cent]->Scale(1./Nevents);

                    //	hMcMatchedpT[pTlead][k][cent]=(TH1D*)list1->FindObject(Form("hMcmatchedpT_pTl%i_R0%.0f_cent%i", pTlead, R*10,cent));
                    hMcMatchedpT[pTlead][k][cent] = (TH1D *) list1->FindObject(
                            Form("hMCmatchedpT_pTl%i_R0%.0f_centbin%i", pTlead, R * 10, cent));
                    hMcMatchedpT[pTlead][k][cent]->Scale(1. / Nevents);
//				heffi[pTlead][k][cent]=(TH1D*)hMcMatchedpT[pTlead][k][cent]->Clone(Form("heffi_pTl%i_R0%.0f_cent%i", pTlead, R*10,cent));
//				heffi[pTlead][k][cent]->Divide(hMcpT[pTlead][k][cent]);		

                    hMcMatchedMCpT[pTlead][k][cent] = (TH1D *) list1->FindObject(
                            Form("hMCmatchedpT_MCpTl%i_R0%.0f_centbin%i", pTlead, R * 10, cent));
                    hMcMatchedMCpT[pTlead][k][cent]->Scale(1. / Nevents);

                    hRcMatchedpT[pTlead][k][cent] = (TH1D *) list1->FindObject(
                            Form("hRCmatchedpT_pTl%i_R0%.0f_centbin%i", pTlead, R * 10, cent));
                    hRcMatchedpT[pTlead][k][cent]->Scale(1. / Nevents);


                    c[pTlead][k][cent] = Form("hResponseMatrix_pTl%i_R0%.0f_centbin%i", pTlead, R * 10, cent);
                    //cf[j][k] = Form("hfjetpTlead_R0%.0f_centbin%i", R*10, centbin);
                    //} else {c[j][k] = Form("hjetpT_R0%.0f_centbin%i_corr", R*10, centbin);cf[j][k] = Form("hfjetpT_R0%.0f_centbin%i_corr", R*10, centbin);}


                    cout << "name: " << c[pTlead][k][cent] << endl;
                    hResponseM_tmp[pTlead][k][cent] = (TH2D *) list1->FindObject(c[pTlead][k][cent]);
                    if (hResponseM_tmp[pTlead][k][cent]->GetEntries() == 0) {
                        cout << "empty histogram: " << hResponseM_tmp[pTlead][k][cent]->GetName() << ", skipping! "
                             << endl;
                        continue;
                    }
                    hResponseM_tmp[pTlead][k][cent]->SetTitle("");
                    //hResponseM_tmp[pTlead][k][cent]->RebinY(10);
                    //hrescharged[j][k]->Rebin(5);
                    hResponseM_tmp[pTlead][k][cent]->Scale(1. / Nevents);


                    TH2D *rmatrix = new TH2D("hresponse", "hresponse", newbins, pTbinArray, newbins2, pTbinArray2);
                    TH2D *rmatrix_tmp = (TH2D *) hResponseM_tmp[pTlead][k][cent]->Clone(Form("hResponseMatrix_pTl%i_R0%.0f_%s", pTlead, R * 10, centrality[cent].Data()));
                    TH1D *hMCreco = new TH1D("hmcreco", "hmcreco", newbins, pTbinArray);
                    TH1D *hMCreco_tmp = (TH1D *) rmatrix_tmp->ProjectionX("hMCreco", 0, -1, "e"); //MC measured spectrum
                    TH1D *hMCtrue = new TH1D("hmctrue", "hmctrue", newbins2, pTbinArray2);
                    TH1D *hMCtrue_tmp = (TH1D *) rmatrix_tmp->ProjectionY("hMCtrue", 0, -1,"e"); //MC measured spectrum - matched

                    rmatrix->Sumw2();
                    hMCreco->Sumw2();
                    hMCtrue->Sumw2();
                    rmatrix = rebin_histogram2D(rmatrix_tmp, hMCreco, hMCtrue,Form("hResponseMatrix_pTl%i_R0%.0f_%s", pTlead, R * 10,centrality[cent].Data()));

//                delete rmatrix_tmp;
                    //   rmatrix->Write("hresponse");


                    rmatrix->SetXTitle("p_{T}^{det} [GeV/#it{c}]");
                    rmatrix->Draw("colz");
                    rmatrix->Write(Form("hResponseMatrix_pTl%i_R0%.0f_%s", pTlead, R * 10, centrality[cent].Data()));

                    //normalize each column in pTtrue to 1
                    hResponseM[pTlead][k][cent] = (TH2D *) hResponseM_tmp[pTlead][k][cent]->Clone(
                            Form("hResponseMatrix_pTl%i_R0%.0f_%s", pTlead, R * 10, centrality[cent].Data()));

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

                    if (!(!strcmp(centrality[cent].Data(), "peripheral") || !strcmp(centrality[cent].Data(), "central")))continue;

               //     if (hResponseM_tmp[pTlead][k][cent-1]->GetEntries()!=0) {
                    //    cout << "adding "<< hResponseM[pTlead][k][cent]->GetName() << " and " << hResponseM[pTlead][k][cent-1]->GetName()<<endl;
                        //adding
                    //    hResponseM[pTlead][k][cent]->Add(hResponseM[pTlead][k][cent-1]);
                        //hResponseM[pTlead][k][cent]->Write();

                        //heffi[pTlead][k][cent]->Add(heffi[pTlead][k][cent-1]);
                        //heffi[pTlead][k][cent]->Write();

                    //if (centbin == 7|| centbin == 8) hrescharged[j][k]->Scale(2);
                    //hResponseM[pTlead][k][cent]->Write();

                    //heffi[pTlead][k][cent]->Write();

                //    hMcpT[pTlead][k][cent]->Add(hMcpT[pTlead][k][cent-1]);
                    hMcpT[pTlead][k][cent]->Write(Form("hMCpT_pTl%i_R0%.0f_%s", pTlead, R * 10, centrality[cent].Data()));
                //    hMcMatchedpT[pTlead][k][cent]->Add(hMcMatchedpT[pTlead][k][cent-1]);
                    hMcMatchedpT[pTlead][k][cent]->Write(Form("hMcMatchedpT_pTl%i_R0%.0f_%s", pTlead, R * 10, centrality[cent].Data()));
                //    hMcMatchedMCpT[pTlead][k][cent]->Add(hMcMatchedMCpT[pTlead][k][cent-1]);
                    hMcMatchedMCpT[pTlead][k][cent]->Write(Form("hMcMatchedMCpT_pTl%i_R0%.0f_%s", pTlead, R * 10, centrality[cent].Data()));
                //    hRcMatchedpT[pTlead][k][cent]->Add(hRcMatchedpT[pTlead][k][cent-1]);
                    hRcMatchedpT[pTlead][k][cent]->Write(Form("hRcMatchedpT_pTl%i_R0%.0f_%s", pTlead, R * 10, centrality[cent].Data()));
            //    }
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
