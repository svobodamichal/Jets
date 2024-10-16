void plot(double R=0.3, int pTlead=5, TString evo="omicron", TString ext="pdf", int nfiles=1)
{
    int rebin=1;
    double TAA=1;
	TString pythia_file=Form("pythia_BGdete_%i.root",nfiles);
	TString dpT_file=Form("response_matrix_BG_sp_R%.1lf_pTlead%i.root",R,pTlead);
	TFile* fin_pyt=new TFile(pythia_file.Data(), "OPEN");
	TFile* fin_dpT=new TFile(dpT_file.Data(), "OPEN");

    TString suffix=""; //MB pythia spectrum
    if(nfiles==1) suffix="_sum_fudge"; //pythia spectrum merged from pThardbins
    
	//load pythia histograms with various detector effect
	TH1D* hpythia0=(TH1D*) fin_pyt->Get(Form("hchjet_pT_R0%.0lf_pTl0%s",R*10,suffix.Data()));
    TH1D* hpythia=(TH1D*) fin_pyt->Get(Form("hchjet_pT_R0%.0lf_pTl%i%s",R*10,pTlead,suffix.Data()));
	TH1D* hpythia_effi=(TH1D*) fin_pyt->Get(Form("hchjet_effi_pT_R0%.0lf_pTl%i%s",R*10,pTlead,suffix.Data()));
	TH1D* hpythia_effi_pTsmear=(TH1D*) fin_pyt->Get(Form("hchjet_effi_pTsmear_pT_R0%.0lf_pTl%i%s",R*10,pTlead,suffix.Data()));
	//TH1D* hpythia_effi_pTsmear_dpT=(TH1D*) fin_pyt->Get(Form("hchjet_effi_pTsmear_dpT_pT_R0%.0lf_pTl%i",R*10,pTlead); //this histogram is the same as hpythia_effi_pTsmear, dpT-smearing will be applied later

    //rebin and divide by the number of merged files so the histograms are normalized per event
    hpythia0->Rebin(rebin);
    hpythia0->Scale(TAA/(rebin*nfiles));
    hpythia->Rebin(rebin);
    hpythia->Scale(TAA/(rebin*nfiles));
    hpythia_effi->Rebin(rebin);
    hpythia_effi->Scale(TAA/(rebin*nfiles));
    hpythia_effi_pTsmear->Rebin(rebin);
    hpythia_effi_pTsmear->Scale(TAA/(rebin*nfiles));
    
	//load response matrix
	TH2D* hresponse=(TH2D*) fin_dpT->Get("hResponse_1E9");
	

	//smear pythia distribution with delta-pT
	TH1D* hpythia_effi_pTsmear_dpT=(TH1D*) smear_with_dpT(hpythia_effi_pTsmear,hresponse,"hpythia_effi_pTsmear_dpT");

   	//-----------------------------------------------
	//Draw
   	//-----------------------------------------------

	TCanvas* c1=new TCanvas("c1","",10,10,1200,900);
	c1->cd();
    c1->SetLogy();
    
    float xmin=-5;
    float xmax=30;
    float ymin=1E-9;
    float ymax=3E-4;
    
	TH1D* hframe=new TH1D("hframe","",40,xmin-5,xmax+5);
	hframe->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
	hframe->GetYaxis()->SetTitle("1/N_{events} 1/2#pi d^{2}N/dp_{T}^{ch}d#eta (GeV/#it{c})^{-1}");
    hframe->GetYaxis()->SetRangeUser(ymin,ymax);
	hframe->Draw();

    int lwidth=2;
    hpythia0->SetLineColor(kOrange+3);
    hpythia0->SetLineWidth(lwidth);
	//hpythia0->Draw("histo l same");
    
    hpythia->GetXaxis()->SetRangeUser(pTlead,xmax);
    hpythia->SetLineWidth(lwidth);
    hpythia->SetLineStyle(7);
    hpythia->SetLineColor(kGreen+3);
    /*hpythia->SetMarkerColor(kGreen+3);
    hpythia->SetMarkerSize(2);
    hpythia->SetMarkerStyle(21);*/
	hpythia->Draw("histo c same");
    
    hpythia_effi_pTsmear->GetXaxis()->SetRangeUser(pTlead,xmax);
	hpythia_effi_pTsmear->SetLineColor(kBlue);
    hpythia_effi_pTsmear->SetLineWidth(lwidth);
    hpythia_effi_pTsmear->SetLineStyle(1);
    /*hpythia_effi->SetMarkerColor(kBlue);
    hpythia_effi->SetLineColor(kBlue);
    hpythia_effi->SetMarkerSize(2);
    hpythia_effi->SetMarkerStyle(34);*/
	hpythia_effi_pTsmear->Draw("histo c same");
    
    /*
    hpythia_effi_pTsmear->GetXaxis()->SetRangeUser(pTlead,xmax);
	hpythia_effi_pTsmear->SetLineColor(kWhite);
    hpythia_effi_pTsmear->SetLineWidth(0);
    hpythia_effi_pTsmear->SetLineStyle(3);
    hpythia_effi_pTsmear->SetMarkerColor(kMagenta+2);
    hpythia_effi_pTsmear->SetMarkerSize(2);
    hpythia_effi_pTsmear->SetMarkerStyle(29);
	hpythia_effi_pTsmear->Draw("p same");*/
    
    hpythia_effi_pTsmear_dpT->GetXaxis()->SetRangeUser(xmin,xmax);
	hpythia_effi_pTsmear_dpT->SetLineColor(kRed);
    hpythia_effi_pTsmear_dpT->SetLineWidth(lwidth);
    hpythia_effi_pTsmear_dpT->SetLineStyle(2);
	hpythia_effi_pTsmear_dpT->Draw("histo c same");

    //legends
    TLegend *leg_info = new TLegend(0.2, 0.65, 0.45, 0.900);
	leg_info->SetTextSize(0.04);
	leg_info->SetFillStyle(0);
	leg_info->SetBorderSize(0);
	leg_info->SetMargin(0.05);
    leg_info->AddEntry("", "PYTHIA 6.4.28", "");
	leg_info->AddEntry("", "charged jets", "");
    leg_info->AddEntry("", Form("p_{T, lead}^{min} =  %i GeV/c", pTlead), "");
    leg_info->AddEntry("", Form("R = %.1lf", R), "");
	leg_info->DrawClone("same");
    
    TLegend *leg_data = new TLegend(0.6, 0.65, 0.97, 0.900);
	leg_data->SetTextSize(0.04);
	leg_data->SetFillStyle(0);
	leg_data->SetBorderSize(0);
	leg_data->SetMargin(0.2);
    //leg_data->AddEntry(hpythia0," p_{T}^{lead} > 0 GeV/c","l");
    leg_data->AddEntry(hpythia,"particle level","l");
    //leg_data->AddEntry(hpythia_effi," #otimes instr.","l");
    leg_data->AddEntry(hpythia_effi_pTsmear," #otimes instr. eff.","l");
    leg_data->AddEntry(hpythia_effi_pTsmear_dpT," #otimes instr. eff. + #delta_{p_{T}}","l");
    leg_data->DrawClone("same");
    
    c1->SaveAs(Form("../../../plotting_out/obr/%s/effi/PYTHIA_EffSize_R0%.0lf_pTlead%i_central.%s",evo.Data(),R*10,pTlead,ext.Data()));
}

//function for smearing pythia distribution with delta-pT
TH1D* smear_with_dpT(TH1D* hvector, TH2D* hmatrix, TString hsmear_name="hpythia_effi_pTsmear_dpT")
{

	TH1D* hsmeared=(TH1D*) hmatrix->ProjectionX(hsmear_name, 1,1); //create empty histogram with the same binning as the x axis of the response matrix
	hsmeared->Reset("MICE");
	
	for(int bin=1;bin<=hvector->GetNbinsX();bin++)
	{
		double pTtrue=hvector->GetBinCenter(bin);
		double val_true=hvector->GetBinContent(bin);
		if(!val_true>0) continue;
		int bin_true=hmatrix->GetYaxis()->FindBin(pTtrue);
		
		//for each pTtrue find the corresponding delta-pT histogram and use it for smearing
		TH1D* hdeltapT=hmatrix->ProjectionX(Form("hdpT_%i",bin_true),bin_true,bin_true);
        //cout<<"hdeltapT integral:"<<hdeltapT->Integral()<<" with width:"<<hdeltapT->Integral("width")<<endl;
		for (int bn=1;bn<=hsmeared->GetNbinsX();bn++)
		{
			double old_value=hsmeared->GetBinContent(bn);
			double new_value=old_value+val_true*hdeltapT->GetBinContent(bn);
			hsmeared->SetBinContent(bn,new_value);
			/*
            //calculate (approximate) errors
			double new_err=TMath::Sqrt(new_value*nevents)/nevents; //the hvector is normalized per event, we have to multiply it by the number of events to get counts
			hsmeared->SetBinError(bn,new_err);*/
		}	
	}

	//cout<<"hvector integral:"<<hvector->Integral()<<" with width:"<<hvector->Integral("width")<<endl;
	//cout<<"hsmeared integral:"<<hsmeared->Integral()<<" with width:"<<hsmeared->Integral("width")<<endl;
    //hsmeared->Scale(hvector->Integral("width")/hsmeared->Integral()); //the smeared histogram should have the same integral as the original one
	return hsmeared;
}

	
