#include "include/eff_functions.C"
#include "include/hadron_ratio.C"

Color_t color[]={kRed,kBlue,kMagenta,kGreen+1,kGray,kRed+1,kOrange,kCyan+1,kGray+2};
int marker[]={21,22,23,24,25,26,27,28,29};

//=================================
//momentum resolution
//=================================
void momentum_resolution(TString particle="Pi", short central=1, TString cutset="nfit15")
{
	//how many pTMC values will we use
	int show=7;
	if(particle=="P") show=4;
	else if(particle=="jet") show=9;
	const int nshow=show;
	double pTshow[]={1.6,2.0,2.5,4.0,5,6,7.5,10,14};
	double pTshow_glob[]={1.6,2.5,4.0,6,8,10,14}; //for global tracks
	double pTshowMin=1.9; //show only resolution distributions for pTMC above this value
	double fsigma[nshow];  //fitted values of sigma parameter
	double fsigerr[nshow];  //sigma parameter error
	TFile* fin=new TFile(Form("root/%s/MC_%s_cent%i.root",cutset.Data(),particle.Data(),central),"OPEN");
	TString fout=Form("fig/%s/MC_%s_cent%i",cutset.Data(),particle.Data(),central);
	TString fout2=Form("fig/%s/sigma_dependency_%s_cent%i",cutset.Data(),particle.Data(),central);
	TH2D* hpTrec_pTmc=(TH2D*) fin->Get("hpTRec_pTMc");
	TH1D* hres[nshow];
	TF1* fres[nshow];
	TString centrality[]={"peripheral (60-80%)","central (0-10%)",""};
	TString globprim="primary";
	TString system="Au+Au #sqrt{s_{NN}}=200 GeV";
	TString trackTypes="#pi^{+}, #pi^{-}";
	float xmax[]={10,10,10,10,15,20,25,25,25};
	float ymin=0.001;
	float ymax=3;
	if(particle=="K") trackTypes="K^{+}, K^{-}";
	else if(particle=="P") trackTypes="p^{+}, p^{-}";
	else if(particle=="jet") 
	{
		trackTypes="ch. hadrons";
		system="pp #sqrt{s}=200 GeV";
	}
	else if(particle=="Pi0") 
	{
		trackTypes="e^{+},e^{-}";
		system="pp #sqrt{s}=200 GeV";
		globprim="global";
		for (int i=0;i<nshow;i++)
		{
			pTshow[i]=pTshow_glob[i];
		}
		ymax=8;
	}
	
	for(int i=0; i<nshow; i++)
	{
		int bin=hpTrec_pTmc->GetYaxis()->FindBin(pTshow[i]);
		hres[i]=(TH1D*) hpTrec_pTmc->ProjectionX(Form("hslice_%i",i),bin,bin);
		hres[i]->SetLineColor(color[i]);
		hres[i]->SetMarkerStyle(marker[i]);
		hres[i]->SetMarkerColor(color[i]);
		hres[i]->Sumw2();
		hres[i]->Rebin(2);
		hres[i]->Scale(1.0/hres[i]->Integral("width"));
		//cout<<"i:"<<hres[i]->Integral()<<endl;
		//cout<<"iw:"<<hres[i]->Integral("width")<<endl;
		
		fres[i]=new TF1(Form("fres%i",i),"gaus(0)",0,25);
		double res_percent=0.005;
		if(particle=="Pi0")res_percent=0.01;
		double sigma=res_percent*pTshow[i]*pTshow[i]; //this is just initial value for the fit
		//double sigma=0.01*pTshow[i]*TMath::Sqrt(pTshow[i]);
		//cout<<"sigma:"<<sigma<<endl;
		fres[i]->SetParameters(1.0/(TMath::Sqrt(2*TMath::Pi())*sigma),pTshow[i],sigma);
		fres[i]->SetLineColor(color[i]);
		fres[i]->SetLineStyle(2);
		
		//cout<<fres[i]->Integral(0,15)<<endl;
	}
	
	TH1D* hframe=new TH1D("hframe","",1,0,xmax[show-1]);
	TCanvas* c1=new TCanvas("c1","resolution",10,10,800,600);
	c1->cd();
	c1->SetLogy();
	hframe->GetYaxis()->SetRangeUser(2E-3,1E2);
	hframe->GetXaxis()->SetTitle("p_{T}^{reco} (GeV/c)");
	hframe->GetYaxis()->SetTitle("arb. units");
	hframe->DrawCopy("");
	for(int i=0; i<nshow; i++)
	{
		if(pTshow[i]>pTshowMin) hres[i]->Draw("esame");
		hres[i]->Fit(fres[i],"0");
		cout<<"resolution:"<<fres[i]->GetParameter(2)/(pTshow[i]*pTshow[i])<<endl;
		if(pTshow[i]>pTshowMin) fres[i]->DrawCopy("same");
		fsigma[i]=fres[i]->GetParameter(2);
		fsigerr[i]=fres[i]->GetParError(2);
		cout<<"pT:"<<pTshow[i]<<"  sigma:"<<fsigma[i]<<" error: "<<fsigerr[i]<<endl;
	}
	
	fres[0]->SetLineColor(kBlack);
	
	TLegend *legspectra = new TLegend(0.65, 0.45, 0.9, 0.85);
	legspectra->SetTextSize(0.045);
	legspectra->SetFillStyle(0);
	legspectra->SetBorderSize(0);
	legspectra->SetHeader("p_{T}^{MC}:");
	for(int i=0; i<nshow; i++)
	{
		if(pTshow[i]>pTshowMin)legspectra->AddEntry(hres[i], Form("%.1lf GeV",pTshow[i]), "p");
	}
	legspectra->AddEntry(fres[0], "N(p_{T}^{MC},#sigma)", "l");
	//legspectra->AddEntry("", Form("#sigma=%.3lf*p_{T}^{2}",res_percent), "");
	legspectra->DrawClone("same");
  
	TLegend *leginfo = new TLegend(0.15, 0.70, 0.45, 0.90);
	leginfo->SetTextSize(0.045);
	leginfo->SetFillStyle(0);
	leginfo->SetBorderSize(0);
	leginfo->AddEntry("",Form("%s tracks (%s)",globprim.Data(),trackTypes.Data()),"");
	leginfo->AddEntry("",system,"");
	if(central<2)leginfo->AddEntry("",centrality[central],"");
	leginfo->DrawClone("same");
	  
	
	c1->SaveAs(Form("%s.gif",fout.Data()));
	c1->SaveAs(Form("%s.pdf",fout.Data()));
	
	TGraphErrors* gr = new TGraphErrors(nshow,pTshow,fsigma,0, fsigerr);
	TF1* fitpol=new TF1("fitpol","[0]+[1]*x+[2]*x*x",0,15);
	fitpol->SetParameters(0,0.005,0.005);
	fitpol->SetParLimits(2,0.001,0.1);
	//fitpol->SetParLimits(1,0,0.1);
	if(particle=="Pi0") 
	{
		fitpol->FixParameter(0,0);
		fitpol->FixParameter(1,0);
		fitpol->FixParameter(2,0.012);
	}
	TF1* fjet=new TF1("fitpol","[0]+[1]*x+[2]*x*x",0,15);
	fjet->SetParameters(-0.026,0.020,0.003);
	fjet->SetLineStyle(2);
	fjet->SetLineColor(kBlue);
	
	TF1* fjet2=new TF1("fitpol","[0]+[1]*x+[2]*x*x",0,15);
	fjet2->SetParameters(0,0,0.003);
	fjet2->SetLineStyle(2);
	fjet2->SetLineColor(kGray+1);

	TH1D* hframe2=new TH1D("hframe","",1,0,20);
	TCanvas* c2=new TCanvas("c2","resolution pT-dependency",10,10,800,600);
	c2->cd();
	c2->SetLogy();
	hframe2->GetYaxis()->SetRangeUser(ymin,ymax);
	hframe2->GetXaxis()->SetTitle("p_{T}^{MC} (GeV/c)");
	hframe2->GetYaxis()->SetTitle("#sigma (GeV/c)");
	hframe2->DrawCopy("");
	
	gr->Draw("same *P");
	gr->Fit(fitpol,"0");
	
	if(particle!="Pi0" && particle!="jet")
	{
		fjet->Draw("same");
		fjet2->DrawCopy("same");
	}
	if(particle=="jet") fitpol->SetLineColor(kBlue);
	fitpol->DrawCopy("same");

	leginfo->DrawClone("same");
	
	double afit=fitpol->GetParameter(0);
	double bfit=fitpol->GetParameter(1);
	double cfit=fitpol->GetParameter(2);
	
	double aerr=100*TMath::Abs(fitpol->GetParError(0)/afit);
	double berr=100*TMath::Abs(fitpol->GetParError(1)/bfit);
	double cerr=100*TMath::Abs(fitpol->GetParError(2)/cfit);
	
	TLatex *latex = new TLatex();
	latex->SetNDC();
	latex->SetTextSize(0.04);
	if(particle!="Pi0")
	{
		latex->DrawLatex(0.4, 0.32,Form("#sigma=%.3lf + %.3lf p_{T} + %.3lf p_{T}^{2}",afit,bfit,cfit));
		latex->DrawLatex(0.4, 0.25,Form("rel. error: %.0lf %, %.0lf %, %.0lf %",aerr,berr,cerr));
	}
	else
		latex->DrawLatex(0.4, 0.35,Form("#sigma = %.3lf p_{T}^{2}",cfit));
	
	TLegend *legspectra2 = new TLegend(0.6, 0.45, 0.85, 0.70);
	legspectra2->SetTextSize(0.045);
	legspectra2->SetFillStyle(0);
	legspectra2->SetBorderSize(0);
	if(particle!="Pi0" && particle!="jet")
	{
		legspectra2->AddEntry(fjet, "ch. hadrons", "l");
		legspectra2->AddEntry("", "(p+p run12)", "");
		legspectra2->AddEntry(fjet2, "#sigma=0.003*p_{T}^{2}", "l");
	}
	//legspectra->AddEntry("", Form("#sigma=%.3lf*p_{T}^{2}",res_percent), "");
	legspectra2->DrawClone("same");
	
	c2->SaveAs(Form("%s.gif",fout2.Data()));
	c2->SaveAs(Form("%s.pdf",fout2.Data()));	
}

//==============================
//tracking efficiency
//==============================
void efficiency(short central=1, TString cutset="nfit15",TString ratio="AuAu")
{
    bool showBothRatios=1; //show efficiencies for both AuAu and pp like hadron ratios
    bool showAlexOnly=1; //show only Alex's h+jet tracking efficiency
	const int npart=3;
	TString particle[npart]={"Pi","K","P"};
	TString centrality[]={"peripheral (60-80%)","central (0-10%)",""};
	TString system="Au+Au #sqrt{s_{NN}}=200 GeV";
	TString globprim="primary";
	TString fout2=Form("fig/%s/Efficiency_%s_cent%i",cutset.Data(),ratio.Data(),central);

	if(cutset=="cutset4" )globprim="global";
	
    //y-axis range
	 //float ymin=0.6;
    //float ymax=1.2; 
    float ymin=0.55;
    float ymax=0.9; 
    if(!central)
    {
        ymin=0.7;
        ymax=1.1; 
    }
	float xmax=8;
    
	//lines
	const float linewidth=2.0;
	Color_t colorList[]={kGray+1,kMagenta,kBlue,kGreen+2,kBlack};
	int lineStyle[]={1,3,2,4,5};
	
	//fills
	Color_t colorFill[]={kGray,kGray,kGray,kGray,kGray};
	int styleFill[]={3004,0,0,0,0};
	//markers
	const float marksize=1.0;
	int markerList[]={20,21,22,23,29,33,34};
	
	TString legDraw[]={"f","l","l","l","l"};
	TString legName[]={"total","#pi^{+}, #pi^{-}","K^{+}, K^{-}", "p^{+}, p^{-}"};

	TFile* fin[npart];
	TH1D* heff[npart];
	
	for(int pr=0; pr<npart; pr++)
	{
		fin[pr]=new TFile(Form("root/%s/MC_%s_cent%i.root",cutset.Data(),particle[pr].Data(),central),"OPEN");
		//efficiency
		TH1D* hrec=fin[pr]->Get("hpTMC_Rec"); //pTMC of matched tracks
		TH1D* hMC=fin[pr]->Get("hpTMC_MC"); //pTMC of generated tracks

		heff[pr]=(TH1D*) hrec->Clone(Form("heff_%i",pr));
		heff[pr]->Divide(hMC);
		heff[pr]->SetLineWidth(linewidth);
		heff[pr]->SetLineStyle(lineStyle[pr+1]);
		heff[pr]->SetMarkerStyle(marker[pr+1]);
		heff[pr]->SetLineColor(colorList[pr+1]);
		heff[pr]->SetMarkerColor(colorList[pr+1]);
		heff[pr]->SetFillColor(colorFill[pr+1]);
		heff[pr]->SetFillStyle(styleFill[pr+1]);
		delete hrec;
		delete hMC;
	}
	
	//calculate total hadron reconstruction efficiency by adding efficiencies for pi,K and p weighted by the corresponding hadron abundancies
	TH1D* heff_tot=heff[0]->Clone("heff_tot");
	heff_tot->Reset("MICE");
	for(int i=1;i<=heff_tot->GetNbinsX();i++)
	{
		double val=0;
		double pT=heff_tot->GetBinCenter(i);
		for(int pr=0; pr<npart; pr++)
		{
			int part=2*pr+1;
			double weight=hadron_ratio(part,pT,ratio)+hadron_ratio(part+1,pT,ratio); //particle+antiparticle
			if(pT>5 && pr==0){ weight=(hadron_ratio(1,pT,ratio)+hadron_ratio(2,pT,ratio)+hadron_ratio(5,pT,ratio)+hadron_ratio(6,pT,ratio));}
			//else if(pT>5 && pr=2) continue;
			val+=weight*heff[pr]->GetBinContent(i);
		}
		heff_tot->SetBinContent(i,val);
	}
	heff_tot->SetLineWidth(linewidth);
	heff_tot->SetLineStyle(lineStyle[0]);
	heff_tot->SetMarkerStyle(marker[0]);
	heff_tot->SetLineColor(colorList[0]);  
	heff_tot->SetMarkerColor(colorList[0]);
	heff_tot->SetFillColor(colorFill[0]);
	heff_tot->SetFillStyle(styleFill[0]);
	
	TH1D* hframe=new TH1D("hframe","",1,0,xmax);
	TCanvas* c1=new TCanvas("c1","efficency",10,10,800,600);
	c1->cd();
	hframe->GetYaxis()->SetRangeUser(ymin,ymax);
	hframe->DrawCopy("");
	for(int pr=0; pr<npart; pr++)
	{
		heff[pr]->Draw("same");
	}
	//f_Efficiency1->Draw("same");
	//f_Efficiency2->Draw("same");
	
	TCanvas* c2=new TCanvas("c2","total efficency",10,10,800,600);
	c2->cd();
	hframe->GetXaxis()->SetTitle("p_{T}^{MC} (GeV/c)");
	hframe->GetYaxis()->SetTitle("efficiency");
	hframe->DrawCopy("");
	if(!showAlexOnly)
    {
        heff_tot->Draw("histo same");
        for(int pr=0; pr<npart; pr++)
        {
            heff[pr]->Draw("histo same");
        }
    }
	TF1* feff_Alex=efficiency11(central,ratio);
	feff_Alex->SetLineColor(kMagenta);
	feff_Alex->Draw("same");
    
    TF1* feff_Alex_pp=efficiency11(central,"pp");
	feff_Alex_pp->SetLineColor(kRed);
    feff_Alex_pp->SetLineStyle(2);
	if(showBothRatios) feff_Alex_pp->Draw("same");
	
	TFile* feffif=new TFile(Form("include/eff_%s.root",ratio.Data()),"OPEN");
	TF1* feff_Stephen_low=feffif->Get("effhL");
	TF1* feff_Stephen_high=feffif->Get("effhH");
	feff_Stephen_low->SetLineColor(kBlue);
	feff_Stephen_high->SetLineColor(kBlue);
	feff_Stephen_low->SetRange(0,1.2);
	feff_Stephen_high->SetRange(1.2,10);
	//feff_Stephen_low->Draw("same");
	//feff_Stephen_high->Draw("same");
	
	float pTx=0.6;
    if(!central) pTx=0.55;
	TF1* f_Efficiency1 = new TF1("f_Efficiency1",Eff_track_rec_function,0,pTx,3);
	f_Efficiency1->SetParameters(7.45643e-01,1.43725e-01,2.02904e+00);
	f_Efficiency1->SetLineColor(kRed);
	if(!showAlexOnly)heff_tot->Fit(f_Efficiency1,"R");
	if(!showAlexOnly)f_Efficiency1->Draw("same");

	TF1* f_Efficiency2 = new TF1("f_Efficiency2",Eff_track_rec_function,pTx,8.0,3);
	f_Efficiency2->SetParameters(7.45643e-01,1.43725e-01,2.02904e+00);
	f_Efficiency2->SetLineColor(kRed);
	if(!showAlexOnly)heff_tot->Fit(f_Efficiency2,"R");
	if(!showAlexOnly)f_Efficiency2->Draw("same");
		
	TLegend *legspectra = new TLegend(0.20, 0.60, 0.4, 0.95);
	legspectra->SetTextSize(0.042);
	legspectra->SetFillStyle(0);
	legspectra->SetBorderSize(0);
	legspectra->SetHeader("");
    if(!showAlexOnly){
		legspectra->AddEntry(heff_tot, legName[0], legDraw[0]);
        for(int i=0; i<npart; i++)
        {
            legspectra->AddEntry(heff[i], legName[i+1], legDraw[i+1]);
        }
    }
	if(!showAlexOnly) legspectra->AddEntry(f_Efficiency1, "fit function", "l");
	if(!showAlexOnly) legspectra->AddEntry(feff_Alex, "h+jet eff. function", "l");
    else
    {   
        legspectra->AddEntry(feff_Alex, Form("%s-like ratio", ratio.Data()), "l");
        if(showBothRatios) legspectra->AddEntry(feff_Alex_pp, "pp-like ratio", "l");
    }
	//legspectra->AddEntry(feff_Stephen_high, "eff. function from Stephen", "l");
	legspectra->DrawClone("same");
  
	TLatex *latex = new TLatex();
	latex->SetNDC();
	latex->SetTextSize(0.04);
	TString ratio_desc=(ratio=="AuAu")? "p/K/#pi  AuAu-like" : "p/K/#pi pp-like";
	if(!showBothRatios) latex->DrawLatex(0.55,0.35,ratio_desc);
	
	TLegend *leginfo = new TLegend(0.40, 0.70, 0.90, 0.90);
	leginfo->SetTextSize(0.045);
	leginfo->SetFillStyle(0);
	leginfo->SetBorderSize(0);
	//leginfo->AddEntry("",Form("%s tracks (%s)",globprim.Data(),trackTypes.Data()),"");
	leginfo->AddEntry("",system,"");
	if(central<2)leginfo->AddEntry("",centrality[central],"");
	leginfo->DrawClone("same");

	c2->SaveAs(Form("%s.gif",fout2.Data()));
	c2->SaveAs(Form("%s.pdf",fout2.Data()));	
	
	cout<<"cent:"<<central<<" cutset:"<<cutset<<" ratio:"<<ratio.Data()<<endl;
	cout<<"par[6]={"<<f_Efficiency1->GetParameter(0)<<","<<f_Efficiency1->GetParameter(1)<<","<<f_Efficiency1->GetParameter(2)
	<<","<<f_Efficiency2->GetParameter(0)<<","<<f_Efficiency2->GetParameter(1)<<","<<f_Efficiency2->GetParameter(2)<<"};"<<endl;
	
}
