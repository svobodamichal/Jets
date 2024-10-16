//******************************************************
//track qa plots
//******************************************************
void qa_tracks(TString name1="dca3",TString ext="pdf")
{

	TString inPath="../plotting_out/root/qa";
	TString outPath="../plotting_out/obr/intersteps/qa";
	TFile* f1=new TFile(Form("%s/qa_PvsG_AuAu_%s.root",inPath.Data(),name1.Data()));
	
	TH2D* h2npoints=(TH2D*) f1->Get("htrNpoints");
	TH1D* h1npoints=(TH1D*)h2npoints->ProjectionX("hNpoints",1,h2npoints->GetYaxis()->GetNbins());
	h1npoints->SetTitle("");
	h1npoints->GetXaxis()->SetTitle("N_{fit}");
	h1npoints->GetYaxis()->SetTitle("arb. units");
	h1npoints->GetXaxis()->SetRangeUser(15,50);
	h1npoints->GetXaxis()->SetTitleSize(0.06);
	h1npoints->GetYaxis()->SetTitleSize(0.055);
	//h1npoints->Scale(1.0/hdca->Integral("width"));

	TH1D* hdca=(TH1D*) f1->Get("hdca_tr");
	hdca->SetTitle("");
	hdca->GetXaxis()->SetTitle("DCA (cm)");
	hdca->GetYaxis()->SetTitle("arb. units");
	hdca->GetXaxis()->SetRangeUser(15,50);
	hdca->GetXaxis()->SetTitleSize(0.06);
	hdca->GetYaxis()->SetTitleSize(0.055);
	hdca->GetXaxis()->SetRangeUser(0,3);
	//hdca->Scale(1.0/hdca->Integral("width"));

	TH2D* hetaphi=(TH2D*) f1->Get("heta_phi_tr");
	hetaphi->SetTitle("");
	hetaphi->GetYaxis()->SetTitle("#phi (rad)");
	hetaphi->GetXaxis()->SetTitleSize(0.06);
	hetaphi->GetYaxis()->SetTitleSize(0.06);

	TCanvas *c1=new TCanvas("c1","c1",10,10,800,600);
	h1npoints->Draw();
	c1->SaveAs(Form("%s/npoints.%s",outPath.Data(),ext.Data()));

	TCanvas *c2=new TCanvas("c2","c2",10,10,800,600);
	c2->cd();
	c2->SetLogy();
	hdca->Draw();
	c2->SaveAs(Form("%s/dca.%s",outPath.Data(),ext.Data()));

	TCanvas *c3=new TCanvas("c3","c3",10,10,800,600);
	c3->cd();
	hetaphi->Draw("COLZ");
	c3->SaveAs(Form("%s/eta_phi.%s",outPath.Data(),ext.Data()));


}


//******************************************************
//event qa plots
//******************************************************
void qa_event(TString cclass="cent")
{
	TString path="../plotting_out/root/qa";
	TFile* f1=new TFile(Form("%s/qa_%s_%s.root",path.Data(),cclass.Data(),"nocuts"),"OPEN");
	TFile* f2=new TFile(Form("%s/qa_%s_%s.root",path.Data(),cclass.Data(),"cuts"),"OPEN");
	TTree* t1=(TTree*) f1->Get("eventTree");
	TTree* t2=(TTree*) f2->Get("eventTree");
	
	TH2D* hdayrefmult2D_1=(TH2D*) f1->Get("hday_refmult");
	TH2D* hdayrefmult2D_2=(TH2D*) f2->Get("hday_refmult");

	//histogram binning
	int runbins=50000;
	int runmin=125000.5;
	int runmax=175000.5;

	//histogram axis ranges
	float rmymin=(cclass="cent") ? 400 : 10; 
	float rmymax=(cclass="cent") ? 500 : 50; 
   float zymin=-15;
	float zymax=15;

	//mean refmult vs runid
	TH2D* hrefmultVSrunid1=new TH2D("hrefmultVSrunid1","refmult vs runid", 150,350,650,runbins,runmin,runmax);
	TH1D* hmeanrefmult1=new TH1D("hmeanrefmult1","mean refmult vs runid",runbins,runmin,runmax);
	TH2D* hrefmultVSrunid2=new TH2D("hrefmultVSrunid2","refmult vs runid", 150,350,650,runbins,runmin,runmax);
	TH1D* hmeanrefmult2=new TH1D("hmeanrefmult2","mean refmult vs runid",runbins,runmin,runmax);

	//mean z vs runid
	TH2D* hzVSrunid1=new TH2D("hzVSrunid1","z vs runid", 60,-30,30,runbins,runmin,runmax);
	TH1D* hmeanz1=new TH1D("hmeanz1","mean z vs runid",runbins,runmin,runmax);
	TH2D* hzVSrunid2=new TH2D("hzVSrunid2","z vs runid", 60,-30,30,runbins,runmin,runmax);
	TH1D* hmeanz2=new TH1D("hmeanz2","mean z vs runid", runbins,runmin,runmax);

	//mean refmult vs day
	TH1D* hmeanrefmultday1=new TH1D("hmeanrefmultday1","mean refmult vs day", 250,0,250);
	TH1D* hmeanrefmultday2=new TH1D("hmeanrefmultday2","mean refmult vs day", 250,0,250);
	
	TH2D* hmeanrefmult_day=new TH2D("hmeanrefmultday2D","mean refmult vs day", 100,400,500,250,0,250);
	
	int nentries1=t1->GetEntries();
	int nentries2=t2->GetEntries();
	//nentries1=1000000;
	//nentries2=1000000;
	int runid1,runid2;
	float refmult1,z1,zvpd1,refmult2,z2,zvpd2;

	t1->SetBranchAddress("refmult",&refmult1);
	t1->SetBranchAddress("runid",&runid1);
	t1->SetBranchAddress("zVPD",&zvpd1);
	t1->SetBranchAddress("zvertex",&z1);
	
	t2->SetBranchAddress("refmult",&refmult2);
	t2->SetBranchAddress("runid",&runid2);
	t2->SetBranchAddress("zVPD",&zvpd2);
	t2->SetBranchAddress("zvertex",&z2);
	
	cout<<"reading trees"<<endl;
	for(int i=0;i<nentries1;i++)
	{
		t1->GetEntry(i);
		//cout<<refmult<<endl;
		int run1=(runid1%1000000);
		hrefmultVSrunid1->Fill(refmult1,run1);
		hzVSrunid1->Fill(z1,run1);
		
	}
	
	for(int i=0;i<nentries2;i++)
	{
		t2->GetEntry(i);
		//cout<<refmult<<endl;
		int run2=(runid2%1000000);
		if(runid2==12154043)continue;
		//if(runid2==12156030)continue;
		hrefmultVSrunid2->Fill(refmult2,run2);
		hzVSrunid2->Fill(z2,run2);
		
	}
	
	cout<<"calculating mean values"<<endl;
	for(int i=1; i<=hrefmultVSrunid1->GetYaxis()->GetNbins(); i++)
	{
		TH1D* hslicerm1=hrefmultVSrunid1->ProjectionX(Form("hslicerm1_%i",i),i,i);
		double mrefmult=hslicerm1->GetMean();
		hmeanrefmult1->SetBinContent(i,mrefmult);
		
		int day=hmeanrefmult1->GetBinCenter(i)/1000;
		hmeanrefmult_day->Fill(mrefmult,day);

		
		TH1D* hslicerm2=hrefmultVSrunid2->ProjectionX(Form("hslicerm2_%i",i),i,i);
		mrefmult=hslicerm2->GetMean();
		hmeanrefmult2->SetBinContent(i,mrefmult);
		delete hslicerm1;
		delete hslicerm2;

		TH1D* hslicez1=hzVSrunid1->ProjectionX(Form("hslicez1_%i",i),i,i);
		double mz=hslicez1->GetMean();
		if(mz!=0)hmeanz1->SetBinContent(i,mz);
		else hmeanz1->SetBinContent(i,-9999);

		TH1D* hslicez2=hzVSrunid2->ProjectionX(Form("hslicez2_%i",i),i,i);
		mz=hslicez2->GetMean();
		if(mz!=0)hmeanz2->SetBinContent(i,mz);
		else hmeanz2->SetBinContent(i,-9999);
		if(mz>5)cout<<"high <z> in run #"<<hmeanz2->GetBinCenter(i)<<endl; 
		if(mz<-5)cout<<"low <z> in run #"<<hmeanz2->GetBinCenter(i)<<endl; 
		
		delete hslicez1;
		delete hslicez2;

	}

	int nbins=hdayrefmult2D_1->GetXaxis()->GetNbins();
	for(int i=1; i<=nbins; i++)
	{
		TH1D* hslice_RefMult=hdayrefmult2D_1->ProjectionY(Form("hsliceRM_%i",i),i,i);
		TH1D* hslice_meanRefMult=hmeanrefmult_day->ProjectionX(Form("hsliceMRM_%i",i),i,i);

		double mrefmult=hslice_RefMult->GetMean();
		double srefmult=hslice_meanRefMult->GetStdDev();
		hmeanrefmultday1->SetBinContent(i,mrefmult);
		hmeanrefmultday1->SetBinError(i,srefmult);
		
		TH1D* hslice2=hdayrefmult2D_2->ProjectionY(Form("hslice2_%i",i),i,i);
		mrefmult=hslice2->GetMean();
		hmeanrefmultday2->SetBinContent(i,mrefmult);
	}

	nbins=hmeanrefmult_day->GetYaxis()->GetNbins();
	TH1D* hslice_all=hmeanrefmult_day->ProjectionX("hslice_all",1,nbins);
	double refmult_mean=hslice_all->GetMean();
	double refmult_sigma=hslice_all->GetStdDev();
	
	
	
	//Draw histograms
	TCanvas* c2=new TCanvas("c2","refmult",10,10,1200,800);
	c2->cd();
	hmeanrefmultday1->GetYaxis()->SetTitle("<refMult>");
	hmeanrefmultday1->GetXaxis()->SetTitle("day");
	hmeanrefmultday1->Draw("");
	hmeanrefmultday2->SetLineColor(kRed);
	hmeanrefmultday2->Draw("histo same");
	
	  TLine *mean = new TLine(hmeanrefmultday1->GetXaxis()->GetBinCenter(1),refmult_mean, hmeanrefmultday1->GetXaxis()->GetBinCenter(hmeanrefmultday1->GetNbinsX()), refmult_mean);
    mean->SetLineWidth(2);
    mean->SetLineStyle(1);
    mean->SetLineColor(kGray+2);
    mean->DrawClone("same");
	 
	TLine *oneSigma1 = new TLine(hmeanrefmultday1->GetXaxis()->GetBinCenter(1),refmult_mean-refmult_sigma, hmeanrefmultday1->GetXaxis()->GetBinCenter(hmeanrefmultday1->GetNbinsX()), refmult_mean-refmult_sigma);
    oneSigma1->SetLineWidth(1);
    oneSigma1->SetLineStyle(2);
    oneSigma1->SetLineColor(kGray+1);
    oneSigma1->DrawClone("same");


	TCanvas* c1=new TCanvas("c1","refmult",10,10,1200,800);
	c1->cd();
	hmeanrefmult1->SetTitle("");
	hmeanrefmult1->GetYaxis()->SetTitle("<refMult>");
	hmeanrefmult1->GetXaxis()->SetTitle("run index");
	hmeanrefmult1->SetAxisRange(rmymin, rmymax,"Y");
	hmeanrefmult1->SetMarkerColor(kBlue);
	hmeanrefmult1->SetMarkerStyle(7);
	hmeanrefmult1->Draw("P");
	hmeanrefmult2->SetLineColor(kRed);
	hmeanrefmult2->SetMarkerColor(kRed);
	hmeanrefmult2->SetMarkerStyle(7);
	//hmeanrefmult2->SetMarkerSize(1);
	hmeanrefmult2->Draw("same P");

	//fake histograms for Legend
	TH1D* hfake1=new TH1D("hfake1","",1,0,1);
	TH1D* hfake2=new TH1D("hfake2","",1,0,1);
	hfake1->SetMarkerColor(kBlue);
	hfake1->SetMarkerStyle(21);
	hfake2->SetMarkerColor(kRed);
	hfake2->SetMarkerStyle(21);

	//Legend
	TLegend *legspectra = new TLegend(0.6, 0.75, 0.90, 0.90);
	legspectra->SetTextSize(0.045);
	legspectra->SetFillStyle(0);
	legspectra->SetBorderSize(0);
	legspectra->AddEntry(hfake1, "all events", "p");
	legspectra->AddEntry(hfake2, "accepted events", "p");
	legspectra->DrawClone("same");

	TCanvas* c3=new TCanvas("c3","z",10,10,1200,800);
	c3->cd();
	hmeanz1->SetTitle("");
	hmeanz1->GetYaxis()->SetTitle("<z> (cm)");
	hmeanz1->GetXaxis()->SetTitle("run index");
	hmeanz1->SetAxisRange(zymin, zymax,"Y");
	hmeanz1->SetMarkerColor(kBlue);
	hmeanz1->SetMarkerStyle(7);
	hmeanz1->Draw("P");
	hmeanz2->SetLineColor(kRed);
	hmeanz2->SetMarkerColor(kRed);
	hmeanz2->SetMarkerStyle(7);
	hmeanz2->Draw("histo same P");

	TLegend *legspectra = new TLegend(0.6, 0.75, 0.90, 0.90);
	legspectra->SetTextSize(0.045);
	legspectra->SetFillStyle(0);
	legspectra->SetBorderSize(0);
	legspectra->AddEntry(hfake1, "all events", "p");
	legspectra->AddEntry(hfake2, "accepted events", "p");
	legspectra->DrawClone("same");

	TString ext[]={"pdf","gif"};
	for(int i=0;i<2;i++)
	{
		c1->SaveAs(Form("../plotting_out/obr/intersteps/event_qa/meanrefmult_%s.%s",cclass.Data(),ext[i].Data()));
		c3->SaveAs(Form("../plotting_out/obr/intersteps/event_qa/meanz_%s.%s",cclass.Data(),ext[i].Data()));
	}


}
