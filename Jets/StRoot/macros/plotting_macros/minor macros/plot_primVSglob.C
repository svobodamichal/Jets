void plot_primVSglob_track(TString fname="main")
{
	double new_bins[16]={0,1,2,3,4,5,6,7,8,9,11,13,15,20,25,30};
	TFile* f1=new TFile(Form("../plotting_out/root/qa/qa_PvsG_AuAu_%s.root",fname.Data()),"OPEN");
	TH1D* hprim=f1->Get("hpT_tr");
	TH1D* hglob=f1->Get("hpT_tr_glob");
	TH1I* hevents=f1->Get("hevents");
	hprim->Sumw2();
	hglob->Sumw2();

	TH2D* hPG=(TH2*)f1->Get("hpT_prim_glob");
	
	TH1D* hp=hprim->Rebin(15,"hp",new_bins);
	TH1D* hg=hglob->Rebin(15,"hg",new_bins);
	
	TCanvas* c1=new TCanvas("c1","c1",10,10,800,1200);
	c1->Divide(1,2);
	c1->cd(1);
	gPad->SetLogy();
	hprim->SetTitle("primary vs global tracks in Run11");
	hprim->Draw("");
	hglob->SetLineColor(kRed);
	hglob->Draw("same");
	
	TLegend *legspectra = new TLegend(0.65, 0.75, 0.9, 0.90);
	legspectra->SetTextSize(0.03);
	legspectra->SetFillStyle(0);
	legspectra->SetBorderSize(0);
	legspectra->AddEntry(hprim,"primary tracks","l");
	legspectra->AddEntry(hglob,"global tracks","l");
	legspectra->DrawClone("same");
	
	c1->cd(2);
	hp->Divide(hg);
	hp->SetTitle("");
	hp->GetYaxis()->SetRangeUser(0,1.5);
	hp->GetYaxis()->SetTitle("primary/global");
	hp->GetXaxis()->SetTitle("p_{T}^{track}");
	hp->Draw("e");

	TLatex *latex = new TLatex();
   latex->SetNDC();
   latex->SetTextSize(0.035);
	latex->DrawLatex(0.4, 0.7,Form("%.1lfM 0-10%% central events",hevents->GetEntries()/1E6));

	c1->SaveAs(Form("../plotting_out/obr/prim_vs_glob/track_ratio_%s.gif",fname.Data()));
	
	TCanvas* c2=new TCanvas("c2","c2",10,10,800,600);
	c2->cd();
	c2->SetLogz();
	hPG->Draw("COLZ");
	c2->SaveAs(Form("../plotting_out/obr/prim_vs_glob/tracks_PvsG_%s.gif",fname.Data()));
	

}
void plot_primVSglob_jet(float R=0.3,float pTlead=5.0)
{
	double new_bins[26]={-30,-20,-10,-5,-2,0,1,2,3,4,5,6,7,8,9,11,13,15,17,20,25,30,35,40,50,60};
	TFile* f1=new TFile("../plotting_out/root/MB/phi/histos_inclusivejet_normal.root","OPEN");
	TFile* f2=new TFile("../plotting_out/root/MB/gamma/histos_inclusivejet_normal.root","OPEN");

	TH1D* hprim=f1->Get(Form("hpT_pTl%.0lf_R0%.0lf",pTlead,R*10));
	TH1D* hglob=f2->Get(Form("hpT_pTl%.0lf_R0%.0lf",pTlead,R*10));
	TH1I* hevents1=f1->Get("hevents");
	TH1I* hevents2=f2->Get("hevents");
	int nevents1=hevents1->GetEntries();
	int nevents2=hevents2->GetEntries();
	hprim->Sumw2();
	hglob->Sumw2();

	TH1D* hp=hprim->Rebin(25,"hp",new_bins);
	TH1D* hg=hglob->Rebin(25,"hg",new_bins);

	hp->Scale(1.0/nevents1,"width");
	hg->Scale(1.0/nevents2,"width");
	
	TCanvas* c1=new TCanvas("c1","c1",10,10,800,600);
	c1->cd();
	c1->SetLogy();
	hp->SetTitle("primary vs global tracks in Run11");
	hp->SetLineColor(kBlue);
	hg->SetLineColor(kRed);
	hg->Draw("e");
	hp->Draw("e same");

	TLegend *legspectra = new TLegend(0.65, 0.75, 0.9, 0.90);
	legspectra->SetTextSize(0.03);
	legspectra->SetFillStyle(0);
	legspectra->SetBorderSize(0);
	legspectra->AddEntry(hp,"primary track jets","l");
	legspectra->AddEntry(hg,"global track jets","l");
	legspectra->DrawClone("same");
	
	TLatex *latex = new TLatex();
   latex->SetNDC();
   latex->SetTextSize(0.035);
	latex->DrawLatex(0.3, 0.30,Form("%.1lfM 0-10%% central events",hevents1->GetEntries()/1E6));
	latex->DrawLatex(0.3, 0.22,Form("R=%.1lf, p_{T}^{lead}>%.0lf GeV/c",R,pTlead));
	c1->SaveAs("../plotting_out/obr/prim_vs_glob/jet_spectra.gif");

		
	TCanvas* c2=new TCanvas("c2","c2",10,10,800,600);
	c2->cd();
	TH1D* hd=(TH1D*) hp->Clone("hratio");
	hd->Divide(hg);
	hd->SetTitle("primary vs global tracks in Run11");
	hd->GetYaxis()->SetRangeUser(0,1.5);
	hd->GetYaxis()->SetTitle("jets with primary tr./global tr.");
	hd->GetXaxis()->SetTitle("p_{T,raw}^{jet}");
	hd->Draw("e");

	latex->DrawLatex(0.4, 0.75,Form("%.1lfM 0-10%% central events",hevents1->GetEntries()/1E6));
	latex->DrawLatex(0.4, 0.68,Form("R=%.1lf, p_{T}^{lead}>%.0lf GeV/c",R,pTlead));

	c2->SaveAs("../plotting_out/obr/prim_vs_glob/jet_ratio.gif");

}
