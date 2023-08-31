//============================================================================================
//============================================================================================

void pyDraw(int spec=0, float R=0.3, int charged=2, TString frag="g", TString ext="gif")
{
	TString chname[]={"full","charged",""}
	TString chname2=chname[charged];
	if(charged<2)chname2+="_";
	TString specshape[]={"flat","pT5"}
	TString finame=Form("root/pythia_%stest_%s_%.1lf_%s.root",chname2.Data(),specshape[spec].Data(),R,frag.Data());
	TFile *fin=new TFile(finame,"OPEN");
	
	Color_t colorList[]={kRed,kBlue,kMagenta};
	
	TH1D* hjet0=(TH1D*) fin->Get("hpTptlJet_pTl0");
	TH1D* hjet5=(TH1D*) fin->Get("hpTptlJet_pTl5");
	TH1D* h4v0=(TH1D*) fin->Get("hpTptl4v_pTl0");
	TH1D* h4v5=(TH1D*) fin->Get("hpTptl4v_pTl5");
	

	
	/*
	TH1D* hempty=(TH1D*)hjet0->Clone("hempty");
	hempty->Reset("MICE");
	*/
	
	TH2D* hpT4v_pTjet_pTl5=(TH2D*) fin->Get("hpT_pTdete_pTl5");
	int bin5r=hpT4v_pTjet_pTl5->GetXaxis()->FindBin(5.001);
	int bin7l=hpT4v_pTjet_pTl5->GetXaxis()->FindBin(5.999);
	int bin7r=hpT4v_pTjet_pTl5->GetXaxis()->FindBin(6.001);
   int lastbin=hpT4v_pTjet_pTl5->GetNbinsX();

	TString hp_name="hproject57";
	TH1D* hproject57=(TH1D*)hpT4v_pTjet_pTl5->ProjectionY(hp_name,bin5r,bin7l);
	hp_name="hproject7up";
	TH1D* hproject7up=(TH1D*)hpT4v_pTjet_pTl5->ProjectionY(hp_name,bin7r,lastbin);
	
	hproject57->SetLineColor(kGreen);
	hproject57->SetFillColor(kGreen);
	hproject7up->SetFillStyle(3001);
	hproject7up->SetFillColor(kGray+2);
	hproject7up->SetLineColor(kGray+2);

	
	hjet0->SetLineColor(kRed);
	hjet5->SetLineColor(kRed);
	h4v0->SetLineColor(kBlue);
	h4v5->SetLineColor(kBlue);
	
	hjet0->SetLineStyle(2);
	h4v0->SetLineStyle(2);
	hproject57->SetLineStyle(2);
	
	hproject7up->SetLineStyle(2);
	hjet5->SetLineStyle(0);
	h4v5->SetLineStyle(0);
	
	int can_x=600;
	int can_y=600;
	
	float spectraXmin=4;
	float spectraXmax=10;
	float spectraYmin=10;
	float spectraYmax=hjet5->GetBinContent(hjet5->FindBin(5.001))*10;
	
  TCanvas *spectra = new TCanvas("cspectra","cspectra",10,10,2*can_x,can_y);
  cspectra->Divide(2,1);
  cspectra->cd(1);
  //cspectra->SetGrid();
  gPad->SetLogy();
  
  hjet0->SetTitle("generated vs reconstructed");
  hjet0->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hjet0->DrawCopy();
    
  hjet5->DrawCopy("same");
  h4v0->DrawCopy("same");
  h4v5->DrawCopy("same");
  
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.03);
  latex->DrawLatex(0.25, 0.85,Form("%s jets",chname[charged].Data()));
  latex->DrawLatex(0.25, 0.8,Form("R=%.1lf, fragmentation:%s",R, frag.Data()));
  latex->DrawLatex(0.25, 0.75,Form("parton distr. shape: %s",specshape[spec].Data()));
  
  TLegend *legratio = new TLegend(0.6, 0.65, 0.9, 0.90);
  legratio->SetFillStyle(0);
  legratio->SetBorderSize(0);
  legratio->SetTextSize(0.025);
  legratio->AddEntry(h4v0,"gen. jets, p_{T}^{lead}>0GeV", "lp");
  legratio->AddEntry(h4v5,"gen. jets, p_{T}^{lead}>5GeV", "lp");
  legratio->AddEntry(hjet0,"reco. jets, p_{T}^{lead}>0GeV", "lp");
  legratio->AddEntry(hjet5,"reco. jets, p_{T}^{lead}>5GeV", "lp");
  legratio->DrawClone("same");
  
  cspectra->cd(2);
  gPad->SetLogy();

  hjet5->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
  hjet5->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
  hproject57->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
  hproject7up->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
  
  hjet5->SetTitle("zoomed");
  hjet5->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  hjet5->DrawCopy("");
  hproject7up->DrawCopy("bsame");
  hproject57->DrawCopy("bsame");
  hproject7up->DrawCopy("same");
  h4v5->DrawCopy("same");
  hjet5->DrawCopy("same");
  
    TLegend *legratio = new TLegend(0.2, 0.75, 0.9, 0.90);
  legratio->SetFillStyle(0);
  legratio->SetBorderSize(0);
  legratio->SetTextSize(0.025);
  legratio->AddEntry(hproject57,"reco. jets origin. from gen. jets w/ p_{T}^{gen}<6GeV", "f");
  legratio->AddEntry(hproject7up,"reco. jets origin. from gen. jets w/ p_{T}^{gen}>6GeV", "f");
  legratio->DrawClone("same");

  TString str=Form("./obr/pyTest_%s_R0%.0lf_%s_%s.%s",chname[charged].Data(),R*10,frag.Data(),specshape[spec].Data(),ext.Data());
  cspectra->SaveAs(str.Data());

   
}

//============================================================================================
//============================================================================================


void DrawFF(float R=0.6, int charged=0, int spec=0, TString frag="g", bool pyEvent=0, TString ext="gif")
{
		TString chname[]={"full","charged",""}
	TString chname2=chname[charged];
	if(charged<2)chname2+="_";
	TString specshape[]={"flat","pT5"}
	TString finame;
	if(!pyEvent) finame=Form("root/pythia_%stest_%s_%.1lf_%s.root",chname2.Data(),specshape[spec].Data(),R,frag.Data());
	else  finame=Form("root/pythia_%spythia_%.1lf.root",chname2.Data(),R);
	TFile *fin=new TFile(finame,"OPEN");
	
	
	Color_t colorList[]={kRed,kBlue,kMagenta};
	

	TH1D* hnjets=(TH1D*) fin->Get("hnJets");
	TH2D* hzpT=(TH2D*) fin->Get("hz_pT");
	const int npTsl=2;
	TH1D* hz[npTsl];
	float zpTproject[npTsl]={45.0,15.0};
	
	for (int pr=0; pr<npTsl; pr++)
	{
		TString hnm=Form("hz%i",pr);
		int bin=hzpT->GetYaxis()->FindBin(zpTproject[pr]);
		hz[pr]=(TH1D)* hzpT->ProjectionX(hnm,bin,bin);
		hz[pr]->SetLineColor(colorList[pr]);
		int njets=hnjets->GetBinContent(hnjets->FindBin(zpTproject[pr]));
		hz[pr]->Scale(200.0/njets);
	}
	
	int can_x=600;
	int can_y=600;
	
	 TCanvas *cspec = new TCanvas("cspec","cspec",10,10,1*can_x,can_y);
  cspec->SetLogy();
  cspec->SetLogx();
  TString draw_str[3]={"","same","same"};
  hz[0]->GetXaxis()->SetRangeUser(0.1, 1.1);
  hz[0]->GetYaxis()->SetRangeUser(0.001, 35.0);
  hz[0]->SetTitle("FF");
  hz[0]->GetXaxis()->SetTitle("x");
  hz[0]->GetYaxis()->SetTitle("xD_{i}(x,Q^{2})");

  for( int i=0; i<npTsl; i++)
  {
	  hz[i]->Draw(draw_str[i]);
  }
  
  
  TLatex *latex2 = new TLatex();
  latex2->SetNDC();
  latex2->SetTextSize(0.04);
  if (!pyEvent) latex2->DrawLatex(0.4, 0.85,Form("%s jets, %s frag.",chname[charged].Data(),frag.Data()));
  else latex2->DrawLatex(0.3, 0.85,Form("PYTHIA %s jets, R=%.1lf",chname[charged].Data(),R));

  TLegend *legratio = new TLegend(0.7, 0.75, 0.9, 0.90);
  legratio->SetFillStyle(0);
  legratio->SetBorderSize(0);
  legratio->SetTextSize(0.03);
  for( int i=0; i<npTsl; i++)
  {
	  legratio->AddEntry(hz[i],Form("p_{T,jet}=%.0lf GeV",zpTproject[i]), "l");
  }
  legratio->DrawClone("same");


  TString str=Form("./obr/FF_%s_%s_pyEv%i.%s",frag.Data(),chname[charged].Data(),pyEvent,ext.Data());
  cspec->SaveAs(str.Data());
  
}

//============================================================================================
//============================================================================================


Draw_pTleadbias(float R=0.6, int charged=0, TString frag="u", TString ext="gif")
{
	
	TString chname[]={"full","charged",""}
	TString chname2=chname[charged];
	TString specshape[]={"flat","pT5"}
	
	
	const short ncompare=3;
	const float xmax=50.0;
	
	Color_t colorList[]={kRed,kBlue,kMagenta};
	
	TFile *fin[ncompare];
	TH1D* hjet0[ncompare];
	TH1D* hjet5[ncompare];
	TH1D* hratio[ncompare];
	
	for(int i=0;i<ncompare;i++)
	{
		TString finame=Form("root/pythia_%s_test_%s_%.1lf_%s.root",chname2.Data(),specshape[i].Data(),R,frag.Data());
		TString jet0name="hpTptlJet_pTl0";
		TString jet5name="hpTptlJet_pTl5";
		if(i==2) 
		{
			finame=Form("root/histos_pythiajet_%s_R%.1lf.root",chname2.Data(),R);
			jet0name="hpT_pTl0";
			jet5name="hpT_pTl5";
		}
		
		fin[i]=new TFile(finame,"OPEN");
		hjet0[i]=(TH1D*) fin[i]->Get(jet0name);
		hjet5[i]=(TH1D*) fin[i]->Get(jet5name);

		hratio[i]=(TH1D*) hjet5[i]->Clone(Form("hratio_%i",i));
		hratio[i]->Divide(hjet0[i]);
		hratio[i]->SetLineColor(colorList[i]);
		
	}
	

	hratio[0]->SetTitle("");
	hratio[0]->GetXaxis()->SetTitle("p_{T}^{jet}");
	hratio[0]->GetYaxis()->SetTitle("#frac{p_{T}^{lead}>5 GeV/c}{p_{T}^{lead}>0 GeV/c}");
	hratio[0]->GetXaxis()->SetTitleSize(0.04);
	hratio[0]->GetXaxis()->SetTitleOffset(1.2);
	hratio[0]->GetYaxis()->SetTitleOffset(1.4);
	hratio[0]->GetXaxis()->SetRangeUser(0, xmax);
	hratio[0]->GetYaxis()->SetRangeUser(0.005, 5.0);
	hratio[2]->Rebin(8);
	hratio[2]->Scale(1.0/8);
	
	int can_x=600;
	int can_y=600;
	TCanvas *cspec = new TCanvas("cspec","cspec",10,10,1*can_x,can_y);
	cspec->SetLogy();
	TString drawstyle="";
	for(int i=0;i<ncompare;i++)
	{
		if(i>0)drawstyle="same";
		hratio[i]->Draw(drawstyle);
	}
	
	TLatex *latex2 = new TLatex();
	latex2->SetNDC();
	latex2->SetTextSize(0.04);
	latex2->DrawLatex(0.2, 0.85,Form("%s jets, R=%.1lf",chname[charged].Data(),R));
	
	TLegend *legratio = new TLegend(0.5, 0.75, 0.9, 0.90);
	legratio->SetFillStyle(0);
	legratio->SetBorderSize(0);
	legratio->SetTextSize(0.035);
	legratio->AddEntry(hratio[0],Form("%s frag., flat spec.",frag.Data()), "l");
	legratio->AddEntry(hratio[1],Form("%s frag., p_{T}^{-5}",frag.Data()), "l");
	legratio->AddEntry(hratio[2],"PYTHIA jets", "l");
	legratio->DrawClone("same");
  
	TLine *one = new TLine(0, 1, xmax, 1);
	one->SetLineWidth(2);
	one->SetLineStyle(2);
	one->SetLineColor(kBlack);
	one->DrawClone("same");
  
	TString str=Form("./obr/pTleadBias_%s_%s_R0%.0lf.%s",frag.Data(),chname[charged].Data(),R*10,ext.Data());
	cspec->SaveAs(str.Data());
  
}
