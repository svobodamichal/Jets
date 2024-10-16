void main(float stat=1000)
{
	float bins[]={1,2,3,4,5,6,7};
	float bins2[]={1,2,4,7};
	
	TH1D* h1=new TH1D("h1","eq_bins",6,bins);
	TH1D* h2=new TH1D("h2","var_bins",3,bins2);
	TH1D* h3=new TH1D("h3","var_bins_weighted",3,bins2);
	
	TH1D* ha=new TH1D("ha","eq_bins",6,bins);
	TH1D* hb=new TH1D("hb","eq_bins",6,bins);
	TH1D* hc=new TH1D("hc","eq_bins",6,bins);
	
	for(int i=0;i<stat;i++)
	{
		float rnd=gRandom->Uniform(1,7);
		float weight=1.0/(rnd);
		h1->Fill(rnd,weight);
		h2->Fill(rnd,weight);
		h3->Fill(rnd,weight);
	}
	h3->Scale(1.0,"width");
	
	for(int i=0;i<stat/3.0;i++)
	{
		float rnd1=h1->GetRandom();
		float rnd2=h2->GetRandom();
		float rnd3=h3->GetRandom();
		
		ha->Fill(rnd1);
		hb->Fill(rnd2);
		hc->Fill(rnd3);
	}
	
	float can_x=400;
	float can_y=400;
	
	TCanvas *cspectra = new TCanvas("cspectra","cspectra",10,10,3*can_x,can_y);
  cspectra->Divide(3);
	
  h1->SetLineColor(kRed);
  h2->SetLineColor(kRed);
  h3->SetLineColor(kRed);
  
  float spectraYmin=1;
  float spectraYmax=stat/7.0;
  h1->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
  h2->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
  h3->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
  
  cspectra->cd(1);
  h1->Draw("");
  ha->Draw("same");
  
    cspectra->cd(2);
  h2->Draw("");
  hb->Draw("same");
  
    cspectra->cd(3);
  h3->Draw("");
  hc->Draw("same");
  
}
