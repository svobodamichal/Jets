Double_t levy(Double_t *x, Double_t *par)
{
	Double_t pT=x[0];
	Double_t mu=par[0];
	Double_t c=par[1];
	Double_t pwr=par[2];
	Double_t s=par[3];
	Double_t x0=par[4];
	Double_t pi=TMath::Pi();
	Double_t y=s*TMath::Sqrt((c/2*pi))*(TMath::Exp(-(c/(2*(pT-mu))))/TMath::Power((pT-mu),pwr));
	y=y*(1.0/(1.0+TMath::Exp(-(pT-x0))));
	return y;
}

Double_t LevyFitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for d2N/(2pi*pT dpT dy)
    // taken from here: http://sampa.if.usp.br/~suaide/blog/files/papers/PLB6372006.pdf
    Double_t pT, y, B, T, n, m0, mu,x0;
    B    = par[0];
    T    = par[1];
    n    = par[2];
    m0   = par[3];
	 mu   = par[4];
	 x0=par[5];
	 pT   = x_val[0];
    Double_t mT = TMath::Sqrt((pT-mu)*(pT-mu)+m0*m0);
    y = B/TMath::Power(1.0+(mT-m0)/(n*T),n);
	 y=y*(1.0/(1.0+TMath::Exp(-(pT-x0))));
    return y;
}

Double_t plaw(Double_t *x, Double_t *par)
{
	Double_t pT=x[0];
	Double_t pow=par[0];
	Double_t A=par[1];
	Double_t x0=par[2];
	Double_t y=A/TMath::Power(pT,pow);
	y=y*TMath::Power((1.0/(1.0+TMath::Exp(-(pT-x0)))),(10.0+pow)/10.0);
	return y;
}

void draw()
   {
		
	TString input="jetonly.root";
	TFile *f =new TFile(input);
	f->cd(); 
   TH1D* h1=f->Get("hDirectSpectrum");
		
	TH1 *frame = new TH1I("frame", "", 1000, 0, +100);
	TCanvas *cspectra = new TCanvas("cspectra","cspectra",10,10,1000,600);
	cspectra->cd();
	cspectra->SetGrid();
	cspectra->SetLogy();
	frame->GetXaxis()->SetRangeUser(0, 50);
	frame->GetYaxis()->SetRangeUser(1E-8, 1E2);
		
	   TF1 *myLevy = new TF1("myLevy",levy,0,50,5);
      myLevy->SetParameters(0,5,5,10,4);
      myLevy->SetParNames("mu","c","power","amplitude","x0");
		myLevy->SetNpx(1000);
      myLevy->Draw();
	
      TF1 *f2 = new TF1("pl",plaw,0,50,3);
		f2->SetParameters(5,10,4);
      f2->SetParNames("power","amplitude","x0");
		f2->SetLineColor(kBlue);
      f2->Draw("same");
		//f2->Draw();
		
		
		
		h1->SetLineColor(kBlack);
		h1->Scale(1.0/1000000);
		h1->Draw("same");
		
		TF1* LevyFit_pT = new TF1("LevyFit_pT",LevyFitFunc,0.0,50.0,6);
		LevyFit_pT->SetParameter(0,10); // B, changes the amplitude, 0.1
		LevyFit_pT->SetParameter(1,0.4); // T, changes slope, 0.4
		LevyFit_pT->SetParameter(2,5.8); // n, changes how fast spectrum drops, 5.8
		LevyFit_pT->SetParameter(3,0.0001); // m0, changes the width, 0.0001		
		LevyFit_pT->SetParameter(4,0); // mu, changes the x-axis shift
		LevyFit_pT->SetParameter(5,4);//x0, position of the low-pT knee
		LevyFit_pT->SetRange(0,50);
		LevyFit_pT->SetLineWidth(2);
		LevyFit_pT->SetLineStyle(1);
		LevyFit_pT->SetLineColor(kGreen);
		LevyFit_pT->SetNpx(1000);
		LevyFit_pT->DrawCopy("same");
		
   }