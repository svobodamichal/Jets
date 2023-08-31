bool skipiter(Int_t iter)
{
  Int_t it=iter+1;
	  if(it==3 || it==3|| it==4 || it==5 ||it==11 ||it==15)
				 return kFALSE;
	  else 
				 return kTRUE;
}

void show_unfold(Float_t pTthresh=0.0, Float_t R=0.4, Int_t priorNo=0, Int_t nbins=200)
{
  
 
  //*******************
  //*     SETUP       *
  //*******************
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetOptDate(1);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.09);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleSize(0.055,"Y");
  gStyle->SetTitleOffset(0.95,"Y");
  gStyle->SetTitleSize(0.06,"X");
  gStyle->SetTitleOffset(0.95,"X");
  gStyle->SetLabelSize(0.03,"X");
  gStyle->SetLabelSize(0.03,"Y");
  
  //canvas size
  Int_t can_x=1200; //1600
  Int_t can_y=680; //900
  
  
  //unfolding iterations
  Int_t firstiter=1;//first iteration (starting from 1!)
  Int_t lastiter=10; //last iteration (starting from 1!)
    //if(SVD1)lastiter=1;
  Int_t firstiter2=1;//first iteration for 2nd dataset
  Int_t lastiter2=10; //last iteration for 2nd dataset
    //if(SVD2)lastiter2=1;
  //Int_t bfIter=3; //which iteration to backfold and to show in comparison histograms (first=1, not 0)
  Int_t bfIter=5; //which iteration to backfold and to show in comparison histograms (first=1, not 0)
  //Int_t bfIter2=2;
  Int_t bfIter2=5;
  
  int statistics=10000; //number of hits in backfolding
  
  TString prior_type[]={/*0*/"pp_scaled",/*1*/"flat",/*2*/"pythia",/*3*/"powlaw5",/*4*/"powlaw6",/*5 */"powlaw4",/*6*/"powlaw3"};
  TString prior_type_name[]={"truth", "flat","biased Pythia","1/pT^{5}","1/pT^{6}","1/pT^{4}","1/(pT)^{3}"};
    
  TString wrkdir = Form("./Unfolded_R%.1lf_Bayes_%ibins",R,nbins);
  TString ext="gif";
  
  //TString responsedir="./root/response_matrix";
  TString responsedir=wrkdir;
    TString outdir = Form("./obr/unfolding_R%.1lf",R);
 

  //TString wrkdir_data = Form("./root/%s",trigger.Data());
  TString wrkdir_true = wrkdir;
  //TString errtype = "Multinomial";
  cout<<"true dir: "<<wrkdir_true.Data()<<endl;
  
  /*set marker styles: 
	      open	full
   circle	24	20
   square	25	21
   triangle up  26	22
   triangle dwn 32	23
   diamond	27	33
   star 	30	29
   */
  Int_t markerU=21; //marker style for unfolding
  Int_t markerP=20; //marker style for prior
  Int_t markerM=29; //marker style for measured
  Int_t markerB=22; //marker style for backfolded
  Int_t markerU2=25; //marker style for unfolding2
  Int_t markerP2=24; //marker style for prior2
  Int_t markerM2=30; //marker style for measured2
  Int_t markerB2=26; //marker style for backfolded2
  
  Float_t marker_size=1.0;
  Float_t line_width=2;
  
  Float_t iscale=1000; //integral scaler - so it can fit the legend window
  TString iscale_str="1E-3";
  Float_t priorRange=100; //up to which value do we want to plot th(no efficiency reconstruction correction)e prior distribution?
  
  Float_t spectraXmin=-20;
  Float_t spectraXmax=100;
  Float_t spectraYmin=1E4;
  Float_t spectraYmax=1E6;
  
  bfIter=bfIter-1; //becasue arrays start from 0
  bfIter2=bfIter2-1;
  firstiter=firstiter-1; //(change start from 1 to 0)
  firstiter2=firstiter2-1; //(change start from 1 to 0)
  const Int_t niter = lastiter;//number of iterations
  const Int_t niter2 = lastiter2;//number of iterations for 2nd dataset
  //*****end of SETUP*****
  
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.04);
  Color_t colorList[]={kBlack,kRed,kGreen+3,kMagenta+2,kBlue+3,kOrange+2,kYellow+2,kBlue-4,kOrange-2,kGreen,kRed+2,7,kGreen+4,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,28};
  Color_t colorList2[30]={kRed,kBlue-4,kOrange,kGreen,13,14,kMagenta+2,kGreen+3,15,kBlue,kOrange+2,kYellow+2,kRed+2,kBlack-2,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};
  TH1 *frame = new TH1I("frame", "", 1000, -100, +100);
 

  // UNFOLDING RESULTS
  TString str = Form("%s/%s/unfolded_SignalSpectrum_R%.1lf_pTthresh%.1lf.root", wrkdir.Data(),prior_type[priorNo].Data(), R, pTthresh);
  //if(doToymodel) str = Form("%s/%s/unfolded_SignalSpectrum_pTthresh%.1lf.root", wrkdir.Data(),prior_type[priorNo].Data(), pTthresh);
  //if (SVD1) str = Form("%s/%s/unfolded_SVD_k%i_R%.1lf_pTthresh%.1lf_%i.root", wrkdir.Data(),prior_type[priorNo].Data(),kTerm, R, pTthresh, nSVDbins);
   TFile *funfolding = new TFile(str.Data(), "OPEN");


  //RESPONSE MATRIX
  /*str=Form("%s/response_matrix_deltapT_R%.1lf.root", responsedir.Data(),R);
  if(doToymodel)str=Form("%s/response_matrix_deltapT_R%.1lf_TOY.root", responsedir.Data(),R);
  TFile *frmatrix = new TFile(str.Data(), "OPEN");
  str=Form("%s/response_matrix_deltapT_R%.1lf.root", responsedir2.Data(),R);
  TFile *frmatrix2 = new TFile(str.Data(), "OPEN");*/

  TH1D *hunfolded[niter];
  TH1D* hratut[niter];
  
  TH2D *hcovariance[niter];
  TH2D *hcorrelation[niter];

  Double_t integral[niter]; //integral of measured distribution and unfolded distributions
  Double_t integralm; //integral of measured dist.
  Double_t integralt; //integral of true distribution
  Double_t int_bfold; //integral of backfolded distributions
  TDirectoryFile *dir;
  
  
  //******Measured distribution *********
  TH1D *hmeasured;
  hmeasured= (TH1D*)((TDirectoryFile*)funfolding->Get("input"))->Get("hmeasured");
  //else hmeasured =(TH1D*) funfolding->Get("hmeasured");
  hmeasured->Sumw2();
  hmeasured->SetMarkerStyle(markerM);
  hmeasured->SetLineColor(colorList[0]);
  hmeasured->SetLineWidth(line_width);
  hmeasured->SetMarkerColor(colorList[0]);
  hmeasured->SetMarkerSize(marker_size);


  
  //******Prior distribution******
  TH1D *hprtmp;
  hprtmp = (TH1D*)((TDirectoryFile*)funfolding->Get("input"))->Get("hprior");
  //else hprtmp=(TH1D*) funfolding->Get("hprior");
  TH1D *hprior = (TH1D*)hmeasured->Clone("hprior");
  hprior->Sumw2();
  hprior->Reset("MICE");
  for(Int_t bin = 1; bin <= hprtmp->GetNbinsX(); bin++)
    {
      Double_t pTtmp = hprtmp->GetBinCenter(bin);
      Double_t yield = hprtmp->GetBinContent(bin);
      hprior->Fill(pTtmp, yield);
    }
  delete hprtmp;
  for(Int_t bin = 1; bin <= hprior->GetNbinsX(); bin++)
    {
      Double_t yield = hprior->GetBinContent(bin);
      hprior->SetBinError(bin,TMath::Sqrt(yield));
    }
  hprior->SetLineColor(colorList[1]);
  hprior->SetLineWidth(line_width);
  delete hprtmp;
  hprior->Scale(hmeasured->Integral("width")/hprior->Integral("width"));
 
  
  
  //******unfolded spectra******
  for(Int_t iter = firstiter; iter < lastiter; iter++)
    {
      //if(skipiter(iter) && SVD1) continue;
      cout<<"iter"<<iter<<endl;
      str = Form("iter%d", iter);
      dir = (TDirectoryFile*)funfolding->Get(str.Data());
      hunfolded[iter] = (TH1D*)dir->Get("hunfolded");
      //hunfolded[iter]->SetLineColor(colorList[iter+3]);
      //hunfolded[iter]->SetMarkerColor(colorList[iter+3]);
      hunfolded[iter]->SetLineColor(colorList[iter+3]-firstiter);
      hunfolded[iter]->SetMarkerColor(colorList[iter+3]-firstiter);
      hunfolded[iter]->SetMarkerStyle(markerU);
      hunfolded[iter]->SetMarkerSize(marker_size);
      hunfolded[iter]->SetLineWidth(line_width);
    	str = "hcovariance";
	hcovariance[iter] = (TH2D*)dir->Get(str.Data());

      TH1D *htemp = (TH1D*)hunfolded[iter]->Clone("htemp");
      hunfolded[iter]->Reset("MICE");

      for(Int_t bin = 1; bin <= htemp->GetNbinsX(); bin++)
	{
	  Double_t yield = htemp->GetBinContent(bin);
	  Double_t error;
	  //error = TMath::Sqrt(hcovariance[iter]->GetBinContent(bin, bin));
	  error = TMath::Sqrt(yield);
	  hunfolded[iter]->SetBinContent(bin, yield);
	  hunfolded[iter]->SetBinError(bin, error);
	}
      delete htemp;
	
    }
    
  
  //*******Response Matrix*************
  
   //TH2D * hresponse=(TH2D*) frmatrix->Get("hResponse_1E9");
  TH2D * hresponse;
  str = "input";
  dir = (TDirectoryFile*)funfolding->Get(str.Data());
  hresponse= (TH2D*)dir->Get("hresponse");

  
  
  //******backfolded distribution******
   /*
  TH1D * hbfold[niter];
  for(Int_t iter = firstiter; iter < lastiter; iter++)
    {
  if(iter!=bfIter) continue;
  str = Form("hbfold_%d", iter);
  hbfold[iter]=(TH1D*)hunfolded[iter]->Clone(str.Data());
  //hbfold[iter]->SetLineColor(colorList[iter+4]);
  //hbfold[iter]->SetMarkerColor(colorList[iter+4]);
  //if(iter==bfIter){
    hbfold[iter]->SetLineColor(colorList[2]);
    hbfold[iter]->SetLineWidth(line_width);
    hbfold[iter]->SetMarkerColor(colorList[2]);
  //}
  hbfold[iter]->SetMarkerStyle(markerB);
  hbfold[iter]->SetMarkerSize(marker_size);
  hbfold[iter]->Sumw2();
  hbfold[iter]->Reset("MICE");
  cout<<"creating backfolded spectrum:"<<endl;
  //for(Int_t bin = 1; bin <= hunfolded[iter]->GetNbinsX(); bin++)
  
  for(Int_t bin = 1; bin <= hunfolded[iter]->GetNbinsX(); bin++)
    {
      Double_t yield = hunfolded[iter]->GetBinContent(bin);      
      if(yield==0)continue;
      yield=yield/statistics;
      Double_t pT=hunfolded[iter]->GetBinCenter(bin);
      if(pT>pTcutoff)continue;
      Int_t xbin=hresponse->GetYaxis()->FindBin(pT);
      Double_t rnd=gRandom->Uniform(0,1);
     // if(rnd>0.5+binCenter) xbin=xbin-1;
      TH1D *hsmear = hresponse->ProjectionX("hsmear", xbin, xbin);
      Double_t norm=hsmear->Integral("width");
      hsmear->Scale(1/norm);
      for(int i=0;i<statistics;i++){
      Double_t pTsmear = hsmear->GetRandom();
      hbfold[iter]->Fill(pTsmear, yield);
      }
      delete hsmear;
      if(bin%10==0)cout<<"."<<endl;
    }//loop over bins  
    cout<<"done"<<endl;
  int_bfold=hbfold[iter]->Integral("width"); 
  }//loop over iterations

  
  double intbn=hbfold[bfIter]->Integral(1,hbfold[bfIter]->FindBin(pTthresh));
  double intmn=hmeasured->Integral(1,hmeasured->FindBin(pTthresh));
  
  double ratio=(intmn-intbn)/hmeasured->Integral();
  cout<<"ratio:"<<ratio<<endl;
  */
  

  
  //*****************************************
  // DRAWING SPECTRA
  //*****************************************


  TString suffix=Form("R%.1lf_%sPrior_pTthresh%.1lf_bfIter%i",R, prior_type[priorNo].Data(), pTthresh, bfIter+1);
   
  //MEASURED, PRIOR, UNFOLDED, BACKFOLDED DISTRIBUTIONS
  TCanvas *spectra = new TCanvas("cspectra","cspectra",10,10,can_x,can_y);
  cspectra->cd();
  cspectra->SetGrid();
  cspectra->SetLogy();
  frame->GetXaxis()->SetTitle("p_{T, corr}^{charged} (GeV/c)");
  frame->GetYaxis()->SetTitle("1/N_{events} 1/2#pi d^{2}N/dp_{T}^{ch}d#eta (GeV/c)^{-1}");
  frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
  frame->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
  //frame->GetYaxis()->SetRangeUser(1E-3, 1E5);
  str=Form("p_{T}^{leading} > %.1lf GeV/c, prior: %s",pTthresh, prior_type_name[priorNo].Data());
  frame->SetTitle(str);
  frame->DrawCopy("");
  hmeasured->DrawCopy("E same");
  hprior->GetXaxis()->SetRangeUser(pTthresh,priorRange);
  hprior->DrawCopy("HIST SAME");

  for(Int_t iter = firstiter; iter < lastiter; iter++){
    if(skipiter(iter))continue;
    hunfolded[iter]->DrawCopy("E same");
    }


  
   //Legend
  TLegend *legspectra = new TLegend(0.6467, 0.60, 0.89, 0.90);
  legspectra->SetTextSize(0.03);
  legspectra->SetFillStyle(0);
  legspectra->SetBorderSize(0);
  legspectra->AddEntry(hmeasured, "measured", "lp");
  legspectra->AddEntry(hprior,Form("prior: %s",prior_type_name[priorNo].Data()), "l");
  //if(showBF) legspectra->AddEntry(hbfold[bfIter],Form("backfolded dist. (iter%i)", bfIter+1), "lp");
  for(Int_t iter = firstiter; iter < lastiter; iter++){
    if(skipiter(iter))continue;
    legspectra->AddEntry(hunfolded[iter], Form("Unfolded, N_{iter} = %i", iter+1), "lp");
  }
  legspectra->DrawClone("same");
  
  latex->DrawLatex(0.35, 0.3,"STAR Preliminary");
  str = Form("%s/spectra_%s.%s", outdir.Data(),suffix.Data(),ext.Data());
  cspectra->SaveAs(str.Data());
  
  // DRAWING RATIO OF SUCCESIVE ITERATIONS
 
  TCanvas *cratio = new TCanvas("cratio","cratio",10,10,can_x,can_y);
  cratio->cd();
  cratio->SetGrid();
  //cratio->SetLogy();
  frame->GetXaxis()->SetRangeUser(0, 25);
  frame->GetYaxis()->SetRangeUser(0.5, 1.5); 
  frame->GetXaxis()->SetTitle("p_{T, corr}^{charged} (GeV/c)");
  //frame->GetYaxis()->SetTitle("Ratio = distribution /truth");
  frame->GetYaxis()->SetTitle("Ratio = n / (n-1) iteration");
  frame->DrawCopy("");
  for(Int_t iter = firstiter+1; iter < lastiter; iter++)
    {if(skipiter(iter))continue;
      //hunfolded[iter]->Divide(htruth);
      TH1D* hdiv=(TH1D*)hunfolded[iter]->Clone(Form("hdiv_%i",iter));
      hdiv->Divide(hunfolded[iter-1]);
      hdiv->DrawCopy("E same");
    }
  TLine *one = new TLine(0, 1, 25, 1);
  one->SetLineWidth(2);
  one->SetLineStyle(2);
  one->SetLineColor(kBlack);
  one->DrawClone("same");
  TLegend *legratio = new TLegend(0.6612, 0.7291, 0.8896, 0.90);
  legratio->SetFillStyle(0);
  legratio->SetBorderSize(0);
  for(Int_t iter = firstiter+1; iter < lastiter; iter++){
    if(skipiter(iter))continue;
    legratio->AddEntry(hunfolded[iter], Form("Unfolded, N_{iter}=%i/N_{iter}=%i", iter+1, iter), "lp");
    }
  legratio->DrawClone("same");
  latex->DrawLatex(0.4, 0.75,"STAR Preliminary");
  str = Form("%s/ratio_%s.%s", outdir.Data(),suffix.Data(),ext.Data());
  cratio->SaveAs(str.Data());
  
  
  
  
   
}//END
