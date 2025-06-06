
TH1D* rebinhisto(TH1D* hOLD, TH1D* hTEMPLATE, TString name)
{	
      Double_t yield[10];
      Double_t error[10];
      TH1D* hNEW=(TH1D*) hTEMPLATE->Clone(name.Data());
      hNEW->Reset("MICE");
      for(int bnvar=1;bnvar<hNEW->GetNbinsX()+1;bnvar++){
	int bincount=0;
	for(int bn=1; bn<hOLD->GetNbinsX()+1;bn++){
	  if(hNEW->FindBin(hOLD->GetBinCenter(bn))!=bnvar)continue;
	  //cout<<"bin new "<<bnvar<<"bin old "<<bn<<endl;
	  yield[bincount]=hOLD->GetBinContent(bn);
	  error[bincount]=hOLD->GetBinError(bn);
	  bincount++;
	}//loop over bins of original histogram
	
	double yieldt=0;
	double errt2=0;
	for(int i=0;i<bincount;i++)
	{
	  yieldt+=yield[i];
	  errt2+=error[i]*error[i];
	}
	hNEW->SetBinContent(bnvar,yieldt/bincount);
	hNEW->SetBinError(bnvar,TMath::Sqrt(errt2)/bincount);
      }//loop over bins of new histogram 
  return hNEW;
}


void plot_RAA(Float_t pTthresh=0.0, Float_t R=0.4, Int_t priorNo=0, TString ident="", Int_t nbins=200, TString trigger="MB")
{
  
 
  //*******************
  //*     SETUP       *
  //*******************
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetOptDate(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.09);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleSize(0.055,"Y");
  gStyle->SetTitleOffset(0.95,"Y");
  gStyle->SetTitleSize(0.06,"X");
  gStyle->SetTitleOffset(0.95,"X");
  gStyle->SetLabelSize(0.035,"X");
  gStyle->SetLabelSize(0.035,"Y");
  
  //canvas size
  Int_t can_x=600; //1600
  Int_t can_y=900; //900
  
  Int_t SVD=0; //first data set is SVD unfolding
  Int_t quad=0;
  Int_t semilog=1;
  
  const Int_t ninputs=5;
  Int_t compIter[ninputs]={4,4,4,6,4}; //to which iteration compare all other iterations
  
  TString prior_type[]={/*0*/"pp_scaled",/*1*/"flat",/*2*/"pythia",/*3*/"powlaw5",/*4*/"powlaw6",/*5 */"powlaw4",/*6*/"powlaw3"};
  TString prior_type_name[]={"truth", "flat","biased Pythia","1/pT^{5}","1/pT^{6}","1/pT^{4}","1/(pT)^{3}"};
    

  TString UnfType="BG_dete";
  TString desctype="BG+detector effects";
  TString desc=Form(", corrected for %s",desctype.Data());

   TString title="Charged jets";
   TString titlex="p_{T, jet}^{charged} (GeV/c)";
   TString titley="1/N_{evt} 1/(2#\pi) dN/dp_{T, jet}^{charged}d#\eta";

   TString ext="png"; //output file extension
  //TString comment=UnfType[unftype]; //will be attached to the end of the output filename and also detrmines the name of the subdirectory with figures
  TString comment=ident; //will be attached to the end of the output filename 
  TString comment2=ident;
  TString prefix="unf_charged"; //beginning of the output filenames
  TString unf_type="Bayes"; //will be in the title of histograms
  if(SVD==1) unf_type="SVD";


  TString wrkdir[ninputs];


     wrkdir[0]=Form("./root/%s/R%.1lf",trigger.Data(),R);
     wrkdir[1] = Form("./root/%s/R%.1lf",trigger.Data(),R);
     wrkdir[2] = Form("./root/%s/R%.1lf",trigger.Data(),R);
     wrkdir[3] = Form("./root/%s/R%.1lf",trigger.Data(),R);
     wrkdir[4]="root/toymodel";

  TString str;
  if(SVD)str="SVD";
  else str="Bayes";
  wrkdir[0] = Form("%s/Unfolded_R%.1lf_%s_%ibins_%s",wrkdir[0].Data(),R,str.Data(),nbins,UnfType.Data());
  wrkdir[1] = Form("%s/Unfolded_R%.1lf_%s_%ibins_%s_m5p",wrkdir[1].Data(),R,str.Data(),nbins,UnfType.Data());
  wrkdir[2] = Form("%s/Unfolded_R%.1lf_%s_%ibins_%s_p5p",wrkdir[2].Data(),R,str.Data(),nbins,UnfType.Data());
  wrkdir[3] = Form("%s/Unfolded_R%.1lf_%s_%ibins_%s",wrkdir[3].Data(),R,str.Data(),nbins,UnfType.Data());
  wrkdir[4] = Form("%s/Unfolded_R%.1lf_%s_%ibins",wrkdir[4].Data(),R,str.Data(),nbins);

  TString outdir ="./obr/RAA";
 
  TString wrkdir_true = "root/toymodel/jetonly";
  
  /*set marker styles: 
	      open	full
   circle	24	20
   square	25	21
   triangle up  26	22
   triangle dwn 32	23
   diamond	27	33
   star 	30	29
   */
  Int_t marker=29; //marker style for raa

  
  Float_t marker_size=1.5;
  Float_t line_width=2;
  
  Float_t iscale=1000; //integral scaler - so it can fit the legend window
  TString iscale_str="1E-3";
  Float_t priorRange=35; //up to which value do we want to plot the prior distribution?
  
  Float_t Xmin=0.0;
  Float_t Xmind=pTthresh;
  Float_t Xmaxd=30;
  Float_t Xmax=35;
  Float_t Ymin=1E-08; //Jana
  Float_t Ymax=1E-02; //Jana
  
  
  //*****end of SETUP*****

  const int nloops=3;
  double toteyl[50][nloops]; //total uncertainty
  double toteyh[50][nloops]; //total uncertainty
  
  
  //for(int loop=0; loop<nloops; loop++){
  
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.035);
  Color_t colorList[]={kBlack,kRed,kGreen+3,kMagenta+2,kBlue,kOrange+2,kYellow+2,kBlue-9,kOrange-2,kGreen,kRed+2,7,kGreen+4,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,28};
  Color_t colorList2[30]={kRed,kBlue-4,kOrange,kGreen,13,14,kMagenta+2,kGreen+3,15,kBlue,kOrange+2,kYellow+2,kRed+2,kBlack-2,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};
  TH1 *frame = new TH1I("frame", "", 1000, -100, +100);
 

  str = Form("%s/../histos_inclusivejet_R%.1lf.root", wrkdir[0].Data(),R);
  TFile *f = new TFile(str.Data(), "OPEN");
  TH1I* hevents;
  hevents= (TH1I*) f->Get("hevts");
  Int_t nevents=hevents->GetEntries();
  f->Close();
  
  str = Form("%s/jetonly.root", wrkdir_true.Data());
  f = new TFile(str.Data(), "OPEN");
  hevents= (TH1I*) f->Get("hevents");
  Int_t toy_events=hevents->GetEntries();
  f->Close();

  Double_t hole=(1/12)+(R*2)/(2*TMath::Pi()); //fraction of acceptance which was droped due to the bad sectors
  cout<<"hole: "<<hole<<endl;
  Double_t scale_jets = 1./(2*(1-R)*2.*TMath::Pi()*nevents*(1-hole));
  Double_t scale_toy=1./(2*(1-R)*2.*TMath::Pi()*toy_events);
  
  Double_t TAA=22;
  Double_t scale_true=TAA/(2*TMath::Pi());
  Double_t scale_true_toy=8*(TAA)/(2*TMath::Pi());
  //Double_t binWidth;

  // UNFOLDING RESULTS
  
    TFile *funfolding[ninputs];
for(int i=0;i<5;i++){
  str = Form("%s/%s/unfolded_SignalSpectrum_R%.1lf_pTthresh%.1lf.root", wrkdir[i].Data(),prior_type[priorNo].Data(), R, pTthresh);
  if (i==3 || i==4) str = Form("%s/%s/unfolded_SignalSpectrum_R%.1lf_pTthresh%.1lf.root", wrkdir[i].Data(),prior_type[6].Data(), R, pTthresh);
  if (SVD && i!=4) str = Form("%s/%s/unfolded_SVD_R%.1lf_pTthresh%.1lf.root", wrkdir[i].Data(),prior_type[priorNo].Data(), R, pTthresh);
  funfolding[i] = new TFile(str.Data(), "OPEN");
}
   
  TH1D *hunfolded[ninputs];
  TH2D *hcovariance[ninputs];
  TH1D *hunfoldedrb[ninputs];
  
  //VARIABLE BINNING
 
  Float_t pTrange=50; //pTrange in input histograms (has to be same in all histos!)
  const Int_t newbins=35;
  Double_t pTbinArray[newbins+1];
  pTbinArray[0]=0;
  Double_t width=1;
    cout<<"binning:"<<endl;

  for(int i=1; i<=newbins; i++){
   if(pTbinArray[i-1]==20) width=2;
   pTbinArray[i]=pTbinArray[i-1]+width;
  }
   
  cout<<"binning:"<<endl;
  for(int i=1; i<=newbins; i++){
  cout<<i<<": "<<pTbinArray[i-1]<<" - "<<pTbinArray[i]<<endl;}
    cout<<"binning3"<<endl;

  TH1D *hvarbin=new TH1D("hvarbin","hvarbin",newbins,pTbinArray);

  Double_t integral[ninputs]; //integral of measured distribution and unfolded distributions
  TDirectoryFile *dir;

  //******Measured distribution *********
  TH1D *hmeasured;
  hmeasured= (TH1D*)((TDirectoryFile*)funfolding[0]->Get("input"))->Get("hmeasured");
  //else hmeasured =(TH1D*) funfolding->Get("hmeasured");
  hmeasured->Sumw2();
  hmeasured->SetLineColor(colorList[0]);
  hmeasured->SetLineWidth(line_width);
  hmeasured->SetMarkerColor(colorList[0]);
  hmeasured->Scale(scale_jets,"width");
  float integralm=hmeasured->Integral("width");
  
  TH1D *hmeasuredtoy;
  hmeasuredtoy= (TH1D*)((TDirectoryFile*)funfolding[4]->Get("input"))->Get("hmeasured");
  hmeasuredtoy->Scale(scale_toy,"width");
  Double_t integralmt=hmeasuredtoy->Integral("width");
  
   //******toymodel true distribution******
  str = Form("%s/jetonly.root", wrkdir_true.Data());
  TFile *ftruth = new TFile(str.Data(), "OPEN");
  TH2D *hPtRecpTleadingPrior = (TH2D*)ftruth->Get("fhPtRecpTleading");
  Int_t firstb = hPtRecpTleadingPrior->GetXaxis()->FindBin(pTthresh);
  Int_t lastb = hPtRecpTleadingPrior->GetNbinsX();
  TH1D *htemp = hPtRecpTleadingPrior->ProjectionY("htemp", firstb, lastb);
  TH1D *htruth = (TH1D*)hmeasured->Clone("htruth");
  htruth->Sumw2();
  htruth->Reset("MICE");
  for(Int_t bin = 1; bin <= htemp->GetNbinsX(); bin++)
    {
      Double_t pTtmp = htemp->GetBinCenter(bin);
      Double_t yield = htemp->GetBinContent(bin);
      Double_t oldWidth=htemp->GetBinWidth(bin);
      Int_t newBin=htruth->FindBin(pTtmp);
      Double_t newWidth=htruth->GetBinWidth(newBin);
      //cout<<"bin:"<<binx<<" pT:"<<pT<<" yield:"<<yield<<" width: "<<binWidth<<" -> "<<newWidth<<" new yield:";
   //yield=yield*(oldWidth/newWidth);
      htruth->Fill(pTtmp, yield);
    }
  delete htemp;
  for(Int_t bin = 1; bin <= htruth->GetNbinsX(); bin++)
    {
      Double_t yield = htruth->GetBinContent(bin);
      htruth->SetBinError(bin,TMath::Sqrt(yield));
    }
  htruth->SetMarkerColor(kRed+1);
  htruth->SetLineColor(kRed+1);
  htruth->SetLineWidth(line_width);
  htruth->Scale(scale_toy,"width");
  //htruth->Scale(integralmt/htruth->Integral("width"));
  float integralt=htruth->Integral("width");

  //TAA*pythia
    str = "./root/MB/PythiaHistosChargedJets.root";
  ftruth = new TFile(str.Data(), "OPEN");
  str=Form("h_PythiaJetPt0p%.0lfTot",R*10);
  TH1D *htemp =(TH1D*) ftruth->Get(str.Data());
  
  TH1D *htruthp = (TH1D*)hmeasured->Clone("htruthp");
  htruthp->Sumw2();
  htruthp->Reset("MICE");
  for(Int_t bin = 1; bin < htemp->GetNbinsX()+1; bin++)
    {
      Double_t pTtmp = htemp->GetBinCenter(bin);
      Double_t yield = htemp->GetBinContent(bin);
      Double_t oldWidth=htemp->GetBinWidth(bin);
      Int_t newBin=htruthp->FindBin(pTtmp);
      Double_t newWidth=htruthp->GetBinWidth(newBin);
      //cout<<"bin:"<<binx<<" pT:"<<pT<<" yield:"<<yield<<" width: "<<binWidth<<" -> "<<newWidth<<" new yield:";
      yield=yield*(oldWidth/newWidth);
      htruthp->Fill(pTtmp, yield);
    }
  for(Int_t bin = 1; bin <= htruthp->GetNbinsX(); bin++)
    {
      Double_t yield = htruthp->GetBinContent(bin);
      Int_t bint=htemp->FindBin(htruthp->GetBinCenter(bin));
      Double_t error= (htemp->GetBinError(bint))/(htemp->GetBinContent(bint));
      htruthp->SetBinError(bin,error*yield);
    }
      delete htemp;

  htruthp->SetMarkerColor(kGray+1);
  htruthp->SetLineColor(kGray+1);
  htruthp->SetLineWidth(line_width);
  htruthp->Sumw2();
  TH1D* htruthp2=(TH1D*) htruthp->Clone("htruthp2");
  //htruthp->Scale(scale_true,"width"); //Jana
  //htruthp2->Scale(scale_true_toy,"width"); //Jana
  float integraltp=htruthp->Integral("width");
  
 
  //******unfolded spectra******
  for(Int_t iter = 0; iter < 5; iter++)
    {
      cout<<"iter"<<iter<<endl;
      if(!SVD)	str = Form("iter%d", compIter[iter]);
      else str = Form("kterm%d",  compIter[iter]);
      dir = (TDirectoryFile*)funfolding[iter]->Get(str.Data());
      hunfolded[iter] = (TH1D*)dir->Get("hunfolded");
      hunfolded[iter]->SetLineColor(colorList[iter+3]);
      hunfolded[iter]->SetMarkerColor(colorList[iter+3]);
      hunfolded[iter]->SetLineWidth(line_width);

	str = "hcovariance";
	if(iter!=4)
	hcovariance[iter] = (TH2D*)dir->Get(str.Data());

      TH1D *htemp = (TH1D*)hunfolded[iter]->Clone("htemp");
      hunfolded[iter]->Reset("MICE");

      for(Int_t bin = 1; bin <= htemp->GetNbinsX(); bin++)
	{
	  Double_t yield = htemp->GetBinContent(bin);
	  Double_t error;
	  if(iter!=4) error = TMath::Sqrt(hcovariance[iter]->GetBinContent(bin, bin));
	  else error = TMath::Sqrt(yield);
	  hunfolded[iter]->SetBinContent(bin, yield);
	  hunfolded[iter]->SetBinError(bin, error);
	}
      delete htemp;
	
      if(iter!=4)
	hunfolded[iter]->Scale(scale_jets,"width");
      else
	hunfolded[iter]->Scale(scale_toy,"width");
      integral[iter]=hunfolded[iter]->Integral("width");
      cout<<"integral unfolded:"<<integral[iter]<<endl;
      
      //creating new, rebinned histograms
      str=Form("hunfoldedrb_%i",iter);
      hunfoldedrb[iter]=(TH1D*) rebinhisto(hunfolded[iter],hvarbin,str);
    }//loop over input files
      
      str="htruthrb";
      TH1D* htruthrb=(TH1D* )rebinhisto(htruth,hvarbin,str);
      //htruthrb->Draw();
      
      str="htruthprb";
      TH1D* htruthprb=(TH1D*)rebinhisto(htruthp,hvarbin,str);
      
      str="htruthp2rb";
      TH1D* htruthp2rb=rebinhisto(htruthp2,hvarbin,str);
    
  //*****************************************
  // DRAWING SPECTRA
  //*****************************************

    
  TString suffix=Form("R%.1lf_%sPrior_pTthresh%.1lf",R, prior_type[priorNo].Data(), pTthresh); 
  suffix=Form("%s_%s",suffix.Data(),comment.Data()); 
  
  /*
  TCanvas *ct= new TCanvas("ct","ct",10,10,can_x,can_y); 
  htruth->Draw();
  //hmeasuredtoy->Draw("same");
  hunfolded[4]->Draw("same");
  htruthp2->Draw("same");
  */
  
  // RAA1 VS RAA2 VS RAA3

  TCanvas *cratiobmbm = new TCanvas("cratiobmbm","cratiobmbm",10,10,can_x,can_y);
  cratiobmbm->cd();
  cratiobmbm->SetGrid();
  // cratiobb->SetLogy(); 
  frame->GetXaxis()->SetRangeUser(Xmin, Xmax);
  frame->GetYaxis()->SetRangeUser(Ymin, Ymax);
  frame->GetXaxis()->SetTitle("p_{T, corr}^{charged} (GeV/c)");
  frame->GetYaxis()->SetTitle("R_{AA}");
  frame->SetTitle(Form("R_{AA} comparison, %s, p_{T}^{leading} > %.1lf GeV/c",desc.Data(),pTthresh));
  frame->DrawCopy("");
  TH1D* hraa[4];
  for(int i=0;i<4;i++){
  hraa[i]=(TH1D*)hunfoldedrb[i]->Clone(Form("hraa_%i",i));
  // hraa[i]->Divide(htruthprb); //Jana
  hraa[i]->DrawCopy("E same");
  }
  
  TH1D* hraa2=(TH1D*)hunfoldedrb[4]->Clone("hraa2");
  // hraa2->Divide(htruthp2rb); //Jana
  hraa2->DrawCopy("E same");
  
  TH1D* hraa3=(TH1D*)htruthrb->Clone("hraa3");
  // hraa3->Divide(htruthp2rb); //Jana
  hraa3->DrawCopy("E same");
  
  //Filling arrays with systematic uncertainties
  //--------------------------------------------
  int zerobin=hraa[0]->FindBin(Xmind);
  int lastbin=hraa[0]->FindBin(Xmaxd);
  const int Nbinsx=lastbin-zerobin;

  double ax[Nbinsx];
  double ay[Nbinsx];
  double aexl[Nbinsx];
  double aexh[Nbinsx];
  double aeyh[Nbinsx]; //efficiency uncertainty
  double aeyl[Nbinsx]; //efficiency uncertainty
  double beyl[Nbinsx]; //unfolding uncertainty
  double beyh[Nbinsx]; //unfolding uncertainty
  double ceyl[Nbinsx]; //prior uncertainty
  double ceyh[Nbinsx]; //prior uncertainty
  double cey[Nbinsx];


  
  
  double bey[Nbinsx];
  /*{0,0,0,0,0,0,0,0.013313,-0.0175193,-0.0271469,-0.0226856,-0.00175555,0.0157637,0.032898,
      0.0641068,0.0805159,0.0620794,0.0519752,0.0637049,0.0134455,0.031961,-0.05556,0.00937,-0.112936,-0.30744,
      -0.112863,-0.0408639,-0.236758,-0.184855,-0.2009,-0.282365,-0.240261,-0.554647,-0.148695,-0.138247,0.444211,
      0.0237264,0.303296,0.46687,0.415405,0.599224};  //relative, wigles
    */  
   cout<<"wigle error"<<endl;
  for(int binx=zerobin; binx<lastbin;binx++ )
  {
    Int_t newbin=binx-zerobin;
    Double_t hval=hraa[0]->GetBinContent(binx);
    Double_t pT=hraa[0]->GetBinCenter(binx);
    ax[newbin]=hraa[0]->GetBinCenter(binx);
    ay[newbin]=hraa[0]->GetBinContent(binx);
    aexl[newbin]=hraa[0]->GetBinWidth(binx)/2;
    aexh[newbin]=hraa[0]->GetBinWidth(binx)/2;
    
    
    aeyl[newbin]=TMath::Abs(hraa[0]->GetBinContent(binx) - hraa[2]->GetBinContent(binx));
    aeyh[newbin]=TMath::Abs(hraa[0]->GetBinContent(binx) - hraa[1]->GetBinContent(binx));

    cey[newbin]=hraa[3]->GetBinContent(binx)- hraa[0]->GetBinContent(binx);

    
  
    double maxer1;
    double start1;
    double center1;
    double center1b;
    double end1;
    double maxer2;
    double start2;
    double center2;
    double center2b;
    double end2;
    
if(R>0.35)//R=0.4
    {
      maxer1=0.1;
      start1=15;
      center1=20;
      center1b=22;
      end1=30;
      maxer2=0.20;
      start2=25;
      center2=35;
      center2b=35;
      end2=50;
    }
    else if(R>0.25){//R=0.3
      if(pTthresh<6) //pTlead=5
      {  
	maxer1=0.05;
	start1=8;
	center1=14;
	center1b=18;
	end1=30;
	maxer2=0.20;
	start2=20;
	center2=35;
	center2b=35;
	end2=50;
      }
      else //pTlead=7
      {
	maxer1=0.05;
	start1=10;
	center1=18;
	center1b=18;
	end1=30;
	maxer2=0.20;
	start2=20;
	center2=27;
	center2b=30;
	end2=50;
      }
    }
    else{ //R=0.2
      if(pTthresh<6) //pTlead=5
	{
        maxer1=0.05;
	start1=8;
	center1=14;
	center1b=18;
	end1=30;
	maxer2=0.20;
	start2=12;
	center2=25;
	center2b=35;
	end2=50;
	}
	else
	{
	maxer1=0.05;
	start1=10;
	center1=18;
	center1b=18;
	end1=30;
	maxer2=0.15;
	start2=12;
	center2=30;
	center2b=35;
	end2=50;
	}
    }
    
    maxer2=maxer2*hval;
    maxer1=maxer1*hval;
    
      beyl[newbin]=0;
      beyh[newbin]=0;
      //toteyl[newbin][loop]=0;
      //toteyh[newbin][loop]=0;
      
    if(pT>start1 && pT<=center1)
      beyl[newbin]=maxer1*(pT-start1)/(center1-start1);
    else if(pT>center1b && pT<end1)
      beyl[newbin]=maxer1*(end1-pT)/(end1-center1b);
    else if(pT>center1 && pT<=center1b)
      beyl[newbin]=maxer1;
    
    
    if(pT>start2 && pT<=center2)
      beyh[newbin]=maxer2*(pT-start2)/(center2-start2);
    else if(pT>center2b && pT<end2)
      beyh[newbin]=maxer2*(end2-pT)/(end2-center2b);
   else if(pT>center2 && pT<=center2b)
      beyh[newbin]=maxer2;
   
   //toteyl[newbin][loop]+=beyl[newbin]*beyl[newbin];
   //toteyh[newbin][loop]+=beyh[newbin]*beyh[newbin];

   //beyl[newbin]=beyl[newbin]+aeyl[newbin];
   //beyh[newbin]=beyh[newbin]+aeyh[newbin];

      ceyl[newbin]=0;
      ceyh[newbin]=0;
      
    if(cey[newbin]<0){
      ceyl[newbin]=TMath::Abs(cey[newbin]);
      }
    else{
      ceyh[newbin]=TMath::Abs(cey[newbin]);
      }
      
    //toteyl[newbin][loop]+=ceyl[newbin]*ceyl[newbin];
    //toteyh[newbin][loop]+=ceyh[newbin]*ceyh[newbin];
    
    ceyl[newbin]=TMath::Sqrt(ceyl[newbin]*ceyl[newbin]+beyl[newbin]*beyl[newbin])+aeyl[newbin];
    ceyh[newbin]=TMath::Sqrt(ceyh[newbin]*ceyh[newbin]+beyh[newbin]*beyh[newbin])+aeyh[newbin];

    //toteyl[newbin][loop]=TMath::Sqrt(toteyl[newbin][loop]);
    //toteyh[newbin][loop]=TMath::Sqrt(toteyh[newbin][loop]);
  }//loop over bins

  
  //Legend
  TLegend *leg = new TLegend(0.6467, 0.75, 0.89, 0.90);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hraa[0], "R_{AA}", "lp");
  leg->AddEntry(hraa[1], "R_{AA} #epsilon-5\%", "lp");
  leg->AddEntry(hraa[2], "R_{AA} #epsilon+5\%", "lp");
  leg->AddEntry(hraa[3], "R_{AA} 2nd prior", "lp");
  leg->AddEntry(hraa2,"R_{AA} unf/pp", "lp");
  leg->AddEntry(hraa3,"R_{AA} true/pp", "lp");
  leg->DrawClone("same");

  TLine *one = new TLine(0, 1, 40, 1);
  one->SetLineWidth(2);
  one->SetLineStyle(2);
  one->SetLineColor(kBlack);
  one->DrawClone("same");

  //}//loop over R or pTlead
  
  //=========================================
  
   // RAA
  TCanvas *craa = new TCanvas("craa","craa",10,10,can_x,can_y);
  craa->cd();
  //craa->SetGrid();
  
  if(semilog)craa->SetLogy();
  frame->GetXaxis()->SetRangeUser(0, 35);
  frame->GetYaxis()->SetRangeUser(0, 1.5); 
  frame->GetXaxis()->SetTitle("");
  //frame->GetYaxis()->SetTitle("Ratio = distribution /truth");
  //frame->GetYaxis()->SetTitle("Ratio = #frac{unfolded}{T_{AA}*d#sigma_{pp}/dp_{T}}");
  //frame->GetYaxis()->SetTitle("Ratio = #frac{unfolded}{true}");
  frame->GetYaxis()->SetTitle("");
  frame->Draw("");
  
  craa->SetLogy(); 
  TH1D* hrat=(TH1D*)hunfoldedrb[0]->Clone("hrat");
  // hrat->Divide(htruthprb);
      hrat->SetLineColor(kBlue);
      hrat->SetMarkerStyle(marker);
      hrat->SetMarkerColor(kBlue);
      hrat->SetMarkerSize(marker_size);
      hrat->SetLineWidth(2);
      
   TGraphAsymmErrors* geffi = new TGraphAsymmErrors(Nbinsx, ax, ay, aexl, aexh, aeyl, aeyh);
   geffi->SetFillColor(kOrange+1);
   geffi->SetFillStyle(1001);
   geffi->SetName("geffi");
   
   TGraphAsymmErrors* gunfolding = new TGraphAsymmErrors(Nbinsx, ax, ay, aexl, aexh, ceyl, ceyh);
   gunfolding->SetFillColor(kGreen-4);
   gunfolding->SetFillStyle(1001);
   gunfolding->SetName("gunfolding");
   
   TGraphAsymmErrors* gprior = new TGraphAsymmErrors(Nbinsx, ax, ay, aexl, aexh, ceyl, ceyh);
   gprior->SetFillColor(kAzure-9);
   gprior->SetFillStyle(1001);
   gprior->SetName("gprior");
   
   double taa_x[1]={Xmax-1};
   double taa_y[1]={1};
   double taa_exl[1]={0.5};
   double taa_exh[1]={0.5};
   double taa_eyl[1]={0.05};
   double taa_eyh[1]={0.05};
   
   TGraphAsymmErrors* gtaa = new TGraphAsymmErrors(1, taa_x, taa_y, taa_exl, taa_exh, taa_eyl, taa_eyh);
   gtaa->SetFillColor(kRed+1);
   gtaa->SetFillStyle(3001);
   gtaa->SetName("gtaa");
   
   
   TAxis *axis = geffi->GetXaxis();
   axis->SetLimits(Xmin,Xmax);                 // along X
   geffi->GetHistogram()->SetMaximum(Ymax);   // along          
   geffi->GetHistogram()->SetMinimum(Ymin);  //   Y 
   /*
    TAxis *cxis = gprior->GetXaxis();
   cxis->SetLimits(Xmin,Xmax);                 // along X
   gprior->GetHistogram()->SetMaximum(Ymax);   // along          
   gprior->GetHistogram()->SetMinimum(Ymin);
   */
    TAxis *bxis = gunfolding->GetXaxis();
   bxis->SetLimits(Xmin,Xmax);                  // along X
   gunfolding->GetHistogram()->SetMaximum(Ymax);   // along          
   gunfolding->GetHistogram()->SetMinimum(Ymin);
 
    TAxis *dxis = gtaa->GetXaxis();
   dxis->SetLimits(Xmin,Xmax);                  // along X
   gtaa->GetHistogram()->SetMaximum(Ymax);   // along          
   gtaa->GetHistogram()->SetMinimum(Ymin);
   
   geffi->SetTitle(title.Data());
   gunfolding->SetTitle(title.Data());
   //gprior->SetTitle(title.Data());
   gtaa->SetTitle(title.Data());
   
   geffi->GetXaxis()->SetTitle(titlex.Data());
   gunfolding->GetXaxis()->SetTitle(titlex.Data());
   //gprior->GetXaxis()->SetTitle(titlex.Data());
   gtaa->GetXaxis()->SetTitle(titlex.Data());

   geffi->GetYaxis()->SetTitle(titley.Data());
   gunfolding->GetYaxis()->SetTitle(titley.Data());
  // gprior->GetYaxis()->SetTitle(titley.Data());
   gtaa->GetYaxis()->SetTitle(titley.Data());
      
   /*gprior->Draw("a2");
   
   craa->cd();
   TPad *overlay = new TPad("overlay","",0,0,1,1);
   overlay->SetFillStyle(4000);
   overlay->SetFillColor(0);
   overlay->SetFrameFillStyle(4000);
   overlay->Draw();
   overlay->cd();
   if(semilog)overlay->SetLogy();*/
   gunfolding->Draw("a2");
   
   craa->cd();
   TPad *overlay2 = new TPad("overlay","",0,0,1,1);
   overlay2->SetFillStyle(4000);
   overlay2->SetFillColor(0);
   overlay2->SetFrameFillStyle(4000);
    overlay2->Draw();
    overlay2->cd();
      if(semilog)overlay2->SetLogy();

   geffi->Draw("a2");
   //geffi->Draw("p");
   
   /*craa->cd();
   TPad *overlay3 = new TPad("overlay","",0,0,1,1);
   overlay3->SetFillStyle(4000);
   overlay3->SetFillColor(0);
   overlay3->SetFrameFillStyle(4000);
    overlay3->Draw();
   overlay3->cd();
   if(semilog)overlay3->SetLogy();

   gtaa->Draw("a2");*/
   
   hrat->SetAxisRange(Xmind,Xmaxd-0.001,"X");
   hrat->DrawCopy("esame");
  
  TLine *one = new TLine(Xmin, 1, Xmax, 1);
  one->SetLineWidth(2);
  one->SetLineStyle(2);
  one->SetLineColor(kGray+3);
  //one->DrawClone("same");
  
  
  double posx1=0.62;
  double posy1=0.64;
  double posx2=0.88;
  double posy2=0.88;
  
  if(R>0.35){
    posy1=0.2;
    posy2=0.4;
    posx1=0.55;
    
  }
  // TLegend *legraa = new TLegend(posx1,posy1,posx2,posy2);
  TLegend *legraa = new TLegend(0.15,0.18,0.45,0.38);
  legraa->SetTextSize(0.03);
  legraa->SetFillStyle(0);
  legraa->SetBorderSize(0);
  //legraa->AddEntry(hrat, "AuAu/Pythia","lp"); //Jana
    legraa->AddEntry(geffi, "tracking efficiency uncertainty", "f");
      legraa->AddEntry(gunfolding, "unfolding uncertainty", "f");
    //legraa->AddEntry(gprior, "prior choice uncertainty", "f");
    //  legraa->AddEntry(gtaa, "T_{AA} uncertainty", "f");
    legraa->AddEntry("","Uncertainties added linearly","");
 
 
  legraa->DrawClone("same");
  str = Form("%s/%s_spectrum_%s.%s", outdir.Data(),prefix.Data(),suffix.Data(),ext.Data());
  
    //Info table
  TLegend *model_info = new TLegend(0.35, 0.70, 0.55, 0.87);
  model_info->SetTextSize(0.035);
  model_info->SetFillStyle(0);
  model_info->SetBorderSize(0);
  model_info->SetMargin(0.05);
  model_info->SetHeader("Run 11 Au+Au #sqrt{s_{NN}}=200 GeV, 60 #mub^{-1}");
  model_info->AddEntry("", "0-10% Central Collisions", "");
  model_info->AddEntry("", Form("Anti-k_{T} \t \t R = %.1lf",R), "");
  model_info->AddEntry("", "p_{T}^{const} > 0.2 GeV/c", "");
  model_info->AddEntry("", Form("p_{T}^{leading} > %.1lf GeV/c",pTthresh), "");
   
  if(R>0.35) model_info->AddEntry("", "A_{reco jet} > 0.4 sr", "");
  else if(R>0.25) model_info->AddEntry("", "A_{reco jet} > 0.2 sr", "");
  else model_info->AddEntry("", "A_{reco jet} > 0.09 sr", "");
      model_info->DrawClone("same");

  
  craa->SaveAs(str.Data());
  
  /*
  //COMPARE RAA
  TCanvas *craaR = new TCanvas("craaR","craaR",10,10,can_x,can_y);
  craaR->cd();

  for(int loop=0;loop<nloops;loop++){
  hraaR[loop]=(TH1D*)hunfoldedrb[loop]->Clone(Form("hraaR_%i",loop));
  hraaR[loop]->Divide(htruthprb);
  hraaR[loop]->SetLineColor(kBlue);
  hraaR[loop]->SetMarkerStyle(marker);
  hraaR[loop]->SetMarkerColor(kBlue);
  hraaR[loop]->SetMarkerSize(marker_size);
  hraaR[loop]->SetLineWidth(2);
  
}*/
}