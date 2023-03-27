#include "../utils/fit_functions.C"
#include "../utils/utils.C"
#include "../utils/statistic_tests.C"

bool skipiter(Int_t iter, bool svd=0)
{
  Int_t it=iter+1;
	if((it==2 || it==2|| it==3 || it==4 ||it==5 ||it==6) && !svd)
		 return kFALSE;
	else if((it==2 || it==3 || it==4|| it==5 || it==6 || it==6 || it==8) && svd)
		return kFALSE;
	else 
		return kTRUE;
}

void show_unfold(Int_t doToymodel=0, short peripheral=0, Float_t pTthresh=5.0, Float_t R=0.2, Int_t priorNo=2, Int_t bfIter=4, Int_t SVD1=0, TString insuf="eta"/*input dir suffix*/, TString typesuf="normal", bool uneqbin=1,  int biningch=1, TString ident="",Int_t doCompare=0, Int_t priorNo2=4, Int_t bfIter2=4, TString insuf2="eta", Int_t SVD2=0, TString trigger="MB")
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
  
  //Int_t SVD1=0; //first data set is SVD unfolding
  //bool peripheral=0; //central or peripheral collisons
  Int_t doToymodel2=0;
  Int_t measured2=0; //do we compare also different measured distributions?
  Int_t rebin2=2; //how many bins to merge 
  Int_t showBF=1; //show backfolded distributions
  Int_t showPP=1; //show pp distribution
  Int_t showHardJet=1; //show reconstructed hard jet distribution (toymodel only)
  Int_t compare3=0;  //RAA comparison: 0 nothing, 1 different efficiencies, 2 different R's, 3 different pTleadings 4:different priors 5:different iterations
  Int_t doRebin2=0; //Rebin second set of histograms
  if(SVD1==1 && doCompare==1 && SVD2==0 && !uneqbin) doRebin2=1;
  Double_t pTcutoff=100.0;
  Double_t RAA=0.5; //(RAA of TOYMODEL)
  if(peripheral)RAA=0.8;

  Color_t col_unf=kMagenta; //color for unfolded histograms
  
  Int_t unftype=1; //0: BG effect only, 1: BG+detector effects, 2: detector effects only
  Int_t unftype2=1;
  
//  TString typesuf="_normal";
//  TString typesuf="_u";
//  TString typesuf="_g";
  
  TString insuf_toy=insuf;
  TString frag_type="u";
  TString nbins1="200";
  TString nbins2="200";
  if(SVD1) nbins1="100";
  if(SVD2) nbins2="100";
  if(uneqbin){
	  nbins1="VAR";
	  nbins2="VAR";
  }

   //unfolding iterations
  Int_t firstiter=1;//first iteration (starting from 1!)
  Int_t lastiter=7; //last iteration (starting from 1!)
  //if(doToymodel)lastiter=10;
    //if(SVD1)lastiter=1;
  Int_t firstiter2=1;//first iteration for 2nd dataset
  Int_t lastiter2=7; //last iteration for 2nd dataset
    //if(SVD2)lastiter2=1;
  //Int_t bfIter=3; //which iteration to backfold and to show in comparison histograms (first=1, not 0)
  //if(SVD1)bfIter=4;
  //Int_t bfIter2=2;
  
  Int_t compIter=bfIter; //to which iteration compare all other iterations
  Int_t compIter2=bfIter2;
  Int_t compIter3=bfIter2;
  Int_t priorNo3=priorNo2;
  
  
  int statistics=10000; //number of hits in backfolding
  
  TString prior_type[]={"measured","flat","pythia","powlaw3","powlaw45","powlaw5","powlaw55","tsalis_1","tsalis_2","tsalis_3","tsalis_4","tsalis_5","tsalis_6","tsalis_7","tsalis_8","tsalis_9"};
  TString prior_type_name[]={"measured", "flat","biased Pythia","1/pT^{3}","1/pT^{4.5}","1/pT^{5}","1/(pT)^{5.5}", "tsalis_1","tsalis_2","tsalis_3","tsalis_4","tsalis_5","tsalis_6","tsalis_7","tsalis_8","tsalis_9"};
  if(!doToymodel) prior_type_name[0]="unfolded-Bayes";
    
  TString UnfType[]={"BG_sp","BGD","D"};
  TString desctype[]={"BG", "BG+detector effects","detector effects"};
  TString desc;
  if(doToymodel) desc=Form(", frag.: %s, corr.: %s",frag_type.Data(),desctype[unftype].Data());
else desc=Form(",corr.: %s",desctype[unftype].Data());
  TString desc2;
  if(priorNo!=priorNo2)
  desc2=Form("%s vs %s",prior_type_name[priorNo].Data(),prior_type_name[priorNo2].Data()); //title in comparison histograms
  else if(SVD1!=SVD2) desc2="SVD vs Bayesian unf.";
  if(compare3==1) desc2="#epsilon +-5\%.";
  else if(compare3==2) desc2="R";
  else if(compare3==3) desc2=Form("p_{T}^{leading}, R=%.1lf",R);
  else if(compare3==4) desc2=Form("prior comparison",R);
  else if(compare3==5) desc2=Form("different iterations",R);

  if(peripheral==1) trigger+="_peripheral";
	else if(peripheral==2)trigger+="_pp";
  TString ext="gif"; //output file extension
  //TString comment=UnfType[unftype]; //will be attached to the end of the output filename and also detrmines the name of the subdirectory with figures
  TString comment=ident; //will be attached to the end of the output filename 
  TString comment2=ident;
  if(priorNo!=priorNo2)
  comment2=Form("%svs%s_%s",prior_type[priorNo].Data(),prior_type[priorNo2].Data(),ident.Data()); //will be attached to the end of the output filename (in comparison histograms)
  TString leg_name1="distribution1:"; //will be in legend as a small title
  TString leg_name2="distribution2:";//will be in legend as a small title
  //TString comment2="quarterVSmid"; //will be attached to the end of the output filename (in comparison histograms)
  TString prefix="unf_charged"; //beginning of the output filenames
  TString unf_type="Bayes"; //will be in the title of histograms
  if(SVD1==1) unf_type="SVD";

  //LEVY FUNCTION
  TF1* LevyFit_pT = new TF1("LevyFit_pT",LevyFitFunc,10.0,40.0,9);
	LevyFit_pT->SetParameter(0,2.24); // B, changes the amplitude, 0.1
	LevyFit_pT->SetParameter(1,0.99); // T, changes slope, 0.4
	LevyFit_pT->SetParameter(2,18.1); // n, changes how fast spectrum drops, 5.8
	LevyFit_pT->SetParameter(3,0.000001); // m0, changes the width, 0.0001		
	LevyFit_pT->SetParameter(4,-4.7); // mu, changes the x-axis shift
	LevyFit_pT->SetParameter(5,30);//pT0
	LevyFit_pT->SetParameter(6,2);//R
	LevyFit_pT->SetParameter(7,1);//a
	LevyFit_pT->SetParameter(8,1);//b
	LevyFit_pT->SetRange(0,50);
	LevyFit_pT->SetLineWidth(2);
	LevyFit_pT->SetLineStyle(1);
	LevyFit_pT->SetLineColor(kGreen);
	LevyFit_pT->SetNpx(1000);
	
	TF1* fjet=new TF1("fjet","[0]*TMath::Power(x,-[1])",10,40);

   fjet->SetParameter(0,0.1);   
   fjet->SetParameter(1,5);   

	fjet->SetLineWidth(2);
	fjet->SetLineStyle(1);
	fjet->SetLineColor(kGreen);
	fjet->SetNpx(1000);
	//-----------------------------

  TString wrkdira;
  TString wrkdira2;
  TString wrkdira3;
  TString wrkdir;
  TString wrkdir2;
  TString wrkdir3;
  TString str;
  
  Double_t R2=R;
  Double_t pTthresh2=pTthresh;
  Double_t R3=R;
  Double_t pTthresh3=pTthresh;
  if(doCompare) {
    Double_t R3=R;
    Double_t pTthresh3=pTthresh;
    if(compare3==2)
    {
      R2=0.3;
      R3=0.2;
    }
    else if(compare3==3)
    {
      pTthresh2=5;
      pTthresh3=4;
    }
    else if(compare3==4)
    {
      priorNo3=8;
    }
    else if(compare3==5)
    {
      compIter2=compIter-1;
		compIter3=compIter+1;
    }
  }
  
  
  TString RAAsuf="";
  if(doToymodel)
  {
	 RAAsuf=Form("_RAA%.1lf",RAA);
    wrkdir="../../plotting_out/root/toymodel";
	 if(peripheral)wrkdir="../../plotting_out/root/toymodel_peripheral";
  }
    else{
     wrkdir=Form("../../plotting_out/root/%s",trigger.Data());
	}
  if(doToymodel2)
  {
    wrkdir2="../../plotting_out/root/toymodel";
	 if(peripheral)wrkdir2="../../plotting_out/root/toymodel_peripheral";
    wrkdir3="../../plotting_out/root/toymodel";
	 if(peripheral)wrkdir3="../../plotting_out/root/toymodel_peripheral";
  }
  else{
     wrkdir2 = Form("../../plotting_out/root/%s",trigger.Data());
     wrkdir3 = Form("../../plotting_out/root/%s",trigger.Data());
}
	TString unfname[]={"Bayes","SVD"};
  wrkdir = Form("%s/%s/Unfolded_R%.1lf_%s_%sbins_bining%i_%s%s_%s_%s",wrkdir.Data(),insuf.Data(),R,unfname[SVD1].Data(),nbins1.Data(),biningch,UnfType[unftype].Data(),RAAsuf.Data(),insuf.Data(),typesuf.Data());

  wrkdir2 = Form("%s/%s/Unfolded_R%.1lf_%s_%sbins_bining%i_%s%s_%s_%s",wrkdir2.Data(),insuf2.Data(),R2,unfname[SVD2].Data(),nbins2.Data(),biningch,UnfType[unftype2].Data(),RAAsuf.Data(),insuf2.Data(),typesuf.Data());
  wrkdir3 = Form("%s/%s/Unfolded_R%.1lf_%s_%sbins_bining%i_%s%s_%s_%s",wrkdir3.Data(),insuf2.Data(),R3,unfname[SVD2].Data(),nbins2.Data(),biningch,UnfType[unftype2].Data(),RAAsuf.Data(),insuf2.Data(),typesuf.Data());
    
   if(!doToymodel){
    
    if(compare3==1){
    wrkdir3=Form("%s_m5_abs",wrkdir3.Data());
    wrkdir2=Form("%s_p5_abs",wrkdir2.Data());
    }
    
    }

  
//cout<<"workingdir: "<<wrkdir<<endl;
//cout<<"workingdir2: "<<wrkdir2<<endl;
  
  //TString responsedir="../../plotting_out/root/response_matrix";
  TString responsedir=wrkdir;
  TString responsedir2=wrkdir2;
  
  TString outdir = Form("../../plotting_out/obr/%s/%s/R%.1lf/pTlead%0.lf/%s",trigger.Data(),UnfType[unftype].Data(),R,pTthresh,unfname[SVD1].Data());
  if(doToymodel) outdir = Form("../../plotting_out/obr/unfolding_toymodel/%s/R%.1lf/pTlead%0.lf/%s",UnfType[unftype].Data(),R,pTthresh,unfname[SVD1].Data());

  //TString wrkdir_data = Form("../../plotting_out/root/%s",trigger.Data());
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
  Float_t priorRange=40; //up to which value do we want to plot th(no efficiency reconstruction correction)e prior distribution?
  
  Float_t spectraXmin=-20;
  Float_t spectraXmax=50;
  Float_t spectraYmin=1E-10;
  Float_t spectraYmax=1E-2;
  
  bfIter=bfIter-1; //becasue arrays start from 0
  bfIter2=bfIter2-1;
  compIter=compIter-1;
  compIter2=compIter2-1;
  compIter3=compIter3-1;
  firstiter=firstiter-1; //(change start from 1 to 0)
  firstiter2=firstiter2-1; //(change start from 1 to 0)
  const Int_t niter = lastiter;//number of iterations
  const Int_t niter2 = lastiter2;//number of iterations for 2nd dataset
  //*****end of SETUP*****
  
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.04);
  Color_t colorList[]={kBlack,kRed,kGreen+3,kMagenta+2,kBlue+3,kOrange,kRed+3,kBlue-4,kOrange-1,kGreen,kYellow+2,7,kGreen+4,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,28};
  Color_t colorList2[30]={kRed,kBlue-4,kOrange,kGreen,13,14,kMagenta+2,kGreen+3,15,kBlue,kOrange+2,kYellow+2,kRed+2,kBlack-2,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};
  
  TH1 *frame = new TH1I("frame", "", 1000, -100, +100);
  TH1 *frame2 = new TH1I("frame", "", 1000, -100, +100);


  str = Form("%s/../histos_inclusivejet_%s.root", wrkdir.Data(),typesuf.Data());
  if(doToymodel) str = Form("%s/../root_RAA%.1lf_%s/histos_jets_R%.1lf_pTcut0.2.root", wrkdir_true.Data(),RAA,insuf_toy.Data(),R);
  TFile *f = new TFile(str.Data(), "OPEN");
  TH1I* hevents;
  if(!doToymodel)hevents= (TH1I*) f->Get("hevts");
  else hevents= (TH1I*) f->Get("hevents");

  Int_t nevents=hevents->GetEntries();
  //Int_t nevents=5E6;
  Int_t nevents2=hevents->GetEntries();
  //Int_t nevents2=5E6;
 
  Int_t toy_events_jetonly;
  if(doToymodel)
  {
	  str = Form("%s/../root_RAA%.1lf_%s/jetonly_R%.1lf_pTcut0.2.root", wrkdir_true.Data(),RAA,insuf_toy.Data(),R);
	  TFile *fj = new TFile(str.Data(), "OPEN");
	  TH1I* hevents=(TH1I*) fj->Get("hevents");
	  toy_events_jetonly=hevents->GetEntries();
  }
  
  Int_t toy_events=nevents; //number of events for toymodel-jetonly-distribution
  cout<<"M events: "<<((double) nevents)/1000000.0<<endl;
  //Int_t toy_events=10E6;
  //if(doToymodel) toy_events=2E7; //
  //Double_t hole=(1.0/12.)+((R*2.)/(2.*TMath::Pi())); //fraction of acceptance which was droped due to the bad sectors
  Double_t hole=0;
  Double_t hole2=hole;
  cout<<"hole: "<<hole<<endl;
  if(doToymodel)hole=0;
  if(doToymodel2)hole2=0;
  Double_t scale_jets = 1./(2*(1.-R)*2.*TMath::Pi()*nevents*(1-hole));
  //Double_t scale_jets = 1;
  Double_t scale_jets2 =1./(2*(1.-R)*2.*TMath::Pi()*nevents2*(1-hole2));
  //Double_t scale_jets =1;
  Double_t scale_true_toy=1./(2.*(1.-R)*2.*TMath::Pi()*toy_events);
  Double_t scale_true_toy_jetonly=1./(2.*(1.-R)*2.*TMath::Pi()*toy_events_jetonly);

  Double_t TAA=22.75;
  if(peripheral==1) TAA=0.49;
  if(peripheral==2) TAA=1;
  Double_t ppXsection=42;
  Double_t scale_true=TAA;//*ppXsection/*/(2*TMath::Pi())*/;
  Double_t scale_true2=(TAA+1.56);//*ppXsection/*/(2*TMath::Pi())*/;
  Double_t scale_true3=(TAA-1.56);//*ppXsection/*/(2*TMath::Pi())*/;
  //Double_t binWidth;


  // UNFOLDING RESULTS
  str = Form("%s/%s/unfolded_SignalSpectrum_R%.1lf_pTthresh%.1lf.root", wrkdir.Data(),prior_type[priorNo].Data(), R, pTthresh);
  //if(doToymodel) str = Form("%s/%s/unfolded_SignalSpectrum_pTthresh%.1lf.root", wrkdir.Data(),prior_type[priorNo].Data(), pTthresh);
  //if (SVD1) str = Form("%s/%s/unfolded_SVD_k%i_R%.1lf_pTthresh%.1lf_%i.root", wrkdir.Data(),prior_type[priorNo].Data(),kTerm, R, pTthresh, nSVDbins);
  //if (SVD1) str = Form("%s/%s/unfolded_SVD_R%.1lf_pTthresh%.1lf.root", wrkdir.Data(),prior_type[priorNo].Data(), R, pTthresh);
  TFile *funfolding = new TFile(str.Data(), "OPEN");
  cout<<"opening file: "<<str<<endl;
  
  //2nd and 3rd distributions for comparison
  TFile *funfolding2;
  TFile *funfolding3;
  if(doCompare){
    
    str = Form("%s/%s/unfolded_SignalSpectrum_R%.1lf_pTthresh%.1lf.root", wrkdir2.Data(),prior_type[priorNo2].Data(), R2, pTthresh2);
    funfolding2= new TFile(str.Data(), "OPEN");
    str = Form("%s/%s/unfolded_SignalSpectrum_R%.1lf_pTthresh%.1lf.root", wrkdir3.Data(),prior_type[priorNo3].Data(), R3, pTthresh3);
    funfolding3 = new TFile(str.Data(), "OPEN");
  }

  //RESPONSE MATRIX
  /*str=Form("%s/response_matrix_deltapT_R%.1lf.root", responsedir.Data(),R);
  if(doToymodel)str=Form("%s/response_matrix_deltapT_R%.1lf_TOY.root", responsedir.Data(),R);
  TFile *frmatrix = new TFile(str.Data(), "OPEN");
  str=Form("%s/response_matrix_deltapT_R%.1lf.root", responsedir2.Data(),R);
  TFile *frmatrix2 = new TFile(str.Data(), "OPEN");*/

  TH1D *hunfolded[niter];
  TH1D *hunfolded2[niter2];
  TH1D *hunfolded3[niter2];
  TH1D* hratut[niter];
  TH1D* hratut2[niter];

  TH2D *hcovariance[niter];
  TH2D *hcovariance2[niter2];
  TH2D *hcorrelation[niter];
  TH2D *hcorrelation2[niter2];
	
  TH1D* hdvec[niter];
  
  Double_t integral[niter]; //integral of measured distribution 
  Double_t integral10[niter]; //integral of measured distribution starting from 10GeV
  Double_t integral2[niter2]; //integral of unfolded2 distributions
  Double_t integralm; //integral of measured dist.
  Double_t integralm2; //integral of measured2 dist.
  Double_t integralt; //integral of true distribution
  Double_t int_bfold; //integral of backfolded distributions
  Double_t int_bfold2; //integral of backfolded2 distributions
  TDirectoryFile *dir;
  
  //efficiency histograms
	/*
  if (unftype == 2 || unftype == 1){
  str = Form("../../plotting_out/root/%s/epsilon/epsilon_R%.1lf_pTlead%.0lf.root",trigger.Data(),R,pTthresh);
  f = new TFile(str.Data(), "OPEN");
  TH1D * hepsilon=(TH1D*)f->Get("hepsilon_unfolded");
	}
  */
  

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
//cout<<"scaler: "<<scale_jets<<endl;
  
  /*Double_t wdth1=hmeasured->GetBinWidth(hmeasured->FindBin(0));
  Double_t wdth2=hmeasured->GetBinWidth(hmeasured->FindBin(80));
  
  if(wdth2>(wdth1+0.0001)) //input histos have unequal binning
  {
    uneqbin=true;
    if(doToymodel)binWidth=0.1;  
    else binWidth=0.5;
  }
  else
    binWidth=1;
  scale_jets=scale_jets/binWidth;
  scale_jets2=scale_jets2/binWidth;*/
  hmeasured->Scale(scale_jets,"width");
  integralm=hmeasured->Integral("width");
  Double_t int10_50=hmeasured->Integral(hmeasured->FindBin(10),hmeasured->FindBin(50),"width");
  cout<<"integral measured:"<<integralm<<endl;
  cout<<"integral measured (10-50GeV): "<<int10_50<<endl;
  
  //*****2nd measured dist.*****
  if(measured2 && doCompare)
  {
  TH1D *hmeasured2;
  hmeasured2= (TH1D*)((TDirectoryFile*)funfolding2->Get("input"))->Get("hmeasured");
  //else hmeasured2 =(TH1D*) funfolding2->Get("hmeasured");
  hmeasured2->Sumw2();
  hmeasured2->SetMarkerStyle(markerM2);
  hmeasured2->SetLineColor(colorList2[0]);
  hmeasured2->SetLineWidth(line_width);
  hmeasured2->SetMarkerColor(colorList2[0]);
  hmeasured2->SetMarkerSize(marker_size);
  hmeasured2->Scale(scale_jets2,"width");
  integralm2=hmeasured2->Integral("width");
  }
  
  //******unfolded spectra******
  for(Int_t iter = firstiter; iter < lastiter; iter++)
    {
      //if(skipiter(iter,SVD1) && SVD1) continue;
      cout<<"iter"<<iter<<endl;
      /*if(!SVD1)*/	str = Form("iter%d", iter);
      //else str = Form("kterm%d", iter+1);
      dir = (TDirectoryFile*)funfolding->Get(str.Data());
      hunfolded[iter] = (TH1D*)dir->Get("hunfolded");
		hunfolded[iter]->Sumw2();
      //hunfolded[iter]->SetLineColor(colorList[iter+3]);
      //hunfolded[iter]->SetMarkerColor(colorList[iter+3]);
      hunfolded[iter]->SetLineColor(colorList[iter+3]-firstiter);
      hunfolded[iter]->SetMarkerColor(colorList[iter+3]-firstiter);
		if(iter==bfIter)
		{
			hunfolded[iter]->SetLineColor(col_unf);
			hunfolded[iter]->SetMarkerColor(col_unf);
		}
      hunfolded[iter]->SetMarkerStyle(markerU);
      hunfolded[iter]->SetMarkerSize(marker_size);
      hunfolded[iter]->SetLineWidth(line_width);
     /* if(!SVD1 && !Roo1){
	str = Form("HCovariance_%s", errtype.Data());
	hcovariance[iter] = (TH2D*)dir->Get(str.Data());

	str = Form("HCorrelation_%s", errtype.Data());
	hcorrelation[iter] = (TH2D*)dir->Get(str.Data());

	str = Form("N_{iter} = %d", iter+1);
	hunfolded[iter]->SetTitle(str.Data());}*/
     // else if(SVD1 && !Roo1){
	

	hcovariance[iter] = (TH2D*)dir->Get("hcovariance");
	if(SVD1)
	hdvec[iter]=(TH1D*) dir->Get("hdvec");
      //}
    
   /*
      TH1D *htemp = (TH1D*)hunfolded[iter]->Clone("htemp");
      hunfolded[iter]->Reset("MICE");
		hunfolded[iter]->Sumw2();
      for(Int_t bin = 1; bin <= htemp->GetNbinsX(); bin++)
	{
	  Double_t yield = htemp->GetBinContent(bin);
	  if((htemp->GetBinCenter(bin))>pTcutoff)continue;
	  Double_t error;
	  error = TMath::Sqrt(hcovariance[iter]->GetBinContent(bin, bin));
	  hunfolded[iter]->SetBinContent(bin, yield);
	  hunfolded[iter]->SetBinError(bin, error);
	}
      delete htemp;
	
*/
      //if(SVD1)hunfolded[iter]->Scale(integralm/hunfolded[iter]->Integral("width"));
      
     /* TH1D* herror=hunfolded[iter]->Clone("herror");
		herror->Reset("MICE");
		for(int bin=1; bin<herror->GetNbinsX();bin++)
		{
			Double_t relerr=0;
			if(hunfolded[iter]->GetBinContent(bin)>0)relerr=hunfolded[iter]->GetBinError(bin)/hunfolded[iter]->GetBinContent(bin);
			herror->SetBinContent(bin,relerr);
		}*/
      hunfolded[iter]->Scale(scale_jets,"width");
		/*for(int bin=1; bin<herror->GetNbinsX();bin++)
		{
			hunfolded[iter]->SetBinError(bin,herror->GetBinContent(bin)*hunfolded[iter]->GetBinContent(bin));
		}*/
      integral[iter]=hunfolded[iter]->Integral("width");
      integral10[iter]=hunfolded[iter]->Integral(hunfolded[iter]->FindBin(10),hunfolded[iter]->GetNbinsX(),"width");

      cout<<iter<<" integral unfolded:"<<integral[iter]<<endl;
		//cout<<iter<<" integral unfolded/integral pp (from 10GeV):"<<integral10[iter]/integralt10<<endl;
		
		/*
		
		Int_t fbin=hunfolded[iter]->GetXaxis()->FindBin(pTthresh+0.001);
		Double_t yield = hunfolded[iter]->GetBinContent(fbin);
      hunfolded[iter]->SetBinContent(fbin, yield*1.0);
		//hunfolded[iter]->Scale(integral[iter]/hunfolded[iter]->Integral("width"));
		
		fbin=fbin+1;
		yield = hunfolded[iter]->GetBinContent(fbin);
      hunfolded[iter]->SetBinContent(fbin, yield*1.0);
		//hunfolded[iter]->Scale(integral[iter]/hunfolded[iter]->Integral("width"));
		
		fbin=fbin+1;
		yield = hunfolded[iter]->GetBinContent(fbin);
      hunfolded[iter]->SetBinContent(fbin, yield*1.0);
		hunfolded[iter]->Scale(integral[iter]/hunfolded[iter]->Integral("width"));
		
		integral[iter]=hunfolded[iter]->Integral("width");
      cout<<"integral unfolded:"<<integral[iter]<<endl;
		*/
    }
    
       //******toymodel true distribution******
    if(doToymodel){
		 //PYTHIA charged jet spectrum with pTlead>5GeV 
		Double_t jet_norm=TAA;
		Double_t B[3] = {2.41349e+01,1.75255e+01,1.27430e+01}; //R=0.2,0.3,0.4
		Double_t T[3] = {2.69190e+00,2.53251e+00,2.40521e+00};
		Double_t n[3] = {1.42495e+02,1.03738e+02,8.14043e+01};
		Double_t m0[3] = {-3.00000e+00,-7.83741e+00,-8.40985e+00};
		Double_t mu[3] = {-2.97634e+01,-2.01221e+01,-1.62948e+01};
		Double_t A[3] = {1.12732e+00,5.45759e-01,3.25500e-01};
		Double_t pwr[3] = {4.82216e+00,4.37727e+00,4.08808e+00};

		Int_t ridx=R*10-2; //0:R=0.2 1:R=0.3 2:R=0.4

		TF1 *fjettrue=new TF1("fjettrue",hardjet,3,40,9);
		fjettrue->SetParameters(jet_norm,B[ridx],T[ridx],n[ridx],m0[ridx],mu[ridx],A[ridx],pwr[ridx],RAA);
		fjettrue->SetNpx(10000);
		fjettrue->SetLineColor(kMagenta);
		/*
		Double_t A_tsal=3.0*TAA*RAA;
		Double_t n_tsal=10.7;
		Double_t T_tsal=0.44;
  
		TF1 *fjettrue=new TF1("fjettrue",TsalisFitFunc,3,40,3);
		fjettrue->SetParameters(A_tsal,n_tsal,T_tsal);
		fjettrue->SetNpx(10000);
		fjettrue->SetLineColor(kMagenta);*/
		
	 }
	
	//True jet distribution (TOYMODEL only)
	TH1D* htruth;
	TH1D *hjettruth;
	if(doToymodel && showHardJet){
		//particle level
		str = Form("%s/../root_RAA%.1lf%s/control_histos.root", wrkdir_true.Data(),RAA,insuf_toy.Data());
		TFile *fcontrol = new TFile(str.Data(), "OPEN");
		Double_t scale_jets_toy_gen = 1./(2*2.*TMath::Pi()*toy_events);
		TH2D*  hgen2d=(TH2D*) fcontrol->Get("hpTpTleadGen");
		TH1D* hgen=(TH1D*)hgen2d->ProjectionY("hgen",hgen2d->GetXaxis()->FindBin(pTthresh),hgen2d->GetXaxis()->GetNbins());
		hgen->Scale(scale_jets_toy_gen,"width");
		htruth=(TH1D*)rebinhisto(hgen, hunfolded[0], "htruth_gen");
		htruth->SetMarkerColor(kRed+1);
		htruth->SetLineColor(kRed+1);
		htruth->SetMarkerStyle(kFullSquare);
		htruth->SetMarkerSize(marker_size);
		htruth->SetLineWidth(line_width);
	
	//reconstructed jet level
  str = Form("%s/../root_RAA%.1lf%s/jetonly_R%.1lf_pTcut0.2.root", wrkdir_true.Data(),RAA,insuf_toy.Data(),R);
  TFile *ftruth = new TFile(str.Data(), "OPEN");
  TH2D *hPtRecpTleadingPrior = (TH2D*)ftruth->Get("fhPtRecpTleading");
  Int_t firstbin = hPtRecpTleadingPrior->GetXaxis()->FindBin(pTthresh);
  Int_t lastbin = hPtRecpTleadingPrior->GetNbinsX();
  TH1D *htemp = hPtRecpTleadingPrior->ProjectionY("htemp", firstbin, lastbin);
  //float test_cont=htemp->GetBinContent(htemp->FindBin(10));
  //cout<<"test content:"<<test_cont<<endl;
  htemp->Scale(scale_true_toy_jetonly,"width");
  hjettruth = (TH1D*) rebinhisto(htemp,hunfolded[0],"hjettruth");
  delete htemp;
  hjettruth->SetMarkerColor(kRed-1);
  hjettruth->SetLineColor(kRed-1);
  hjettruth->SetMarkerStyle(kFullSquare);
  hjettruth->SetMarkerSize(marker_size);
  hjettruth->SetLineWidth(line_width);
  //hjettruth->Scale(integralm/hjettruth->Integral("width"));
  integralt=hjettruth->Integral("width");
}//toymodel hard jet distribution
 
	

	//******Prior distribution******
  TH1D *hprtmp;
  hprtmp = (TH1D*)((TDirectoryFile*)funfolding->Get("input"))->Get("hprior");
  //else hprtmp=(TH1D*) funfolding->Get("hprior");
  TH1D *hprior = rebinhisto(hprtmp, hunfolded[0], "hprior");
  /*TH1D *hprior = (TH1D*)hmeasured->Clone("hprior");
  hprior->Sumw2();
  hprior->Reset("MICE");
  for(Int_t bin = 1; bin <= hprtmp->GetNbinsX(); bin++)
    {
      Double_t pTtmp = hprtmp->GetBinCenter(bin);
      Double_t yield = hprtmp->GetBinContent(bin);
      hprior->Fill(pTtmp, yield);
    }*/
  delete hprtmp;
  /*for(Int_t bin = 1; bin <= hprior->GetNbinsX(); bin++)
    {
      Double_t yield = hprior->GetBinContent(bin);
      hprior->SetBinError(bin,TMath::Sqrt(yield));
    }*/
  hprior->SetLineColor(colorList[1]);
  hprior->SetLineWidth(line_width);
  //hprior->SetMarkerStyle(markerP);
  //hprior->SetMarkerColor(colorList[1]);
  //hprior->SetMarkerSize(marker_size);
  hprior->Scale(scale_jets,"width");
  Double_t int_prior=hprior->Integral("width"); 
  if(int_prior>0)hprior->Scale(integralm/int_prior);
  cout<<"prior integral:"<<int_prior<<endl;
    
  //******2nd prior******
  if(doCompare){
  TH1D *hprtmp;
  hprtmp = (TH1D*)((TDirectoryFile*)funfolding2->Get("input"))->Get("hprior");
  //else hprtmp=(TH1D*) funfolding2->Get("hprior");
  TH1D *hprior2 = rebinhisto(hprtmp, hunfolded[0], "hprior2");
  delete hprtmp;

  hprior2->SetLineColor(colorList2[1]);
  hprior2->SetLineWidth(line_width);
  //hprior2->SetMarkerStyle(markerP2);
  //hprior2->SetMarkerColor(colorList2[1]);
  //hprior2->SetMarkerSize(marker_size);
  hprior2->Scale(scale_jets2,"width");
  Double_t int_prior2=hprior2->Integral("width"); 
  if(int_prior2>0)hprior2->Scale(integralm/int_prior2);
  }//doCompare
    
      //pp baseline
  if(showPP){
	  /*
  str = "../../plotting_out/root/pp_matt/PythiaHistosChargedJets2.root";
  ftruth = new TFile(str.Data(), "OPEN");*/
  str="../../plotting_out/root/pp_Fuqiang/Jan_fastjet3/starjet_pythia_xsection_charged_jets.root";
  //str = Form("../../plotting_out/root/pp_pythia8/histos_pythiajet_R%.1lf.root",R);
  //cout<<"loading file: "<<str.Data()<<endl;
  ftruth = new TFile(str.Data(), "OPEN");
  str=Form("hpT_R0%.0lf_pTl%.0lf",R*10,0/*pTthresh*/);
  //cout<<"loading histogram: "<<str.Data()<<endl;
  TH1D *htemp =(TH1D*) ftruth->Get(str.Data());


  str="hevts";
  TH1I* hpevents=(TH1I*) ftruth->Get(str.Data());
  Int_t npevents=1; //Pythia spectrum is already normalized; //hpevents->GetEntries()/14.0;

  TH1D *hppbase = (TH1D*)hunfolded[0]->Clone("hppbase");
  hppbase->Reset("MICE");
  hppbase->Sumw2();
  for(Int_t bin = 1; bin < htemp->GetNbinsX()+1; bin++)
    {
      Double_t pTtmp = htemp->GetBinCenter(bin);
      Double_t yield = htemp->GetBinContent(bin);
      Double_t oldWidth=htemp->GetBinWidth(bin);
      Int_t newBin=hppbase->FindBin(pTtmp);
      Double_t newWidth=hppbase->GetBinWidth(newBin);
      //cout<<"bin:"<<binx<<" pT:"<<pT<<" yield:"<<yield<<" width: "<<binWidth<<" -> "<<newWidth<<" new yield:";
      yield=yield*(oldWidth/newWidth);
      hppbase->Fill(pTtmp, yield);
    }
  for(Int_t bin = 1; bin <= hppbase->GetNbinsX(); bin++)
    {
      Double_t yield = hppbase->GetBinContent(bin);
      Int_t bint=htemp->FindBin(hppbase->GetBinCenter(bin));
      Double_t error= (htemp->GetBinError(bint))/(htemp->GetBinContent(bint));
      hppbase->SetBinError(bin,error*yield);
    }
      delete htemp;


  hppbase->SetMarkerColor(kGray+1);
  hppbase->SetLineColor(kGray+1);
  hppbase->SetMarkerStyle(kFullSquare);
  hppbase->SetMarkerSize(marker_size);
  hppbase->SetLineWidth(line_width);
  hppbase->Sumw2();
  TH1D* hppbase2=(TH1D*) hppbase->Clone("hppbase2");
  TH1D* hppbase3=(TH1D*) hppbase->Clone("hppbase3");
  hppbase->Scale(scale_true/npevents);
  hppbase2->Scale(scale_true2/npevents);
  hppbase3->Scale(scale_true3/npevents);
  integralt=hppbase->Integral("width");
  Double_t integralt10=hppbase->Integral(hppbase->FindBin(10),hppbase->GetNbinsX(),"width");

  }
      //cout<<"DEBUG 4"<<endl;

  //******unfolded spectra for 2nd dataset******
  if(doCompare){
  for(Int_t iter = firstiter2; iter < lastiter2; iter++)
    {
      //if(skipiter(iter,SVD2) && SVD2) continue;
      str = Form("iter%d", iter);
        dir = (TDirectoryFile*)funfolding2->Get(str.Data());
	hunfolded2[iter] = (TH1D*)dir->Get("hunfolded");
       //hunfolded2[iter]->SetLineColor(colorList2[iter+3]);
      //hunfolded2[iter]->SetMarkerColor(colorList2[iter+3]);
      hunfolded2[iter]->SetLineColor(colorList2[iter+3]-firstiter2);
      hunfolded2[iter]->SetMarkerColor(colorList2[iter+3]-firstiter2);
      hunfolded2[iter]->SetMarkerStyle(markerU2);
      hunfolded2[iter]->SetMarkerSize(marker_size);
      hunfolded2[iter]->SetLineWidth(line_width);
      

     /* if(!SVD2 && !Roo1){
	str = Form("HCovariance_%s", errtype.Data());
	hcovariance2[iter] = (TH2D*)dir->Get(str.Data());

	str = Form("HCorrelation_%s", errtype.Data());
	hcorrelation2[iter] = (TH2D*)dir->Get(str.Data());

	str = Form("N_{iter} = %d", iter+1);
	hunfolded2[iter]->SetTitle(str.Data());}*/
      //else if(!Roo1){
	str = "hcovariance";
	  hcovariance2[iter] = (TH2D*)dir->Get(str.Data());
     // }
      
     // if(!Roo1){
      TH1D *htemp = (TH1D*)hunfolded2[iter]->Clone("htemp");
      hunfolded2[iter]->Reset("MICE");

      for(Int_t bin = 1; bin <= htemp->GetNbinsX(); bin++)
	{
	  Double_t yield = htemp->GetBinContent(bin);
	  if((htemp->GetBinCenter(bin))>pTcutoff)continue;
	  if(!doToymodel)Double_t error = TMath::Sqrt(hcovariance2[iter]->GetBinContent(bin, bin));
	  else error = TMath::Sqrt(yield);
	  hunfolded2[iter]->SetBinContent(bin, yield);
	  hunfolded2[iter]->SetBinError(bin, error);
	}
      delete htemp;
     //}
      hunfolded2[iter]->Scale(scale_jets2,"width");
      //if(SVD2)hunfolded2[iter]->Scale(integralm/hunfolded2[iter]->Integral("width"));
      //integral2[iter]=hunfolded2[iter]->Integral("width");  
    }
    
    //3rd unfolded
    for(Int_t iter = firstiter2; iter < lastiter2; iter++)
    {
      //if(skipiter(iter,SVD2) && SVD2) continue;
   	str = Form("iter%d", iter);
      dir = (TDirectoryFile*)funfolding3->Get(str.Data());
	hunfolded3[iter] = (TH1D*)dir->Get("hunfolded");
       //hunfolded3[iter]->SetLineColor(colorList2[iter+3]);
      //hunfolded3[iter]->SetMarkerColor(colorList2[iter+3]);
      hunfolded3[iter]->SetLineColor(colorList2[iter+4]-firstiter2);
      hunfolded3[iter]->SetMarkerColor(colorList2[iter+4]-firstiter2);
      hunfolded3[iter]->SetMarkerStyle(markerU2);
      hunfolded3[iter]->SetMarkerSize(marker_size);
      hunfolded3[iter]->SetLineWidth(line_width);


      TH1D *htemp = (TH1D*)hunfolded3[iter]->Clone("htemp");
      hunfolded3[iter]->Reset("MICE");

      for(Int_t bin = 1; bin <= htemp->GetNbinsX(); bin++)
	{
	  Double_t yield = htemp->GetBinContent(bin);
	  if((htemp->GetBinCenter(bin))>pTcutoff)continue;
	  if(!doToymodel) error = TMath::Sqrt(hcovariance2[iter]->GetBinContent(bin, bin));
	  else error = TMath::Sqrt(yield);
	  hunfolded3[iter]->SetBinContent(bin, yield);
	  hunfolded3[iter]->SetBinError(bin, error);
	}
      delete htemp;

      hunfolded3[iter]->Scale(scale_jets2,"width");
      //integral3[iter]=hunfolded3[iter]->Integral("width");
      
    }
      //cout<<"DEBUG 5"<<endl;

    
  }//doCompare
  
  //*******Response Matrix*************
  
   //TH2D * hresponse=(TH2D*) frmatrix->Get("hResponse_1E9");
  TH2D * hresponse;
  str = "input";
  dir = (TDirectoryFile*)funfolding->Get(str.Data());
  hresponse= (TH2D*)dir->Get("hresponse");

  if(doCompare){
  TH2D * hresponse2;
  str = "input";
  dir = (TDirectoryFile*)funfolding2->Get(str.Data());
  hresponse2= (TH2D*)dir->Get("hresponse");
  
  
  }
  //******backfolded distribution******
     //cout<<"DEBUG 6"<<endl;

  TH1D * hbfold[niter];
if(showBF){
  for(Int_t iter = firstiter; iter < lastiter; iter++)
    {
		str = Form("iter%d", iter);
      dir = (TDirectoryFile*)funfolding->Get(str.Data());
      hbfold[iter] = (TH1D*)dir->Get("hbackfolded");
		hbfold[iter]->Sumw2();
		hbfold[iter]->Scale(scale_jets,"width");
		
		hbfold[iter]->SetLineColor(colorList[iter+4]);
		hbfold[iter]->SetMarkerColor(colorList[iter+4]);
       if(iter==bfIter){
			hbfold[iter]->SetLineColor(colorList[2]);
			hbfold[iter]->SetLineWidth(line_width);
			hbfold[iter]->SetMarkerColor(colorList[2]);
			}
		hbfold[iter]->SetMarkerStyle(markerB);
		hbfold[iter]->SetMarkerSize(marker_size);
  
		   /*
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
		Double_t width=hunfolded[iter]->GetBinWidth(bin);
      Double_t yield = hunfolded[iter]->GetBinContent(bin);      
      if(yield==0)continue;
      yield=yield*width/statistics;
      Double_t pT=hunfolded[iter]->GetBinCenter(bin);
      if(pT>pTcutoff)continue;
      Int_t xbin=hresponse->GetYaxis()->FindBin(pT);

			
		if(unftype==2 || unftype==1){
			Double_t epsilon=hepsilon->GetBinContent(hepsilon->FindBin(pT));
		}
		TH1D *hsmear = hresponse->ProjectionX("hsmear", xbin, xbin);
      Double_t norm=hsmear->Integral("width");
      hsmear->Scale(1/norm);
      for(int i=0;i<statistics;i++){
			Double_t pTsmear = hsmear->GetRandom();
			if(unftype==2 || unftype==1){
				Double_t rand=gRandom->Uniform(0,1);
				if (rand>epsilon) continue;
			}
			Double_t wd=hbfold[iter]->GetBinWidth(hbfold[iter]->FindBin(pTsmear));
			hbfold[iter]->Fill(pTsmear, yield/wd);
      }
      delete hsmear;
      if(bin%10==0)cout<<"."<<endl;
    }//loop over bins  
    cout<<"done"<<endl;
  int_bfold=hbfold[iter]->Integral("width"); */
  }//loop over iterations
}//showBF
  
  double intbn=0;
	if(showBF)intbn=hbfold[bfIter]->Integral(1,hbfold[bfIter]->FindBin(pTthresh));
  double intmn=hmeasured->Integral(1,hmeasured->FindBin(pTthresh));
  
  double ratio=(intmn-intbn)/hmeasured->Integral();
  cout<<"ratio:"<<ratio<<endl;
    //cout<<"DEBUG 7"<<endl;

  
  //******backfolded distribution for 2nd prior******
  if(doCompare){
  //TH2D * hresponse2=(TH2D*) frmatrix2->Get("hResponse_1E9");
  TH1D * hbfold2[niter2];
  for(Int_t iter = firstiter2; iter < lastiter2; iter++)
    {
		 str = Form("iter%d", iter);
      dir = (TDirectoryFile*)funfolding2->Get(str.Data());
      hbfold2[iter] = (TH1D*)dir->Get("hbackfolded");
		hbfold2[iter]->Sumw2();
		hbfold2[iter]->Scale(scale_jets2,"width");
		
		hbfold2[iter]->SetLineColor(colorList2[iter+4]);
		hbfold2[iter]->SetMarkerColor(colorList2[iter+4]);
       if(iter==bfIter2){
			hbfold2[iter]->SetLineColor(colorList2[2]);
			hbfold2[iter]->SetLineWidth(line_width);
			hbfold2[iter]->SetMarkerColor(colorList2[2]);
			}
		hbfold2[iter]->SetMarkerStyle(markerB2);
		hbfold2[iter]->SetMarkerSize(marker_size);
		 /*
   if(iter!=bfIter2) continue;
  str = Form("hbfold2_%d", iter);
  hbfold2[iter]=(TH1D*)hunfolded2[iter]->Clone(str.Data());
  //hbfold2[iter]->SetLineColor(colorList2[iter+6]);
  //hbfold2[iter]->SetMarkerColor(colorList2[iter+6]);
  //if(iter==bfIter2){
    hbfold2[iter]->SetLineColor(colorList2[2]);
    hbfold2[iter]->SetLineWidth(line_width);
    hbfold2[iter]->SetMarkerColor(colorList2[2]);
  //}
  hbfold2[iter]->SetMarkerStyle(markerB2);
  hbfold2[iter]->SetMarkerSize(marker_size);
  hbfold2[iter]->Sumw2();
  hbfold2[iter]->Reset("MICE");
  cout<<"creating backfolded spectrum II:"<<endl;
  for(Int_t bin = 1; bin <= hunfolded2[iter]->GetNbinsX(); bin++)
    {
      Double_t yield = hunfolded2[iter]->GetBinContent(bin);      
      if(yield==0)continue;
      yield=yield/statistics;
      Double_t pT=hunfolded2[iter]->GetBinCenter(bin);
      if(pT>pTcutoff)continue;
      Int_t xbin;
      xbin=hresponse2->GetYaxis()->FindBin(pT);
       Double_t rnd=gRandom->Uniform(0,1);
	   if(unftype==2 || unftype==1){
			Double_t epsilon=hepsilon->GetBinContent(hepsilon->FindBin(pT));
		}
     // if(rnd>0.5+binCenter) xbin=xbin-1;
      TH1D *hsmear;
      hsmear = hresponse2->ProjectionX("hsmear", xbin, xbin);
      
      Double_t norm=hsmear->Integral("width");
      //cout<<"norm "<<norm<<endl;
      hsmear->Scale(1/norm,"width");
      
      for(int i=0;i<statistics;i++){
      Double_t pTsmear = hsmear->GetRandom();
	  if(unftype==2 || unftype==1){
				Double_t rand=gRandom->Uniform(0,1);
				if (rand>epsilon) continue;
			}
      hbfold2[iter]->Fill(pTsmear, yield);
      }
      delete hsmear;
      if(bin%10==0)cout<<"."<<endl;
    }//loop over bins  
    cout<<"done"<<endl;
    int_bfold2=hbfold2[bfIter2]->Integral("width"); */
  }//loop over iterations
  }//doCompare
    //cout<<"DEBUG 8"<<endl;

  
  //*****************************************
  // DRAWING SPECTRA
  //*****************************************

 if(pTcutoff>99)
  TString suffix=Form("R%.1lf_%sPrior_pTthresh%.1lf_bfIter%i",R, prior_type[priorNo].Data(), pTthresh, bfIter+1);
 //TString suffix=Form("R%.1lf_%sPrior_pTthresh%.1lf_bfIter5",R, prior_type[priorNo].Data(), pTthresh);
 else
  TString suffix=Form("R%.1lf_%sPrior_pTthresh%.1lf_bfIter%i_pTcutoff%.1lf",R, prior_type[priorNo].Data(), pTthresh, bfIter+1, pTcutoff);


  if(!doCompare){
   
  suffix=Form("%s_%s",suffix.Data(),comment.Data()); 

  //MEASURED, PRIOR, UNFOLDED, BACKFOLDED DISTRIBUTIONS
  TCanvas *cspectra = new TCanvas("cspectra","cspectra",10,10,can_x,can_y);
  cspectra->cd();
  cspectra->SetGrid();
  cspectra->SetLogy();
  frame->GetXaxis()->SetTitle("p_{T, corr}^{charged} (GeV/c)");
  frame->GetYaxis()->SetTitle("1/N_{events} 1/2#pi d^{2}N/dp_{T}^{ch}d#eta (GeV/c)^{-1}");
  frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
  frame->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
  //frame->GetYaxis()->SetRangeUser(1E-3, 1E5);
  str=Form("Unf.: %s, p_{T}^{leading} > %.1lf GeV/c, prior: %s%s",unf_type.Data(), pTthresh, prior_type_name[priorNo].Data(),desc.Data());
  frame->SetTitle(str);
  frame->DrawCopy("");
  hmeasured->DrawCopy("E same");
  //hunfolded[4]->Fit("LevyFit_pT","R");
  //hunfolded[4]->Fit("fjet","R");
  /*
  hppbase->Fit("fjet","R");
  //cout<<"levy integral"<<LevyFit_pT->Integral(10,50)<<endl;
  fjet->DrawClone("same");*/
  hprior->GetXaxis()->SetRangeUser(pTthresh,priorRange);
  hprior->DrawCopy("HIST SAME C");
  if(doToymodel && showHardJet)htruth->DrawCopy("E same");
  if(doToymodel && showHardJet)hjettruth->DrawCopy("E same");
  //if(doToymodel)fjettrue->DrawClone("same");
  if(showPP)hppbase->DrawCopy("E same");
  if(showBF) hbfold[bfIter]->DrawCopy("E same");
  for(Int_t iter = firstiter; iter < lastiter; iter++){
    if(skipiter(iter,SVD1))continue;
	 if(iter!=bfIter) continue;
    hunfolded[iter]->DrawCopy("E same");
    }

  //Info table
  TLegend *model_info = new TLegend(0.12037, 0.52, 0.3348, 0.900);
  model_info->SetFillStyle(0);
  model_info->SetBorderSize(0);
  model_info->SetMargin(0.05);
  model_info->SetHeader(Form("Run11 AuAu 200 GeV, %s",trigger.Data()));
  if(doToymodel) model_info->SetHeader("TOYMODEL: RHIC kinematics");
  //model_info->AddEntry("", "Single Particle", "");
  if (unftype<2){
		if(!peripheral) model_info->AddEntry("", "0-10% Central Collisions", "");
		else if(peripheral==1) model_info->AddEntry("", "60-80% Central Collisions", "");
		else if(peripheral==2) model_info->AddEntry("", "p+p collisions", "");
  }
  else if (unftype==2) model_info->AddEntry("", "hard jets only", "");
  if(doToymodel)model_info->AddEntry("", Form("N_{events} = %.1lfM", nevents/1E6), "");
  else model_info->AddEntry("", "Integrated luminosity: 6#mub^{-1}","");
  //model_info->AddEntry("", Form("N_{events} = %.1lfM", nevents/1E6), "");
  model_info->AddEntry("", Form("Anti-k_{T} \t \t R = %.1lf",R), "");
  model_info->AddEntry("", "p_{T}^{const} > 0.2 GeV/c", "");
   
  if(R>0.35) model_info->AddEntry("", "A_{reco jet} > 0.4 sr", "");
  else if(R>0.25) model_info->AddEntry("", "A_{reco jet} > 0.2 sr", "");
  else model_info->AddEntry("", "A_{reco jet} > 0.09 sr", "");
  //model_info->AddEntry("", "binning: 0.5GeV", "");
  //model_info->AddEntry("", "removed: bad sector + 1*R edge", "");
  model_info->DrawClone("same");
  
   //Legend
  TLegend *legspectra = new TLegend(0.6467, 0.60, 0.89, 0.90);
  legspectra->SetTextSize(0.03);
  legspectra->SetFillStyle(0);
  legspectra->SetBorderSize(0);
  legspectra->AddEntry(hmeasured, Form("Measured, I=%.3lf*%s",integralm*iscale, iscale_str.Data()), "lp");
  
  //if(doToymodel) legspectra->AddEntry(htruth,"T_{AA}*d#sigma_{pp}/dp_{T}", "lp");
  //if(doToymodel) legspectra->AddEntry(fjettrue,"generating function", "lp");
  if(doToymodel && showHardJet) legspectra->AddEntry(htruth,"hard spectrum (PL)", "lp");
  if(doToymodel && showHardJet) legspectra->AddEntry(hjettruth,"hard jets only", "lp");
  if(showPP) legspectra->AddEntry(hppbase,Form("T_{AA}*p+p (Pythia), I=%.3lf*%s",integralt*iscale, iscale_str.Data()), "lp");
  //legspectra->AddEntry(hprior,Form("prior dist., I=%.3lf*%s", int_prior*iscale, iscale_str.Data()), "l");
  legspectra->AddEntry(hprior,Form("prior: %s",prior_type_name[priorNo].Data()), "l");
  //legspectra->AddEntry(hprior,"truth", "l");
  //legspectra->AddEntry(hbfold[bfIter],Form("backfolded dist. (iter%i), I=%.3lf*%s", bfIter+1,int_bfold[bfIter]*iscale, iscale_str.Data()), "lp");
  if(showBF) legspectra->AddEntry(hbfold[bfIter],Form("backfolded dist. (iter%i)", bfIter+1), "lp");
  for(Int_t iter = firstiter; iter < lastiter; iter++){
    if(skipiter(iter,SVD1))continue;
	if(iter!=bfIter)continue;
    if(!SVD1)legspectra->AddEntry(hunfolded[iter], Form("Unfolded, N_{iter} = %i, I=%.3lf*%s ", iter+1, integral[iter]*iscale,iscale_str.Data()), "lp");
    else legspectra->AddEntry(hunfolded[iter], Form("Unfolded, kterm = %i, I=%.3lf*%s ", iter+1, integral[iter]*iscale,iscale_str.Data()), "lp");
    //legspectra->AddEntry(hunfolded[iter], Form("Unfolded, N_{iter} = %d", iter+1), "lp");
  }
  
  legspectra->DrawClone("same");
  
  latex->DrawLatex(0.35, 0.3,"STAR Preliminary");
  str = Form("%s/%s_spectra_%s.%s", outdir.Data(),prefix.Data(),suffix.Data(),ext.Data());
  cspectra->SaveAs(str.Data());
  
  // DRAWING RATIO OF SUCCESIVE ITERATIONS
  
  TCanvas *cratio = new TCanvas("cratio","cratio",10,10,can_x,can_y);
  cratio->cd();
  cratio->SetGrid();
  //cratio->SetLogy();
  frame->GetXaxis()->SetRangeUser(0, 30);
  frame->GetYaxis()->SetRangeUser(0.5, 1.5); 
  frame->GetXaxis()->SetTitle("p_{T, corr}^{charged} (GeV/c)");
  //frame->GetYaxis()->SetTitle("Ratio = distribution /truth");
  frame->GetYaxis()->SetTitle("Ratio = n / (n-1) iteration");
  frame->DrawCopy("");
  for(Int_t iter = firstiter+1; iter < lastiter; iter++)
    {if(skipiter(iter,SVD1))continue;
      //hunfolded[iter]->Divide(htruth);
      TH1D* hdiv=(TH1D*)hunfolded[iter]->Clone(Form("hdiv_%i",iter));
      hdiv->Divide(hunfolded[iter-1]);
      hdiv->DrawCopy("E same");
    }
  TLine *one = new TLine(0, 1, 30, 1);
  one->SetLineWidth(2);
  one->SetLineStyle(2);
  one->SetLineColor(kBlack);
  one->DrawClone("same");
  TLegend *legratio = new TLegend(0.6612, 0.7291, 0.8896, 0.90);
  legratio->SetFillStyle(0);
  legratio->SetBorderSize(0);
  for(Int_t iter = firstiter+1; iter < lastiter; iter++){
    if(skipiter(iter,SVD1))continue;
    legratio->AddEntry(hunfolded[iter], Form("Unfolded, N_{iter}=%i/N_{iter}=%i", iter+1, iter), "lp");
    }
  legratio->DrawClone("same");
  model_info->DrawClone("same");
  latex->DrawLatex(0.4, 0.75,"STAR Preliminary");
  str = Form("%s/%s_ratio_%s.%s", outdir.Data(),prefix.Data(),suffix.Data(),ext.Data());
  cratio->SaveAs(str.Data());
  
  
  // DRAWING RATIO UNFOLDED VS. TRUE (parton level)
  if(doToymodel && showHardJet){
  TCanvas *cratiot = new TCanvas("cratiot","ratio unf/true (PL)",10,10,can_x,can_y);
  cratiot->cd();
  cratiot->SetGrid();
  //cratiot->SetLogy();
  frame->GetXaxis()->SetRangeUser(0, 40);
  frame->GetYaxis()->SetRangeUser(0.1, 2); 
  frame->GetXaxis()->SetTitle("p_{T, corr}^{charged} (GeV/c)");
  //frame->GetYaxis()->SetTitle("Ratio = distribution /truth");
  //frame->GetYaxis()->SetTitle("Ratio = #frac{unfolded}{T_{AA}*d#sigma_{pp}/dp_{T}}");
  //frame->GetYaxis()->SetTitle("Ratio = #frac{unfolded}{true}");
  frame->GetYaxis()->SetTitle("unfolded/true (parton level)");
  frame->DrawCopy("");
  for(Int_t iter = firstiter; iter < lastiter; iter++)
    {
      if(skipiter(iter,SVD1))continue;
      //hunfolded[iter]->Divide(htruth);
      hratut[iter]=(TH1D*)hunfolded[iter]->Clone(Form("hratut_%i",iter));
      hratut[iter]->Divide(htruth);
      hratut[iter]->DrawCopy("E same");
    }
  TLine *one = new TLine(0, 1, 40, 1);
  one->SetLineWidth(2);
  one->SetLineStyle(2);
  one->SetLineColor(kBlack);
  one->DrawClone("same");
  
  TLegend *legrat = new TLegend(0.62, 0.7291, 0.8896, 0.90);
  legrat->SetFillStyle(0);
  legrat->SetBorderSize(0);
  for(Int_t iter = firstiter; iter < lastiter; iter++){
    if(skipiter(iter,SVD1))continue;
    legrat->AddEntry(hratut[iter], Form("Unfolded N_{iter}=%i vs. true", iter+1), "lp");
    }
  legrat->DrawClone("same");
  
  model_info->DrawClone("same");
  
  Double_t chi2t[niter];
  Double_t chi2t_trunk[niter];
  for(int iter=firstiter; iter<lastiter;iter++)
  {
  Int_t ndf=0;
  Int_t ndf_trunk=0;
   for(Int_t bin = 1; bin <= hratut[bfIter]->GetNbinsX(); bin++)
	{
		if(hratut[bfIter]->GetBinCenter(bin)>40)continue;
	  Double_t yld = hratut[bfIter]->GetBinContent(bin);
	  Double_t error=hratut[bfIter]->GetBinError(bin);
	  
	  if(yld>0){
	  chi2t[iter]+=(yld-1.0)*(yld-1.0)/(error*error);
	  ndf++;
	  
	  if(hratut[bfIter]->GetBinCenter(bin)<15) continue;
	  chi2t_trunk[iter]+=(yld-1.0)*(yld-1.0)/(error*error);
	   ndf_trunk++;
	  }
	}
  chi2t[iter]=chi2t[iter]/ndf;
  chi2t_trunk[iter]=chi2t_trunk[iter]/ndf_trunk;
  }
  
  latex->DrawLatex(0.4, 0.7,Form("#chi^{2}/NDF (%i itr) = %.2lf",chi2t[bfIter], bfIter+1));
  latex->DrawLatex(0.4, 0.6,Form("#chi^{2}/NDF (>15GeV) = %.2lf",chi2t_trunk[bfIter]));
  
  str = Form("%s/%s_ratioUnfVsTrue_%s.%s", outdir.Data(),prefix.Data(),suffix.Data(),ext.Data());
  cratiot->SaveAs(str.Data());
  
  //*******************************************
  //ratio unfolded/true (reconstructed jet level)
  //*******************************************
  TCanvas *cratiot2 = new TCanvas("cratiot2","ratio unf/true (JL)",10,10,can_x,can_y);
  cratiot2->cd();
  cratiot2->SetGrid();
  //cratiot->SetLogy();
  frame->GetXaxis()->SetRangeUser(0, 40);
  frame->GetYaxis()->SetRangeUser(0.1, 2); 
  frame->GetXaxis()->SetTitle("p_{T, corr}^{charged} (GeV/c)");
  //frame->GetYaxis()->SetTitle("Ratio = distribution /truth");
  //frame->GetYaxis()->SetTitle("Ratio = #frac{unfolded}{T_{AA}*d#sigma_{pp}/dp_{T}}");
  //frame->GetYaxis()->SetTitle("Ratio = #frac{unfolded}{true}");
  frame->GetYaxis()->SetTitle("unfolded/true (jet level)");
  frame->DrawCopy("");
  for(Int_t iter = firstiter; iter < lastiter; iter++)
    {
      if(skipiter(iter,SVD1))continue;
      //hunfolded[iter]->Divide(htruth);
      hratut2[iter]=(TH1D*)hunfolded[iter]->Clone(Form("hratut_%i",iter));
		hratut2[iter]->Divide(hjettruth);
		float stat_test=chi2_test(hunfolded[iter],hjettruth,30,5);
		cout<<"iteration: "<<iter+1<<" unfolded vs true (jet level) stat. test: "<<stat_test<<endl;
		/*
		hratut2[iter]->Reset("MICE");
		for(int bin=1;bin<hunfolded[iter]->GetNbinsX();bin++)
		{
			float val=hunfolded[iter]->GetBinContent(bin);
			if(!val>0)continue;
			float wdth=hunfolded[iter]->GetBinWidth(bin);
			float err=hunfolded[iter]->GetBinError(bin);
			float rerr=err/val; //relative error
			float x1=hunfolded[iter]->GetBinLowEdge(bin);
			float x2=hunfolded[iter]->GetBinLowEdge(bin+1);
			float intfjet=fjettrue->Integral(x1,x2); //integral of the charged PYTHIA over the bin
			float utrat=0;
			if(intfjet>0) utrat=val*wdth/intfjet; //ratio of bin integrals
			hratut2[iter]->SetBinContent(bin,utrat);
			hratut2[iter]->SetBinError(bin,utrat*rerr);
			
		}*/
      hratut2[iter]->DrawCopy("histo l same");
    }
  TLine *one = new TLine(0, 1, 40, 1);
  one->SetLineWidth(2);
  one->SetLineStyle(2);
  one->SetLineColor(kBlack);
  one->DrawClone("same");
  
  TLegend *legrat = new TLegend(0.62, 0.7291, 0.8896, 0.90);
  legrat->SetFillStyle(0);
  legrat->SetBorderSize(0);
  for(Int_t iter = firstiter; iter < lastiter; iter++){
    if(skipiter(iter,SVD1))continue;
    legrat->AddEntry(hratut2[iter], Form("Unfolded N_{iter}=%i vs. true", iter+1), "lp");
    }
  legrat->DrawClone("same");
  
  model_info->DrawClone("same");
  }
  
   // RAA
  if(showPP){
  TCanvas *craa = new TCanvas("craa","craa",10,10,can_x,can_y*1.3);
  craa->SetGrid();
  /*
  if(doToymodel){
  Double_t eps=0.02;
   TPad* p1 = new TPad("p1","p1",0,0.35-eps,1,1,0); p1->Draw();
   p1->SetBottomMargin(eps);
	p1->SetGrid();
	
	
	TPad* p2 = new TPad("p2","p2",0,0,1,0.35*(1.-eps),0); p2->Draw(); 
	p2->SetTopMargin(0);
	p2->SetBottomMargin(0.25);
	p2->SetGrid();
	p2->SetFillColor(0);
   p2->SetFillStyle(0);
	
	p1->cd();
  }*/
  //cratiot->SetLogy();
  frame->GetXaxis()->SetRangeUser(0, 40);
  frame->GetYaxis()->SetRangeUser(0, 1.4); 
  frame->GetXaxis()->SetTitle("p_{T, corr}^{charged} (GeV/c)");
  frame->GetYaxis()->SetTitle("R_{AA}");
  frame->DrawCopy("");
  
  TH1D* hrat=(TH1D*)hunfolded[compIter]->Clone("hrat");
  hrat->Divide(hppbase);
  hrat->DrawCopy("esame");
    /*
  TH1D* hrat2;
  if(doToymodel){
	 hrat2=(TH1D*)htruth->Clone("hrat2");
	 hrat2->Divide(hppbase);
    hrat2->DrawCopy("esame");
    }
    
    TH1D* hdiff;
    if(doToymodel){
      hdiff=(TH1D*)hrat->Clone("hdiff");
		hdiff->Reset("MICE");
		hdiff->SetTitle("");
		for(int bin=1; bin<=hdiff->GetNbinsX(); bin++){
			Double_t yld=(hrat->GetBinContent(bin))-(hrat2->GetBinContent(bin));
			//yld=yld+1.0;
			//yld=yld/hrat->GetBinContent(bin);
			Double_t err=(hrat->GetBinError(bin))*(hrat->GetBinError(bin))+(hrat2->GetBinError(bin))*(hrat2->GetBinError(bin));
			err=TMath::Sqrt(err);
			hdiff->SetBinContent(bin,yld);
			hdiff->SetBinError(bin,err);
		}
		hdiff->SetLineColor(kMagenta);
		hdiff->SetMarkerColor(kMagenta);
    }*/

  TLine *one = new TLine(0, 1, 40, 1);
  one->SetLineWidth(2);
  one->SetLineStyle(2);
  one->SetLineColor(kBlack);
  one->DrawClone("same");
  
  TLegend *legraa = new TLegend(0.62, 0.7291, 0.8896, 0.90);
  legraa->SetFillStyle(0);
  legraa->SetBorderSize(0);
  legraa->AddEntry(hunfolded[compIter], Form("R_{AA} (N_{iter}=%i)", compIter+1), "lp");
    //legraa->AddEntry(hrat[iter], Form("Unfolded N_{iter}=%i vs. Pythia", iter+1), "lp");
/*
    if(doToymodel){
		 legraa->AddEntry(hrat2,"R_{AA}^{true} ", "lp");
		 legraa->AddEntry(hdiff,"R_{AA} - R_{AA}^{true}", "lp");
	 }*/
	
  legraa->DrawClone("same");
    model_info->DrawClone("same");
	 
	/* if(doToymodel){
  p2->cd();
  frame2->GetXaxis()->SetRangeUser(0, 40);
  frame2->GetYaxis()->SetRangeUser(-0.85, 0.85);
  frame2->GetXaxis()->SetTitle("p_{T, corr}^{charged} (GeV/c)");
  frame2->GetYaxis()->SetTitle("");
  frame2->SetTitle("");
  frame2->GetYaxis()->SetLabelSize(0.06);
  frame2->GetXaxis()->SetLabelSize(0.06);
  //hdiff->GetYaxis()->SetTitleSize(0.06);
  frame2->GetXaxis()->SetTitleSize(0.09);
  frame2->DrawCopy("");
  hdiff->DrawCopy("esame");
  
  Double_t val=0;
  TLine *zer = new TLine(0, val, 40, val);
  zer->SetLineWidth(2);
  zer->SetLineStyle(2);
  zer->SetLineColor(kBlack);
  zer->DrawClone("same");
	 }*/
  str = Form("%s/%s_RAA_%s.%s", outdir.Data(),prefix.Data(),suffix.Data(),ext.Data());
  craa->SaveAs(str.Data());
  }
  /*
  //Ratio iterations vs one given iteration
  TCanvas *cratio2 = new TCanvas("cratio2","cratio2",10,10,can_x,can_y);
  cratio2->cd();
  cratio2->SetGrid();
  //cratio->SetLogy();
  frame->GetXaxis()->SetRangeUser(0, 25);
  frame->GetYaxis()->SetRangeUser(0.5, 1.5); 
  frame->GetXaxis()->SetTitle("p_{T, corr}^{charged} (GeV/c)");
  //frame->GetYaxis()->SetTitle("Ratio = distribution /truth");
  if(!SVD1)frame->GetYaxis()->SetTitle(Form("Ratio: n / %i iteration",compIter+1));
  else frame->GetYaxis()->SetTitle(Form("Ratio: kterm=n / kterm=%i ",compIter+1));
  frame->DrawCopy("");
  for(Int_t iter = firstiter; iter < lastiter; iter++)
    {if(skipiter(iter,SVD1) || iter==compIter)continue;
      //hunfolded[iter]->Divide(htruth);
      TH1D* hdiv=(TH1D*)hunfolded[iter]->Clone(Form("hdiv_%i",iter));
      hdiv->Divide(hunfolded[compIter]);
      hdiv->DrawCopy("E same");
    }
  TLine *uno = new TLine(0, 1, 25, 1);
  uno->SetLineWidth(2);
  uno->SetLineStyle(2);
  uno->SetLineColor(kBlack);
  uno->DrawClone("same");
  TLegend *legratio = new TLegend(0.6612, 0.7291, 0.8896, 0.90);
  legratio->SetFillStyle(0);
  legratio->SetBorderSize(0);
  for(Int_t iter = firstiter; iter < lastiter; iter++){
    if(skipiter(iter,SVD1) || iter==compIter)continue;
    if(!SVD1)legratio->AddEntry(hunfolded[iter], Form("Unfolded, N_{iter}=%i/N_{iter}=%i", iter+1, compIter+1), "lp");
    else legratio->AddEntry(hunfolded[iter], Form("Unfolded, kterm=%i/kterm=%i", iter+1, compIter+1), "lp");
    }
  legratio->DrawClone("same");
  model_info->DrawClone("same");
  latex->DrawLatex(0.4, 0.75,"STAR Preliminary");
  str = Form("%s/ratio/%s_ratio2_%s.%s", outdir.Data(),prefix.Data(),suffix.Data(),ext.Data());
  cratio2->SaveAs(str.Data());
  */
  // DRAWING RATIO BACKFOLDED VS. MEASURED
  if(showBF){
  TCanvas *cratiob = new TCanvas("cratiob","cratiob",10,10,can_x,can_y);
  cratiob->cd();
  cratiob->SetGrid();
  //cratiob->SetLogy();
  frame->GetXaxis()->SetRangeUser(-15, 40);
  frame->GetYaxis()->SetRangeUser(1E-1, 2); 
  frame->GetXaxis()->SetTitle("p_{T, corr}^{charged} (GeV/c)");
  //frame->GetYaxis()->SetTitle("Ratio = distribution /truth");
  frame->GetYaxis()->SetTitle("Ratio = #frac{backfolded}{measured}");
  frame->DrawCopy("");
  
  TH1D* hratbf[niter];
  for(Int_t iter = firstiter; iter < lastiter; iter++)
    {
      //if(iter!=bfIter)continue;
      //hunfolded[iter]->Divide(htruth);
      hratbf[iter]=(TH1D*)hbfold[iter]->Clone(Form("hbrat_%i",iter));
      hratbf[iter]->Divide(hmeasured);
      hratbf[iter]->DrawCopy("E same");
    }/*
  str = Form("iter%d", bfIter);
  dir = (TDirectoryFile*)funfolding->Get(str.Data());
  hratbf = (TH1D*)dir->Get("hbfmratio");
  hratbf->DrawCopy("E same");
  */
  
  //chi2/ndf calculation
  Double_t chi2[niter];
  Double_t chi2_trunk[niter];
  int trunk_cut=10;
  
  for(int iter=firstiter; iter<lastiter;iter++)
  {
  Int_t ndf=0;
  Int_t ndf_trunk=0;
  
  for(Int_t bin = 1; bin <= hratbf[iter]->GetNbinsX(); bin++)
	{
		if(hratbf[iter]->GetBinCenter(bin)>40)continue;
	  Double_t yld = hratbf[iter]->GetBinContent(bin);
	  Double_t error=hratbf[iter]->GetBinError(bin);
	 
	  if(yld>0){
	  chi2[iter]+=(yld-1.0)*(yld-1.0)/(error*error);
	  ndf++;
	  
	  if(hratbf[bfIter]->GetBinCenter(bin)<trunk_cut) continue;
	  chi2_trunk[iter]+=(yld-1.0)*(yld-1.0)/(error*error);
	   ndf_trunk++;
	  }
	}//bin loop
  chi2[iter]=chi2[iter]/ndf;
  chi2_trunk[iter]=chi2_trunk[iter]/ndf_trunk;
  }//iter loop
  
  latex->DrawLatex(0.4, 0.7,Form("#chi^{2}/NDF (%i itr) = %.2lf",chi2[bfIter], bfIter+1));
  latex->DrawLatex(0.4, 0.6,Form("#chi^{2}/NDF (>%iGeV) = %.2lf",trunk_cut,chi2_trunk[bfIter]));
  
  TLegend *legratb = new TLegend(0.62, 0.6, 0.8896, 0.90);
  legratb->SetFillStyle(0);
  legratb->SetBorderSize(0);
  for(Int_t iter = firstiter; iter < lastiter; iter++){
    //if(iter!=bfIter)continue;
    legratb->AddEntry(hbfold[iter], Form("Backfolded N_{iter}=%i vs. measured, chi2:%.1lf, chi2tr:%.1lf", iter+1,chi2[iter],chi2_trunk[iter] ), "lp");
    }
  legratb->DrawClone("same");
  
  model_info->DrawClone("same");
  
  TLine *one = new TLine(-15, 1, 40, 1);
  one->SetLineWidth(2);
  one->SetLineStyle(2);
  one->SetLineColor(kBlack);
  one->DrawClone("same");
  
  latex->DrawLatex(0.35, 0.3,"STAR Preliminary");
  
  str = Form("%s/%s_ratio_bfoldVSmeasured_%s.%s", outdir.Data(),prefix.Data(),suffix.Data(),ext.Data());
  cratiob->SaveAs(str.Data());
  
  }//if showBF
  }// !doCompare
  
    

  //--------------------------------
  //COMPARISON HISTOGRAMS
  //--------------------------------
  else{
   
  suffix=Form("%s_%s",suffix.Data(),comment2.Data()); 

  //MEASURED1,2 PRIOR1,2 UNFOLDED1,2 BACKFOLDED1,2 DISTRIBUTIONS
  TCanvas *spectra = new TCanvas("cspectra","cspectra",10,10,can_x,can_y);
  cspectra->cd();
  cspectra->SetGrid();
  cspectra->SetLogy();
  frame->GetXaxis()->SetTitle("p_{T, corr}^{charged} (GeV/c)");
  frame->GetYaxis()->SetTitle("1/N_{events} 1/2#pi d^{2}N/dp_{T}^{ch}d#eta (GeV/c)^{-1}");
  frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
  frame->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
  //frame->GetXaxis()->SetRangeUser(-30, 40);
  //frame->GetYaxis()->SetRangeUser(1E-7, 1E-1);
  str=Form("p_{T}^{leading} > %.1lf GeV/c, comparison: %s",pTthresh, desc2.Data());
  frame->SetTitle(str);
  frame->DrawCopy("");
  hmeasured->DrawCopy("E same");
  if(measured2)hmeasured2->DrawCopy("E same");
  
  hprior->GetXaxis()->SetRangeUser(pTthresh,priorRange);
  hprior->DrawCopy("HIST SAME C");
  //if(doRebin2)hprior2->Rebin(rebin2);
  
  hprior2->GetXaxis()->SetRangeUser(pTthresh,priorRange);
  if(!measured2)hprior2->DrawCopy("HIST SAME C");
  //htruth->DrawCopy("E same");
  if(showBF)hbfold[bfIter]->DrawCopy("E same");
  if(doRebin2){
	  hbfold2[bfIter2]->Rebin(rebin2);
      hbfold2[bfIter2]->Scale(1.0/rebin2);
  }
  if(showBF)hbfold2[bfIter2]->DrawCopy("E same");
  for(Int_t iter = firstiter; iter < lastiter; iter++){
			 if(skipiter(iter,SVD1))continue;
 			 if(iter!=bfIter)continue;

    hunfolded[iter]->DrawCopy("E same");
    }
  for(Int_t iter = firstiter2; iter < lastiter2; iter++){
			 if(skipiter(iter,SVD2))continue;
			 if(iter!=bfIter2)continue;
    if(doRebin2){
		hunfolded2[iter]->Rebin(rebin2);
		hunfolded2[iter]->Scale(1.0/rebin2);
	}
    hunfolded2[iter]->DrawCopy("E same");
    }
    
  //Info table
  TLegend *model_info = new TLegend(0.12037, 0.53, 0.3348, 0.900);
  model_info->SetFillStyle(0);
  model_info->SetBorderSize(0);
  model_info->SetMargin(0.05);
  model_info->SetHeader(Form("Run11 AuAu 200 GeV, %s",trigger.Data()));
  if(doToymodel)model_info->SetHeader("Toymodel: RHIC kinematics");
  //model_info->AddEntry("", "Single Particle", "");
  model_info->AddEntry("", "Charged jets", "");
  if(!peripheral)model_info->AddEntry("", "0-10% central collisions", "");
  else if(peripheral==1) model_info->AddEntry("", "60-80% central collisions", "");
  else if(peripheral==2) model_info->AddEntry("", "p+p collisions", "");
  model_info->AddEntry("", Form("N_{events} = %.1lfM", nevents/1E6), "");
  model_info->AddEntry("", Form("Anti-k_{T} \t \t R = %.1lf",R), "");
  model_info->AddEntry("", "p_{T}^{const} > 0.2 GeV/c", "");
  if(R>0.35) model_info->AddEntry("", "A_{reco jet} > 0.4 sr", "");
  else if(R>0.25) model_info->AddEntry("", "A_{reco jet} > 0.2 sr", "");
  else model_info->AddEntry("", "A_{reco jet} > 0.09 sr", "");
  //model_info->AddEntry("", "binning: 0.5 GeV", "");
  //model_info->AddEntry("", "removed: bad sector + 1*R edge", "");
  model_info->DrawClone("same");

  latex->DrawLatex(0.35, 0.3,"STAR Preliminary");
  
  TString m2type="STAR";
  if(doToymodel2)m2type="TOYMODEL";
   //Legend
  TLegend *legspectr = new TLegend(0.65, 0.45, 0.95, 0.90);
  legspectr->SetTextSize(0.03);
  legspectr->SetFillStyle(0);
  legspectr->SetBorderSize(0);
  legspectr->AddEntry(hmeasured, Form("measured, I=%.1lf*%s",integralm*iscale, iscale_str.Data()), "lp");
  //legspectr->AddEntry(hmeasured, Form("TOYMODEL, I=%.1lf*%s",integralm*iscale, iscale_str.Data()), "lp");
  if(measured2) legspectr->AddEntry(hmeasured2, Form("measured (%s), I=%.1lf*%s",m2type.Data(), integralm*iscale, iscale_str.Data()), "lp");
  if(unftype != unftype2)
    legspectr->AddEntry("",leg_name1,"");
  //legspectr->AddEntry(htruth,"T_{AA}*d#sigma_{pp}/dp_{T}", "lp");
  legspectr->AddEntry(hprior,Form("prior1: %s",prior_type_name[priorNo].Data()), "l");
  //legspectr->AddEntry(""," ", "");
  for(Int_t iter = firstiter; iter < lastiter; iter++){
			 if(skipiter(iter,SVD1))continue;
			 if(iter!=bfIter)continue;
    legspectr->AddEntry(hunfolded[iter], Form("unfolded, N_{iter}= %i",iter+1), "lp");
    //legspectr->AddEntry(""," ", "");
    }
  if(showBF)legspectr->AddEntry(hbfold[bfIter],Form("backfolded, N_{iter}= %i",bfIter+1), "lp");
  //legspectr->AddEntry(""," ", "");
  if(unftype != unftype2)
    legspectr->AddEntry("",leg_name2,"");
  if(!measured2){
    legspectr->AddEntry(hprior2,Form("prior2: %s",prior_type_name[priorNo2].Data()), "l");
    //legspectr->AddEntry(""," ", "");
  }
  for(Int_t iter = firstiter2; iter < lastiter2; iter++){
	if(skipiter(iter,SVD2))continue;
    if(iter!=bfIter2)continue;
    legspectr->AddEntry(hunfolded2[iter], Form("unfolded, N_{iter}= %i",iter+1), "lp");
    //legspectr->AddEntry(""," ", "");
  }
  if(showBF)legspectr->AddEntry(hbfold2[bfIter2],Form("backfolded, N_{iter}= %i",bfIter2+1), "lp");
  //legspectr->AddEntry(""," ", "");
  legspectr->DrawClone("same");
  
 
  
  str = Form("%s/compare/%s_spectra_bfiter%i-%i_%s.%s", outdir.Data(),prefix.Data(),bfIter+1,bfIter2+1,suffix.Data(),ext.Data());
  cspectra->SaveAs(str.Data());


  //RATIO BACKFOLDED1 VS BACKFOLDED2  
if(showBF){
  TCanvas *cratiobb = new TCanvas("cratiobb","cratiobb",10,10,can_x,can_y);
  cratiobb->cd();
  cratiobb->SetGrid();
  //cratiobb->SetLogy(); 
  frame->GetXaxis()->SetRangeUser(-15, 40);
  frame->GetYaxis()->SetRangeUser(5E-1, 1.5);
  frame->GetXaxis()->SetTitle("p_{T, corr}^{charged} (GeV/c)");
  frame->GetYaxis()->SetTitle("Ratio = #frac{backfolded1}{backfolded2}");
  frame->SetTitle(Form("Ratio of backfolded dist.: %s, N_{iter} = %i, %i",desc2.Data(),bfIter+1,bfIter2+1));
    if(bfIter==bfIter2) frame->SetTitle(Form("Ratio of backfolded dist.: %s, N_{iter} = %i",desc2.Data(),bfIter+1));
  frame->DrawCopy("");
  
  TH1D* hratbb=(TH1D*)hbfold[bfIter]->Clone(Form("hbbrat_%i",bfIter));
  hratbb->Divide(hbfold2[bfIter2]);
  hratbb->DrawCopy("E same");
  
  model_info->DrawClone("same");

  TLine *one = new TLine(-15, 1, 40, 1);
  one->SetLineWidth(2);
  one->SetLineStyle(2);
  one->SetLineColor(kBlack);
  one->DrawClone("same");

  str = Form("%s/compare/%s_ratio_backfolded_iter%i-%i_%s.%s", outdir.Data(),prefix.Data(),bfIter+1,bfIter2+1,suffix.Data(),ext.Data());
  cratiobb->SaveAs(str.Data());
}//if(showBF)
  //RATIO UNFOLDED1 VS UNFOLDED2
  TCanvas *cratiouu = new TCanvas("cratiouu","cratiouu",10,10,can_x,can_y);
  cratiouu->cd();
  cratiouu->SetGrid();
  //cratiouu->SetLogy();
  frame->GetXaxis()->SetRangeUser(0, 30); 
  frame->GetYaxis()->SetRangeUser(5E-1, 2);
  frame->GetXaxis()->SetTitle("p_{T, corr}^{charged} (GeV/c)");
  frame->GetYaxis()->SetTitle("Ratio = #frac{unfolded1}{unfolded2}");
  frame->SetTitle(Form("Ratio of unfolded dist.: %s, N_{iter} = %i, %i",desc2.Data(),compIter+1,compIter2+1));
  if(compIter==compIter2)frame->SetTitle(Form("Ratio of unfolded dist.: %s, N_{iter} = %i",desc2.Data(),compIter+1));
  frame->DrawCopy("");
  latex->DrawLatex(0.4, 0.8, Form("p_{T}^{leading}>%0.1lf GeV/c", pTthresh));
  latex->DrawLatex(0.4, 0.7,"STAR Preliminary");
  TH1D* hratuu=(TH1D*)hunfolded[compIter]->Clone(Form("huurat_%i",compIter));
  hratuu->Divide(hunfolded2[compIter2]);
  hratuu->DrawCopy("E same");

  TLine *one = new TLine(0, 1, 40, 1);
  one->SetLineWidth(2);
  one->SetLineStyle(2);
  one->SetLineColor(kBlack);
  one->DrawClone("same");
  
  model_info->DrawClone("same");
  
  str = Form("%s/compare/%s_ratio_unfolded_iter%i-%i_%s.%s", outdir.Data(),prefix.Data(),bfIter+1,bfIter2+1,suffix.Data(),ext.Data());
  cratiouu->SaveAs(str.Data());

  //RATIO BACKFOLDED1/Measured VS BACKFOLDED2/Measured
  
  if(showBF){
  TCanvas *cratiobmbm = new TCanvas("cratiobmbm","cratiobmbm",10,10,can_x,can_y);
  cratiobmbm->cd();
  cratiobmbm->SetGrid();
  //cratiobb->SetLogy(); 
  frame->GetXaxis()->SetRangeUser(-15, 40);
  frame->GetYaxis()->SetRangeUser(5E-1, 1.5);
  frame->GetXaxis()->SetTitle("p_{T, corr}^{charged} (GeV/c)");
  frame->GetYaxis()->SetTitle("Ratio = #frac{backfolded}{measured}");
  frame->SetTitle(Form("backfolded/measured, %s, p_{T}^{leading} > %.1lf GeV/c",desc2.Data(),pTthresh));
  frame->DrawCopy("");
  
  TH1D* hratbm1=(TH1D*)hbfold[bfIter]->Clone(Form("hbmrat_%i",bfIter));
  hratbm1->Divide(hmeasured);
  hratbm1->DrawCopy("E same");
  
  TH1D* hratbm2=(TH1D*)hbfold2[bfIter2]->Clone(Form("hbmrat2_%i",bfIter));
  if(measured2)hratbm2->Divide(hmeasured2);
  else hratbm2->Divide(hmeasured);
  hratbm2->DrawCopy("E same");
  
    //Legend
  TLegend *leg = new TLegend(0.6467, 0.75, 0.89, 0.90);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hratbm1, Form("backfolded1/measured, N_{iter}=%i",bfIter+1), "lp");
  leg->AddEntry(hratbm2, Form("backfolded2/measured, N_{iter}=%i",bfIter2+1), "lp");
  leg->DrawClone("same");
  
  model_info->DrawClone("same");

  TLine *one = new TLine(-15, 1, 40, 1);
  one->SetLineWidth(2);
  one->SetLineStyle(2);
  one->SetLineColor(kBlack);
  one->DrawClone("same");

  str = Form("%s/compare/%s_ratio_backfolded_measured12_iter%i-%i_%s.%s", outdir.Data(),prefix.Data(),bfIter+1,bfIter2+1,suffix.Data(),ext.Data());
  cratiobmbm->SaveAs(str.Data());
  }//if(showBF)
  
  
  // RAA1 VS RAA2 VS RAA3
  if(showPP){
  TCanvas *cratioraa = new TCanvas("cratioraa","cratioraa",10,10,can_x,can_y);
  cratioraa->cd();
  cratioraa->SetGrid();
  //cratiobb->SetLogy(); 
  frame->GetXaxis()->SetRangeUser(0, 40);
  frame->GetYaxis()->SetRangeUser(0.0, 1.5);
  frame->GetXaxis()->SetTitle("p_{T, corr}^{charged} (GeV/c)");
  frame->GetYaxis()->SetTitle("R_{AA}");
  if(compare3==3)frame->SetTitle(Form("R_{AA} comparison, %s",desc2.Data()));
  else frame->SetTitle(Form("R_{AA} comparison, %s, p_{T}^{leading} > %.1lf GeV/c",desc2.Data(),pTthresh));
  frame->DrawCopy("");
  
  TH1D* hraa1=(TH1D*)hunfolded[compIter]->Clone("hraa1");
  hraa1->Divide(hppbase);
  hraa1->DrawCopy("E same");
  
  TH1D* hraa2=(TH1D*)hunfolded2[compIter2]->Clone("hraa2");
  hraa2->Divide(hppbase);
  hraa2->DrawCopy("E same");
  
  TH1D* hraa3=(TH1D*)hunfolded3[compIter3]->Clone("hraa3");
  hraa3->Divide(hppbase);
  hraa3->SetLineColor(2);
  hraa3->SetMarkerColor(2);
  if(compare3>0)hraa3->DrawCopy("E same");
  
  /*
  for(int binx=1; binx<nbins1+1;binx++ )
  {
    if(hraa1->GetBinCenter(binx)<0 || hraa1->GetBinCenter(binx)>41) continue;
     cout<<binx<<",";
  }
  cout<<endl;
  cout<<endl;
  for(int binx=1; binx<nbins1+1;binx++ )
  {
        if(hraa1->GetBinCenter(binx)<0 || hraa1->GetBinCenter(binx)>41) continue;

   Double_t err1= hraa1->GetBinContent(binx) - hraa2->GetBinContent(binx);
   cout<<err1<<",";
   //cout<<hraa1->GetBinWidth(binx)/2<<",";
  }
  cout<<endl;
  cout<<endl;
  for(int binx=1; binx<nbins1+1;binx++ )
  {
        if(hraa1->GetBinCenter(binx)<0 || hraa1->GetBinCenter(binx)>41) continue;

   Double_t err1= hraa1->GetBinContent(binx) - hraa3->GetBinContent(binx);
   cout<<err1<<",";
  }
  cout<<endl;*/
  /*
  TH1D* hraa2=(TH1D*)hunfolded2[compIter2]->Clone("hraa2");
  hraa2->Divide(htruth2);
  hraa2->DrawCopy("E same");
  
  TH1D* hraa3=(TH1D*)hunfolded2[compIter2]->Clone("hraa3");
  hraa3->Divide(htruth3);
  hraa3->SetLineColor(2);
  hraa3->DrawCopy("E same");
  */
  
    //Legend
  TLegend *leg = new TLegend(0.6467, 0.75, 0.89, 0.90);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
/*
  leg->AddEntry(hraa1, Form("R_{AA}1, N_{iter}=%i",compIter+1), "lp");
  if(!SVD2)leg->AddEntry(hraa2, Form("R_{AA}2, N_{iter}=%i",compIter2+1), "lp");
  else leg->AddEntry(hraa2, Form("R_{AA}2, kterm=%i",compIter2+1), "lp");
*/
if(compare3==0){
  leg->AddEntry(hraa1, "R_{AA}1", "lp");
  leg->AddEntry(hraa2, "R_{AA}2", "lp");
  }
else if(compare3==1){
  leg->AddEntry(hraa1, "R_{AA}1, #epsilon", "lp");
  leg->AddEntry(hraa2, "R_{AA}2, #epsilon+5\%", "lp");
  leg->AddEntry(hraa3, "R_{AA}3, #epsilon-5\%", "lp");
}
else if(compare3==2){
  leg->AddEntry(hraa1, Form("R_{AA}1, R=%.1lf",R), "lp");
  leg->AddEntry(hraa2, Form("R_{AA}1, R=%.1lf",R2), "lp");
  leg->AddEntry(hraa3, Form("R_{AA}1, R=%.1lf",R3), "lp");
}
else if(compare3==3){
  leg->AddEntry(hraa1, Form("R_{AA}1, p_{T}^{leading}>%.1lf",pTthresh), "lp");
  leg->AddEntry(hraa2, Form("R_{AA}1, p_{T}^{leading}>%.1lf",pTthresh2), "lp");
  leg->AddEntry(hraa3, Form("R_{AA}1, p_{T}^{leading}>%.1lf",pTthresh3), "lp");
}
else if(compare3==4){
  leg->AddEntry(hraa1, Form("R_{AA}1 (%s)",prior_type_name[priorNo].Data()), "lp");
  leg->AddEntry(hraa2, Form("R_{AA}2, (%s)",prior_type_name[priorNo2].Data()), "lp");
  leg->AddEntry(hraa3, Form("R_{AA}3, (%s)",prior_type_name[priorNo3].Data()), "lp");
}
else if(compare3==5){
  leg->AddEntry(hraa1, Form("R_{AA} N_{iter}=%i",compIter+1), "lp");
  leg->AddEntry(hraa2, Form("R_{AA} N_{iter}=%i",compIter2+1), "lp");
  leg->AddEntry(hraa3, Form("R_{AA} N_{iter}=%i",compIter3+1), "lp");
  }
/*
  leg->AddEntry(hraa1, Form("R_{AA}1, T_{AA}=%.1lf",TAA), "lp");
  leg->AddEntry(hraa2, Form("R_{AA}2, T_{AA}=%.1lf",TAA+1), "lp");
  leg->AddEntry(hraa3, Form("R_{AA}3, T_{AA}=%.1lf",TAA-1), "lp");
  */
  leg->DrawClone("same");
  
  model_info->DrawClone("same");

  TLine *one = new TLine(0, 1, 40, 1);
  one->SetLineWidth(2);
  one->SetLineStyle(2);
  one->SetLineColor(kBlack);
  one->DrawClone("same");

  str = Form("%s/compare/%s_RAAcomp_iter%i-%i_%s.%s", outdir.Data(),prefix.Data(),bfIter+1,bfIter2+1,suffix.Data(),ext.Data());
  cratioraa->SaveAs(str.Data());
  }//if(showPP)
  
 }//doCompare

 
 // DRAWING CORRELATION MATRIX
 /*if(!SVD1 && !Roo1){
   
   for(Int_t iter = firstiter; iter < lastiter; iter++)
    {
      if(skipiter(iter,SVD1))continue;
     str=Form("corr_iter%i",iter);
    TCanvas *ccorr=new TCanvas(str,str,10,10,can_x,can_y);
    ccorr->cd();
    ccorr->SetGrid();
  
      //if(iter!=bfIter)continue;
      
      //str = Form("ccor%d", iter);
      frame->GetXaxis()->SetRangeUser(0, 40);
      frame->GetYaxis()->SetRangeUser(0, 40);
      hcorrelation[iter]->GetZaxis()->SetRangeUser(-1, +1);
      frame->GetXaxis()->SetTitle("p_{T}^{unfolded} (GeV/c)");
      frame->GetYaxis()->SetTitle("p_{T}^{unfolded} (GeV/c)");
      frame->SetTitle("Correlation Matrix");
      frame->DrawCopy("");
      hcorrelation[iter]->DrawClone("COLZ SAME");
      latex->SetTextSize(0.04);
      latex->DrawLatex(0.15, 0.8, Form("Correlation Matrix - %s", errtype.Data()));  
      latex->DrawLatex(0.15, 0.7, "#rho(x_{i}, x_{j}) = #frac{cov(x_{i}, x_{j})}{#sigma_{i}#sigma_{j}}");
      latex->DrawLatex(0.15, 0.6, Form("N_{iterations} = %d", iter + 1));
      latex->DrawLatex(0.15, 0.5, Form("N_{events} = %.1lfM", nevents/1E6));
      latex->DrawLatex(0.15, 0.4,"STAR Preliminary");
      str = Form("%s/corr/%s_corrmatrix_iter%i_%s.%s", outdir.Data(),prefix.Data(), iter, suffix.Data(),ext.Data());
      ccorr->SaveAs(str.Data());
      
    }
      
   
}*/
 

      //Covariance Matrix      
  TCanvas *ccov =new TCanvas("ccov","ccov",10,10,can_x,can_y);
  ccov->cd();
  ccov->SetGrid();
  ccov->SetLogz();
  frame->GetXaxis()->SetRangeUser(-20, 60);
  frame->GetYaxis()->SetRangeUser(0, 60);
  frame->GetXaxis()->SetTitle("p_{T}^{unfolded} (GeV/c)");
  frame->GetYaxis()->SetTitle("p_{T}^{unfolded} (GeV/c)");
  frame->SetTitle("Covariance Matrix");
  frame->DrawCopy("");
  hcovariance[bfIter]->Draw("COLZ SAME");
  latex->DrawLatex(0.12, 0.8, "Covariance Matrix");
  latex->DrawLatex(0.12, 0.7, Form("N_{iterations} = %d", iter + 1));
  latex->DrawLatex(0.12, 0.6, Form("N_{events} = %.1lfM", nevents/1E6));
str = Form("%s/%s_CovMatrix_iter%i_%s.%s", outdir.Data(),prefix.Data(), iter, suffix.Data(),ext.Data());
ccov->SaveAs(str.Data());


  if(SVD1){
//D-vector
	TCanvas *cdv =new TCanvas("cdv","cdv",10,10,can_x,can_y);
  cdv->cd();
  cdv->SetGrid();
  cdv->SetLogy();
  frame->GetXaxis()->SetRangeUser(0, 100);
  frame->GetYaxis()->SetRangeUser(1E-2, 1E3);
  frame->GetXaxis()->SetTitle("k-term");
  frame->GetYaxis()->SetTitle("D");
  frame->SetTitle("D-vector");
  frame->DrawCopy("");
  hdvec[bfIter]->Draw("SAME");
  
str = Form("%s/%s_Dvector_kterm%i_%s.%s", outdir.Data(),prefix.Data(), iter, suffix.Data(),ext.Data());
cdv->SaveAs(str.Data());
 }
 
//Response Matrix      
  TCanvas *crm =new TCanvas("crm","crm",10,10,can_x,can_y);
  crm->cd();
  crm->SetGrid();
  crm->SetLogz();
  frame->GetXaxis()->SetRangeUser(-30, 60);
  frame->GetYaxis()->SetRangeUser(-30, 60);
  frame->GetXaxis()->SetTitle("p_{T, measured}^{charged} (GeV/c)");
  frame->GetYaxis()->SetTitle("p_{T, true}^{charged} (GeV/c)");
  frame->SetTitle("Response Matrix");
  frame->DrawCopy("");
hresponse->Draw("COLZ SAME");
latex->DrawLatex(0.25, 0.7,"STAR Preliminary");
str = Form("%s/Rmatrix_%s.%s", outdir.Data(),suffix.Data(),ext.Data());
crm->SaveAs(str.Data());
}
