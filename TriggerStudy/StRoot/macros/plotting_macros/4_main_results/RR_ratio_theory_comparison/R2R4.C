#include "../../utils/common.h"

void R2R4(float pTlead=5.0,int binning=1, TString evo="GPC2",TString label="", TString ext="pdf"){
	
	const int ncent=2; //number of centrality bins (= number of panels to be plotted)
	//what do we want to show
	bool showKrishnaAll=1; //show all 3 variants of hybrid model calculation
	bool showVitevAuAu=1;
	bool showVitevpp=0;
	bool showSoyez=0;
	bool showKrishna=1;
		bool wake=1; //what effects to include in the hybrid model
		bool wake_positive_only=0; //what effects to include in the hybrid model
	bool showHerwig=1;
	bool showPythia=1;
	bool biasedPP=0; //use pTlead cut also for PYTHIA and HERWIG?
	bool plotLogY=1; //log y-axis
	if(showKrishnaAll)
	{
		showKrishna=0;
	}

  float ptlead=pTlead;
  int Jan_binning=0;
  const int ngrp=13;//used for Jan's binning
  double xbins[ngrp+1]={5,6,7,8,9,10,12,14,16,18,20,25,30,40};

  gStyle->SetLineStyleString(12,"[40 20]");
  gStyle->SetLineStyleString(13,"[13 13]");
  gStyle->SetLineStyleString(14,"[40 13 13 13]");

  TString outDir=Form("../../../plotting_out/obr/%s/results/comparison",evo.Data());
  
  TString cent_dir[ncent]= {"peripheral","central"};
  TString cent_leg[ncent]= {"peripheral (60-80%)", "central (0-10%)" };
  
  //clors
  Color_t cPy=kRed; //PYTHIA
  Color_t cHw=kRed; //Herwig
  Color_t cNLO=kRed-2; //Soyez NLO
  Color_t cSTAR=kBlue; //STAR data
  Color_t cSTAR2=kBlue-7; //STAR data - errorbars
  Color_t cVit1=kAzure-7; //Ivan's NLO
  Color_t cVit2=kCyan-2; //Ivan's SCET model
  Color_t ckrsh=kMagenta+2; //Krishna's HYBRID model
	Color_t ckrsh1=kMagenta+2; //Krishna's HYBRID model
	Color_t ckrsh2=kViolet-9; //Krishna's HYBRID model
	Color_t ckrsh3=kOrange-8; //Krishna's HYBRID model
  
  //line width
  int lwPy=2;
  int lwHw=2;
  int lwNLO=2;
  int lwSTAR=2;
  int lwSTAR_sys=1;
  int lwVit1=2;
  int lwVit2=2;
  
  //line style
  int lsPy=14;
  int lsHw=2;
  int lsNLO=3;
  int lsSTAR=1;
  int lsVit1=1;
  int lsVit2=2;
  
  //marker style and size
  int msSTAR=29;
  float msizeSTAR=3;
  
  //fill style
  int vit1Fill=1001; //Ivans's RRratio fill style
  int vit2Fill=1001; //Ivans's RRratio fill style
  int krshFill=1001; //Krishna's RRratio fill style
	int krshFill1=1001; //Krishna's RRratio fill style
	int krshFill2=1001; //Krishna's RRratio fill style
	int krshFill3=1001; //Krishna's RRratio fill style
  
  	//line at unity
	Color_t lineOneColor=kGray+1;
	int lineOneWidth=1;
	int lineOneStyle=2;
  
  //label
	float G_latex_sz=0.032;
	Color_t G_latex_cl=kGray+1;
  
	//Axis ranges
	float xMin=0;//x-axis range
	float xMax[ncent]={41,41}; //x-axis range
	float ptMax[ncent]={25,30}; //STAR data x-range
	float yMin=(plotLogY)? 0.05:0.05;
	float yMax=(plotLogY)?2.0:1.3;
	
	//biased region mark up
	float xbiasLine=16; 
	float ybiasLine_min=0.3; 
	float ybiasLine_max=1.6;
	float xbiasDesc=0.53;
	float ybiasDesc=0.85;
	Color_t biasColor=kPink-9;
	
	//****************
	//load data
	//****************
	
  //---------------------------------------------------------------------------
  //Soyez paper
  float Xsoyez_full[10]={12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5};
  float Ysoyez_full[10]={0.297,0.353,0.402,0.426,0.446,0.446,0.452,0.446,0.437,0.399};
  float Esoyez_full[10]={0.20,0.16,0.13,0.13,0.12,0.12,0.11,0.11,0.12,0.13};
  TGraph *gsoyez=graph_band(10,Xsoyez_full,Ysoyez_full,Esoyez_full,xMax[1]-2);

  //---------------------------------------------------------------------------
  //Vitev
  //gROOT->LoadMacro("./util.C");
  TGraphErrors *gppR2R4=gReadData("theory/jet_spectra/Vitev/PPXsecRatio_200GeV_R2overR4.txt",1,0,xMax[1]-2);
  TGraphErrors *gAAR2R4dn=gReadData("theory/jet_spectra/Vitev/AuAuXsecRatio_smallCNM_200GeV_DN_R2overR4.txt",1,0,xMax[1]-2);
  TGraphErrors *gAAR2R4up=gReadData("theory/jet_spectra/Vitev/AuAuXsecRatio_smallCNM_200GeV_UP_R2overR4.txt",1,0,xMax[1]-2);

  TGraphErrors *gppR2R4_dEdx=gReadData("theory/jet_spectra/Vitev/RRpp_02to04_200GeVNLO.dat",1,0,xMax[1]-2);
  TGraphErrors *gAAR2R4up_NLO=gReadData("theory/jet_spectra/Vitev/RR_NLO-Au0.2to0.4coldUP.dat",1,0,xMax[1]-2);
  TGraphErrors *gAAR2R4dn_NLO=gReadData("theory/jet_spectra/Vitev/RR_NLO-Au0.2to0.4coldDN.dat",1,0,xMax[1]-2);

  //---------------------------------------------------------------------------
  //Krishna
  //TGraphErrors *gKrishnaR2R4=gReadData3(Form("theory/hybrid_results_2017/cent010jetOver02_R04.dat"),5,2.0);
  TGraphErrors *gKrishnaR2R4=gReadData4(Form("theory/hybrid_results_2019/010_SpecRatio_R2_over_R4_wake_%i_ignore_neg_%i.dat",wake,wake_positive_only),5,2.0);
  
  TGraphErrors *gKrishnaR2R4_1=gReadData4("theory/hybrid_results_2019/010_SpecRatio_R2_over_R4_wake_0_ignore_neg_0.dat",5,2.0);
  TGraphErrors *gKrishnaR2R4_2=gReadData4("theory/hybrid_results_2019/010_SpecRatio_R2_over_R4_wake_1_ignore_neg_1.dat",5,2.0);
  TGraphErrors *gKrishnaR2R4_3=gReadData4("theory/hybrid_results_2019/010_SpecRatio_R2_over_R4_wake_1_ignore_neg_0.dat",5,2.0);
  
  
  //===========================================================================
  //Herwig, Pythia
  float binsize=1;
  float scale=1000;//mb->mub
  TH1D *h_ptPy_R02=getPythiaHisto("theory/jet_spectra/pythia_chrgEta1/histos_pythiajet_R0.2.root",scale*1e-6,binsize,(biasedPP) ? ptlead:0);//1e6 is nevents
  TH1D *h_ptPy_R04=getPythiaHisto("theory/jet_spectra/pythia_chrgEta1/histos_pythiajet_R0.4.root",scale*1e-6,binsize,(biasedPP) ? ptlead:0);

  TH1D *h_ptHw_R02=getPythiaHisto("theory/jet_spectra/herwig_chrgEta1/starjet_herwig_R0.2_chrg.root",scale,binsize,(biasedPP) ? ptlead:0);//Herwig nevents has already been divided
  TH1D *h_ptHw_R04=getPythiaHisto("theory/jet_spectra/herwig_chrgEta1/starjet_herwig_R0.4_chrg.root",scale,binsize,(biasedPP) ? ptlead:0);

  //---------------------------------------------------------------------------
  //NLO theory curve
  //gROOT->LoadMacro("./UtilityGraph.C");
  TGraphErrors *gNLO2=gReadData("theory/jet_spectra/Vogelsang/star_ct14_200_05_R02_sc1.dat",scale*1e-9,0);
  TGraphErrors *gNLO4=gReadData("theory/jet_spectra/Vogelsang/star_ct14_200_05_R04_sc1.dat",scale*1e-9,0);
  TGraphErrors *gNLO=Gdiv(gNLO2,gNLO4);


  //---------------------------------------------------------------------------
  if(Jan_binning){
    h_ptPy_R02->Rebin(ngrp,"hptPy_R02",xbins);
    h_ptHw_R02->Rebin(ngrp,"hptHw_R02",xbins);

    h_ptPy_R04->Rebin(ngrp,"hptPy_R04",xbins);
    h_ptHw_R04->Rebin(ngrp,"hptHw_R04",xbins);

 
    //take care of bin size
    hptPy_R02->Scale(h_ptPy_R02->GetBinWidth(1),"width");
    hptHw_R02->Scale(h_ptHw_R02->GetBinWidth(1),"width");

    hptPy_R04->Scale(h_ptPy_R04->GetBinWidth(1),"width");
    hptHw_R04->Scale(h_ptHw_R04->GetBinWidth(1),"width");


  }else{
    TH1D *hptPy_R02=(TH1D*)h_ptPy_R02->Clone();
    TH1D *hptHw_R02=(TH1D*)h_ptHw_R02->Clone();

    TH1D *hptPy_R04=(TH1D*)h_ptPy_R04->Clone();
    TH1D *hptHw_R04=(TH1D*)h_ptHw_R04->Clone();

 
  }
  hptPy_R02->SetName("hptPy_R02");
  hptHw_R02->SetName("hptHw_R02");
  hptPy_R04->SetName("hptPy_R04");
  hptHw_R04->SetName("hptHw_R04");

  //ratio R=0.2/R=0.4
  TH1D *hptPy=(TH1D*)hptPy_R02->Clone();hptPy->Divide(hptPy_R04);
  TH1D *hptHw=(TH1D*)hptHw_R02->Clone();hptHw->Divide(hptHw_R04);
	

  //-----------
  //STAR data
  //-----------
	TFile *fp[ncent];

  	TGraphAsymmErrors *gjan_unf[ncent];
	TGraphAsymmErrors *gjan_corr[ncent];
	TGraphAsymmErrors *gjan[ncent];
	TGraphAsymmErrors *gjan_sys[ncent];
	
	for(int cent=0; cent<2; cent++)
	{
		TString fname=Form("../../../plotting_out/systematics/%s/%s/bining%i/systematics_R0.4_pT%0.lf.root",cent_dir[cent].Data(),evo.Data(),binning,pTlead);
		fp[cent]=new TFile(fname,"read");
		cout<<"opening file: "<<fname.Data()<<endl;
   
	
		gjan_unf[cent]=(TGraphAsymmErrors*) fp[cent]->Get("RR_unfold_BSL");
		gjan_corr[cent]=(TGraphAsymmErrors*) fp[cent]->Get("RR_corr_BSL");
		gjan[cent]=(TGraphAsymmErrors*) fp[cent]->Get("RR_BSL");
		gjan_sys[cent]=gcombine_syserr(gjan_unf[cent],gjan_corr[cent],1);

		//ShrinkGraph(gjan_unf[cent],pTlead,ptMax[cent]);
		//ShrinkGraph(gjan_corr[cent],pTlead,ptMax[cent]);
		ShrinkGraph(gjan[cent],pTlead,ptMax[cent]);
		ShrinkGraph(gjan_sys[cent],pTlead,ptMax[cent]);

		//gjan_sys[cent]=(TGraphAsymmErrors*) fp[cent]->Get("RR_unfold_BSL");
	}

		//Soyez NLO
		gsoyez->SetLineColor(cNLO);
		gsoyez->SetLineWidth(lwNLO);
		gsoyez->SetMarkerColor(cNLO);
		gsoyez->SetLineStyle(lsNLO);
		
		//Vitev AuAu
		TGraphErrors* gvit1=area_graph(gAAR2R4up_NLO,gAAR2R4dn_NLO); //NLO pQCD
		gvit1->SetLineColor(cVit1);
		//gvit1->SetLineStyle(lsVit1);
		//gvit1->SetLineWidth(lwVit1);
		gvit1->SetFillStyle(vit1Fill);
		gvit1->SetFillColor(cVit1);
		
		TGraphErrors* gvit2=area_graph(gAAR2R4up,gAAR2R4dn); //SCET
		gvit2->SetLineColor(cVit2);
		//gvit2->SetLineStyle(lsVit2);
		//gvit2->SetLineWidth(lwVit2);
		gvit2->SetFillStyle(vit2Fill);
		gvit2->SetFillColor(cVit2);
		
		//Krishna
		gKrishnaR2R4->SetLineColor(ckrsh);
		gKrishnaR2R4->SetFillStyle(krshFill);
		gKrishnaR2R4->SetFillColor(ckrsh);
		
		gKrishnaR2R4_1->SetLineColor(ckrsh1);
		gKrishnaR2R4_1->SetFillStyle(krshFill1);
		gKrishnaR2R4_1->SetFillColor(ckrsh1);
		
		gKrishnaR2R4_2->SetLineColor(ckrsh2);
		gKrishnaR2R4_2->SetFillStyle(krshFill2);
		gKrishnaR2R4_2->SetFillColor(ckrsh2);
		
		gKrishnaR2R4_3->SetLineColor(ckrsh3);
		gKrishnaR2R4_3->SetFillStyle(krshFill3);
		gKrishnaR2R4_3->SetFillColor(ckrsh3);


  //===========================================================================
  //Draw plot
  //===========================================================================

	//gROOT->LoadMacro("./Utility.C");
	gStyle->SetOptStat(0);
	
	TLatex *latexL = new TLatex();
	latexL->SetNDC();
	latexL->SetTextSize(G_latex_sz);
	latexL->SetTextColor(G_latex_cl);
  
	TCanvas *cfig1=new TCanvas("cfig1","fig1",10,10,ncent*600,600);
	cfig1->Divide(ncent,1,0,0);
  
	TH1D *htmp1[ncent];
	TLegend *legraa0;
	TLegend *legraa1;
	
	for(int cent=0; cent<ncent; cent++)
	{
		float ytxt=-1,dytxt=0.055,tsize=0.035;

		cfig1->cd(cent+1);
		gPad->SetMargin(0.15,0.03,0.20,0.03);
		if(cent==0)gPad->SetRightMargin(0.0);
		if(cent==1)
		{
			gPad->SetLeftMargin(0.0);
			xMin=xMin=0.1;
		}
		gPad->SetTicks(1);
		if(plotLogY)gPad->SetLogy();
  
		htmp1[cent]=new TH1D(Form("htmp%i",cent),"",100,xMin,xMax[cent]);
		htmp1[cent]->SetXTitle("p_{T, jet}  or  p_{T, jet}^{ch}  (GeV/#it{c})");
		if(cent==0)htmp1[cent]->SetYTitle("R=0.2/R=0.4");
		htmp1[cent]->SetTitleOffset(1.2,"x");htmp1[cent]->SetTitleOffset(1.1,"y");
		htmp1[cent]->SetTitleSize(0.06,"x");htmp1[cent]->SetTitleSize(0.06,"y");
		htmp1[cent]->SetLabelSize(0.05,"x");htmp1[cent]->SetLabelSize(0.05,"y");
		//htmp1[cent]->SetNdivisions(505,"x");htmp1[cent]->SetNdivisions(505,"y");
		htmp1[cent]->SetMinimum(yMin);htmp1[cent]->SetMaximum(yMax);
		htmp1[cent]->DrawClone();
		
		TLine *one = new TLine(xMin, 1, xMax[cent],1);
		one->SetLineWidth(lineOneWidth);
		one->SetLineStyle(lineOneStyle);
		one->SetLineColor(lineOneColor);
		one->DrawClone("same");

		TLatex *ltx=new TLatex();ltx->SetNDC(1);ltx->SetTextSize(tsize);
		float xtex=0.2;
		if(cent) xtex=0.1;
		ltx->SetTextColor(1);ltx->DrawLatex(xtex,0.35,cent_leg[cent]);
		
		//Pythia and Herwig
		if(showHerwig) plot_h(hptHw,cHw,lsHw,lwHw,"",0.55,ytxt,tsize,0,pTlead,xMax[cent]-2);ytxt-=dytxt;
		if(showPythia) plot_h(hptPy,cPy,lsPy,lwPy,"",0.55,ytxt,tsize,0,pTlead,xMax[cent]-2);ytxt-=dytxt;

		//Soyez NLO
		if(showSoyez && cent==1) gsoyez->Draw("l");
		
		//Vitev AuAu
		if(showVitevAuAu && cent==1) gvit1->Draw("3");
		if(showVitevAuAu && cent==1) gvit2->Draw("3");
		
		//Krishna
		if(showKrishna && cent==1) gKrishnaR2R4->Draw("3");
		if(showKrishnaAll && cent==1) 
		{
			gKrishnaR2R4_1->Draw("3");
			gKrishnaR2R4_2->Draw("3");
			gKrishnaR2R4_3->Draw("3");
		}
	
		//STAR data
		gjan_sys[cent]->SetLineColor(cSTAR);
		gjan_sys[cent]->SetFillStyle(0);
		gjan_sys[cent]->SetLineWidth(lwSTAR_sys);
		gjan_sys[cent]->DrawClone("2");
		gjan[cent]->SetLineWidth(lwSTAR);
		gjan[cent]->SetLineColor(cSTAR2);
		gjan[cent]->SetMarkerColor(cSTAR);
		gjan[cent]->SetMarkerStyle(msSTAR);
		gjan[cent]->SetMarkerSize(msizeSTAR);
		gjan[cent]->DrawClone("p");
		
		
		//biased region - this one is only for legend
		TGraph *gbias=new TGraph();
		gbias->SetLineWidth(lwSTAR);
		gbias->SetLineColor(kGray);
		gbias->SetMarkerColor(kGray+1);
		gbias->SetMarkerStyle(msSTAR);
		gbias->SetMarkerSize(msizeSTAR);
		
		//mark up unbiased region
		TLine *bias = new TLine(xbiasLine,ybiasLine_min,xbiasLine,ybiasLine_max);
		bias->SetLineWidth(2);
		bias->SetLineStyle(2);
		bias->SetLineColor(biasColor);
		bias->DrawClone("same");
      TLatex *latexbias = new TLatex();
		latexbias->SetNDC();
		latexbias->SetTextSize(0.035);
		latexbias->SetTextColor(biasColor);
		float xbdesc=xbiasDesc;
		if(cent) xbdesc=xbiasDesc-0.1;
		latexbias->DrawLatex(xbdesc, ybiasDesc,"--> ~ UNBIASED");


		//legend
		double posx1,posx2,posy1,posy2;		
		if(cent==0)
		{
		posx1=0.56;
		posx2=0.86;
		posy1=(plotLogY) ? 0.25:0.58;
		posy2=(plotLogY) ? 0.62:0.88;
		legraa0 = new TLegend(posx1,posy1,posx2,posy2);
		legraa0->SetTextSize(0.036);
		legraa0->SetFillStyle(0);
		legraa0->SetBorderSize(0);
		legraa0->SetHeader("Charged jets, anti-k_{T}");
		legraa0->AddEntry("","Au+Au #sqrt{s_{NN}}= 200 GeV","");
		//	legraa0->AddEntry(gjan[cent], Form("  STAR, p_{T,lead}^{min}= %.1lf GeV/#it{c}",pTlead),"lp"); //Jana
		legraa0->AddEntry(gjan[cent], Form("  STAR, p_{T,lead}^{min}= %.0lf GeV/#it{c}",pTlead),"lp"); //Jana
		legraa0->AddEntry(gbias, "  STAR, biased region","lp");
		legraa0->AddEntry("","p+p #sqrt{s}=200GeV","");
		if(showHerwig) 
		{
			legraa0->AddEntry(hptHw, "  Herwig","lp");
			if (biasedPP) legraa->AddEntry("", Form("    p_{T,lead}^{min}= %.1lf GeV/#it{c}",pTlead),"");
		}
		//legraa->AddEntry("",Form("    p_{T}^{lead}>%.1lf",pTlead),"");
		if(showPythia) 
		{
			legraa0->AddEntry(hptPy,"  PYTHIA","lp");
			if (biasedPP)  legraa->AddEntry("", Form("    p_{T,lead}^{min}= %.1lf GeV/#it{c}",pTlead),"");
		}
		//legraa->AddEntry("",Form("    p_{T}^{lead}>%.1lf",pTlead),"");
		if(showSoyez) legraa0->AddEntry(gsoyez, "  Soyez NLO (full jets)","lp");
		cfig1->cd(1);
		legraa0->DrawClone("same");
		}

		//legend
		if(cent==1)
		{
			posx1=0.45;
			posx2=0.75;
			posy1=(plotLogY) ?0.25 : 0.55;  //gppR2R4->SetLineColor(cpp1);gppR2R4->SetLineStyle(14);gppR2R4->Draw("l");keyLine(0.55,ytxt,"Vitev pp (Model 2)",cpp1,14,tsize);
			posy2=(plotLogY) ? 0.55: 0.95;
		legraa1 = new TLegend(posx1,posy1,posx2,posy2);
		legraa1->SetTextSize(0.036);
		legraa1->SetFillStyle(0);
		legraa1->SetBorderSize(0);
		//legraa->AddEntry("",Form("    p_{T}^{lead}>%.1lf",pTlead),"");
		legraa1->SetHeader("Au+Au #sqrt{s_{NN}}= 200 GeV");
		if(showKrishna) legraa1->AddEntry(gKrishnaR2R4, "  Hybrid model, ch. jets", "f");
			if(showKrishnaAll)
			{
				legraa1->AddEntry("","Hybrid model, ch. jets:","");
				legraa1->AddEntry(gKrishnaR2R4_1, "    no medium resp.", "f");
				legraa1->AddEntry(gKrishnaR2R4_2, "    pos. resp. from wake ", "f");
				legraa1->AddEntry(gKrishnaR2R4_3, "    full medium resp.", "f");
			}
			legraa1->AddEntry("","Full jets:","");
			if(showVitevAuAu) legraa1->AddEntry(gvit1/*gAAR2R4up_NLO*/, "  NLO pQCD","f");
			if(showVitevAuAu) legraa1->AddEntry(gvit2/*gAAR2R4up*/, "  SCET","f");
			
		cfig1->cd(2);
		legraa1->DrawClone("same");
		latexL->DrawLatex(0.2,0.35,label);
		legraa1->DrawClone("same");
		
		}
	
	//biased region in gray
	ShrinkGraph(gjan[cent],pTlead,xbiasLine);
	ShrinkGraph(gjan_sys[cent],pTlead,xbiasLine);
	gjan[cent]->SetLineColor(kGray);
	gjan[cent]->SetMarkerColor(kGray+1);
	gjan_sys[cent]->SetLineColor(kGray+1);
	gjan[cent]->DrawClone("p");
	gjan_sys[cent]->DrawClone("2");
	
		
		
	}//centrality class
	TString KrishnaSuff="";
	//if(showKrishnaAll)KrishnaSuff="_HybridAll";
	cfig1->SaveAs(Form("%s/R2R4_pTl%.0lf_bin%i%s.%s",outDir.Data(),pTlead,binning,KrishnaSuff.Data(),ext.Data()));
	return;

}
//===========================================================================
//ADDITIONAL FUNCTIONS
//===========================================================================
TGraphErrors *gReadData(char *file="star_6m_200_0208_sc1.dat",float scale=1,float err=0,float xmax=10){
  FILE *fp=fopen(file,"read");
  float x[1000],y[1000],e[1000];
  int n=0;
  while((fscanf(fp,"%e %e",x+n,y+n)==2)&&x[n]<xmax)++n;
  //for(int i=0;i<n;++i)printf("%e %e\n",x[i],y[i]);
  for(int i=0;i<n;++i){y[i]*=scale;e[i]=err*y[i];}
  TGraphErrors *g=new TGraphErrors(n,x,y,0,e);
  return(g);
}
TGraphErrors *gReadData2(char *file="star_6m_200_0208_sc1.dat",float scale=1,float err=0){
  FILE *fp=fopen(file,"read");
  float x[1000],y[1000],ex[1000],ey[1000];
  int n=0;
  while(fscanf(fp,"%e %e",x+n,y+n)==2)++n;
  for(int i=0;i<n-1;++i){ex[i]=(x[i+1]-x[i])/2;}ex[n-1]=(x[n-1]-x[n-2])/2;
  for(int i=0;i<n;++i){y[i]*=scale;ey[i]=err*y[i];}
  TGraphErrors *g=new TGraphErrors(n,x,y,ex,ey);
  return(g);
}
TGraphErrors *gReadData3(char *file="file.dat",float xbin_width=0,float magnify=1.0){
  FILE *fp=fopen(file,"read");
  float x[1000],y[1000],y_low[1000],y_high[1000],ex[1000],ey[1000];
  int n=0;
  while(fscanf(fp,"%e %e %e",x+n,y_high+n,y_low+n)==3)++n;
  for(int i=0;i<n;++i)
  {
	  ex[i]=xbin_width/2;
	  ey[i]=magnify*(y_high[i]-y_low[i])/2;
	  y[i]=y_low[i]+ey[i];
	}
  TGraphErrors *g=new TGraphErrors(n,x,y,ex,ey);
  return(g);
}

TGraphErrors *gReadData4(char *file="file.dat",float xbin_width=0,float magnify=1.0){
  FILE *fp=fopen(file,"read");
  float x[1000],y[1000],y_low[1000],y_high[1000],ex[1000],ey[1000],skip;
  int n=0;
  while(fscanf(fp,"%e %e %e %e %e %e %e",x+n,&skip, &skip, &skip, &skip,y_high+n,y_low+n) != EOF)++n;
  for(int i=0;i<n;++i)
  {
	  ex[i]=xbin_width/2;
	  ey[i]=magnify*(y_high[i]-y_low[i])/2;
	  y[i]=y_low[i]+ey[i];
	}
  TGraphErrors *g=new TGraphErrors(n,x,y,ex,ey);
  return(g);
}
//===============================================================
TH1D *getPythiaHisto(char *file="pythia_fullEta1/histos_pythiajet_R0.4.root",float norm=1,float binsize=1,float ptlead=0){
  TFile *fp=new TFile(file,"read");
  TH1D *hpt=(TH1D*)fp->Get(Form("hpT_pTl%g",ptlead));
  hpt->Scale(norm);
  int nrebin=binsize/hpt->GetBinWidth(1);
  printf("rebin full jet pt spectra by %d.\n",nrebin);
  hpt->Rebin(nrebin);hpt->Scale(1./nrebin);
  return(hpt);
}
//===============================================================

TGraph *graph_band(int n,float *x,float *y,float *e,float xmax=10, float xmin=0){
  for(int i=0;i<n;++i)e[i]*=y[i];
  const int nn=2*n+1;float xx[nn],yy[nn];
  int N=0;
  for(int i=0;i<n;++i){
    if(x[i]>xmax || x[i]<xmin)continue;
    xx[N]=x[i];yy[N]=y[i]+e[i];
    N++;
  }
  for(int i=n;i<2*n;++i){
    if(x[2*n-1-i]>xmax || x[2*n-1-i]<xmin)continue;
    xx[N]=x[2*n-1-i];yy[N]=y[2*n-1-i]-e[2*n-1-i];
    N++;
  }
  xx[N]=xx[0];yy[N]=yy[0];
  TGraph *g=new TGraphErrors(N+1,xx,yy);
  return(g);
}

//===============================================================

void plot_h(TH1D *h,int ic=2,int is=3,int iw=1,char *txt="",float xtxt=0.5,float ytxt=0.5,float tsiz=0.04,int sym=0,float xmin=0,float xmax=10){
  h->SetLineColor(ic);
  h->SetLineStyle(is);
  h->SetLineWidth(iw);
  if(sym)h->SetMarkerStyle(abs(sym));
  if(!sym)draw_hist_graph(h,"l",xmin,xmax);
  else if(sym>0)draw_hist_graph(h,"p",xmin,xmax);
  else draw_hist_graph(h,"lp",xmin,xmax);
  if(txt!="")keyLine(xtxt,ytxt,txt,ic,is,tsiz,iw,1);
}

//===============================================================

void plot_g2(TGraphErrors *g1,TGraphErrors *g2,int ic=2,int is=1,int iw=1,char *txt="",float xtxt=0.5,float ytxt=0.5,float tsiz=0.04){
  int n=g1->GetN()+g2->GetN()+1;
  const int N=n;
  float xx[N],yy[N];
  double x,y;
  for(int i=0;i<g1->GetN();++i){
    g1->GetPoint(i,x,y);
    xx[i]=x;yy[i]=y;
  }
  for(int i=0;i<g2->GetN();++i){
    g2->GetPoint(i,x,y);
    xx[N-2-i]=x;yy[N-2-i]=y;
  }
  xx[N-1]=xx[0];yy[N-1]=yy[0];
  TGraph *g=new TGraph(N,xx,yy);
  g->SetLineColor(ic);g->SetLineStyle(is);g->SetLineWidth(iw);g->Draw("l");
  keyLine(xtxt,ytxt,txt,ic,is,tsiz,iw,1);
  return;
  g1->SetLineColor(ic);g1->SetLineStyle(is);g1->SetLineWidth(iw);g1->Draw("l");
  g2->SetLineColor(ic);g2->SetLineStyle(is);g2->SetLineWidth(iw);g2->Draw("l");
  TLine *line=new TLine();
  line->SetLineColor(ic);line->SetLineStyle(is);line->SetLineWidth(iw);
  double x1,y1,x2,y2;
  g1->GetPoint(0,x1,y1);g2->GetPoint(0,x2,y2);line->DrawLine(x1,y1,x2,y2);
  g1->GetPoint(g1->GetN()-1,x1,y1);g2->GetPoint(g2->GetN()-1,x2,y2);line->DrawLine(x1,y1,x2,y2);
  keyLine(xtxt,ytxt,txt,ic,is,tsiz,iw,1);
}


//===============================================================

//make a band from two graphs with same x values 
TGraphErrors* area_graph(TGraphErrors *g1,TGraphErrors *g2)
{
  int n=g1->GetN();
  //const int N=n;
  const int N=100;
  double xx[N],yy[N],yy_err[N];
  double x1,y1,x2,y2;
  for(int i=0;i<g1->GetN();++i){
    g1->GetPoint(i,x1,y1);
	 g2->GetPoint(i,x2,y2);
    xx[i]=x1;
	 yy[i]=y1+y2/2;
	 yy_err[i]=TMath::Abs(y2-y1)/2;
  }
  
  TGraphErrors *g=new TGraphErrors(n,xx,yy,0,yy_err);
  return g;
}

//===============================================================

TH1D *hsyst(int syst=1,TH1D *h0,TH1D *h1=NULL,TH1D *h2=NULL,TH1D *h3=NULL){
  int n0=h0->GetNbinsX();
  int n1=n0,n2=n0,n3=n0;
  if(h1)int n1=h1->GetNbinsX();
  if(h2)int n2=h2->GetNbinsX();
  if(h3)int n3=h3->GetNbinsX();
  if(n1!=n0||n2!=n0||n3!=n0)printf("Nbins ERROR!\n");
  float xmin0=h0->GetBinLowEdge(1);
  float xmin1=xmin0,xmin2=xmin0,xmin3=xmin0;
  if(h1)float xmin1=h1->GetBinLowEdge(1);
  if(h2)float xmin2=h2->GetBinLowEdge(1);
  if(h3)float xmin3=h3->GetBinLowEdge(1);
  if(xmin1!=xmin0||xmin2!=xmin0||xmin3!=xmin0)printf("xmin ERROR!\n");
  float bsize0=h0->GetBinWidth(1);
  float bsize1=bsize0,bsize2=bsize0,bsize3=bsize0;
  if(h1)float bsize1=h1->GetBinWidth(1);
  if(h2)float bsize2=h2->GetBinWidth(1);
  if(h3)float bsize3=h3->GetBinWidth(1);
  if(bsize1!=bsize0||bsize2!=bsize0||bsize3!=bsize0)printf("bin size ERROR!\n");
  TH1D *h=(TH1D*)h0->Clone();
  for(int i=1;i<=n0;++i){
    float ymin=999,ymax=-999;
    if(h1){float y1=h1->GetBinContent(i)-1;if(y1<ymin)ymin=y1;if(y1>ymax)ymax=y1;}
    if(h2){float y2=h2->GetBinContent(i)-1;if(y2<ymin)ymin=y2;if(y2>ymax)ymax=y2;}
    if(h3){float y3=h3->GetBinContent(i)-1;if(y3<ymin)ymin=y3;if(y3>ymax)ymax=y3;}
    float y0=h0->GetBinContent(i)-1;if(y0<ymin)ymin=y0;if(y0>ymax)ymax=y0;//fqw
    if(ymin<0&&ymax<0)ymax=-ymin;else if(ymin>0&&ymax>0)ymin=-ymax;
    //float y0=h0->GetBinContent(i);if(y0>0)ymax=sqrt(ymax*ymax+y0*y0);else ymin=-sqrt(ymin*ymin+y0*y0);//fqw
    if(syst<0)h->SetBinContent(i,ymin+1);else h->SetBinContent(i,ymax+1);
  }
  return(h);
} 

//===============================================================

void ShrinkGraph(TGraphAsymmErrors* graph, float min, float max)
{
	for(int point=(graph->GetMaxSize()+1); point>=0; point--)
	{
		double xpoint,ypoint;
		graph->GetPoint(point, xpoint, ypoint);
		//cout<<point<<": x="<<xpoint<<endl;
		if(xpoint>max || xpoint<min) graph->RemovePoint(point);
	}
	return;
}

//===============================================================

TGraphAsymmErrors* gcombine_syserr(TGraphAsymmErrors *gr1,TGraphAsymmErrors *gr2,bool setErrX=1)
{
	int npoints=gr1->GetN();
	//const int n=npoints;
	const int n=50;
	double gx[n];
	double gx_l[n];
	double gx_h[n];
	double gy[n];
	double gy_l[n];
	double gy_h[n];
	
	for(int i=0;i<npoints;++i){
		double x1,y1;
		gr1->GetPoint(i,x1,y1);
		double x1_down;
		double x1_up;
		
		if(setErrX==0) //set x errors to 0
		x1_down=0; 
		x1_up=0; 
		else
		{
			x1_down=gr1->GetErrorXlow(i);
			x1_up=gr1->GetErrorXhigh(i);
		}
		double y1_down=gr1->GetErrorYlow(i);
		double y1_up=gr1->GetErrorYhigh(i);

		double y2_down=gr2->GetErrorYlow(i);
		double y2_up=gr2->GetErrorYhigh(i);
	
		gx[i]=x1;
		gx_l[i]=x1_down;
		gx_h[i]=x1_up;
		gy[i]=y1;
		gy_l[i]=TMath::Sqrt(y1_down*y1_down+y2_down*y2_down);
		gy_h[i]=TMath::Sqrt(y1_up*y1_up+y2_up*y2_up);
		
	}
	TGraphAsymmErrors* graph=new TGraphAsymmErrors(n, gx, gy, gx_l, gx_h, gy_l, gy_h);
	return graph;
}		
//===============================================================

TGraphErrors *Gdiv(TGraphErrors *g1,TGraphErrors *g2)
{
  float x[1000],y[1000],ex[1000],ey[1000];
  double x1,y1,e1,x2,y2,e2;
  int N=g1->GetN();
  for(int k=0;k<N;++k){
    g1->GetPoint(k,x1,y1);
    g2->GetPoint(k,x2,y2);
    e1=g1->GetErrorY(k);
    e2=g2->GetErrorY(k);
    //printf("%f\t%f\t%f\t%f\t%f\t%f\n",x1,y1,e1,x2,y2,e2);
    x[k]=x1;
    ex[k]=g1->GetErrorX(k);
    if(y2!=0){
      y[k]=y1/y2;
      ey[k]=sqrt(pow(e1/y1,2)+pow(e2/y2,2))*fabs(y[k]);
    }
    else{y[k]=99999;ey[k]=0;}
  }
  TGraphErrors *g=new TGraphErrors(N,x,y,ex,ey);
  g->SetTitle(g1->GetTitle());
  return(g);
}
//===============================================================

TGraphErrors *draw_hist_graph(TH1D *h,char *option="l",float xmin=0,float xmax=100){
  //const int N=h->GetNbinsX();
  //float x[N],y[N],e[N];
  int N=h->GetNbinsX();
  float x[1000],y[1000],e[1000];
  int n=0;
  for(int i=0;i<N;++i){
    float xtmp=h->GetBinCenter(i+1);
    if(xtmp<xmin)continue;
    if(xtmp>xmax)continue;
    x[n]=xtmp;
    y[n]=h->GetBinContent(i+1);
    e[n]=h->GetBinError(i+1);
    n++;
  }
  if(strstr(option,"p"))TGraphErrors *g=new TGraphErrors(n,x,y,0,e);
  else TGraphErrors *g=new TGraphErrors(n,x,y,0,0);
  g->SetMarkerStyle(h->GetMarkerStyle());
  g->SetMarkerColor(h->GetMarkerColor());
  g->SetLineStyle(h->GetLineStyle());
  g->SetLineWidth(h->GetLineWidth());
  g->SetLineColor(h->GetLineColor());
  g->DrawClone(option);
  return(g);
}
//===============================================================

