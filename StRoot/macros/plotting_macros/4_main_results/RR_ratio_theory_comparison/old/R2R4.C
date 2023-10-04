#include "../../utils/common.h"

void R2R4(C_systems system=cent, float pTlead=5.0,int binning=1, TString evo="omicron",TString label="THIS THESIS", TString ext="pdf"){
	
	
	//what do we want to show
	bool showVitevAuAu=1;
	bool showVitevpp=0;
	bool showSoyez=0;
	bool showKrishna=1;
	bool showHerwig=1;
	bool showPythia=1;
	bool biasedPP=0; //use pTlead cut also for PYTHIA and HERWIG?
	
	if(system==peri)
	{
		showVitevAuAu=0;
		showKrishna=0;
	}
  float ptlead=pTlead;
  int Jan_binning=0;
  const int ngrp=13;//used for Jan's binning
  double xbins[ngrp+1]={5,6,7,8,9,10,12,14,16,18,20,25,30,40};

  float ptmax=(system==cent) ? 40:40; //40:30

  gStyle->SetLineStyleString(12,"[40 20]");
  gStyle->SetLineStyleString(13,"[13 13]");
  gStyle->SetLineStyleString(14,"[40 13 13 13]");

  TString outDir=Form("../../../plotting_out/obr/%s/results/comparison",evo.Data());
  
  TString cent_dir=C_system_name[system];
  TString cent_leg= (system==cent) ? "central (0-10%)" : "peripheral (60-80%)";
  
  //clors
  Color_t cPy=kRed;
  Color_t cHw=kRed;
  Color_t cNLO=kRed-2;
  Color_t cSTAR=kBlue;
  Color_t cVit1=kCyan-8;
  Color_t cVit2=kAzure-8;
  Color_t ckrsh=kAzure+3;
  
  //line width
  int lwPy=2;
  int lwHw=2;
  int lwNLO=2;
  int lwSTAR=2;
  int lwSTAR_sys=1;
  int lwVit1=1;
  int lwVit2=2;
  
  //line style
  int lsPy=14;
  int lsHw=2;
  int lsNLO=3;
  int lsSTAR=1;
  int lsVit1=1;
  int lsVit2=4;
  
  //marker style and size
  int msSTAR=29;
  float msizeSTAR=1.5;
  
  //fill style
  int vit1Fill=1001; //Ivans's RRratio fill style
  int vit2Fill=1001; //Ivans's RRratio fill style
  int krshFill=1001; //Krishna's RRratio fill style
  
  
  //label
  	float G_latex_sz=0.032;
	Color_t G_latex_cl=kGray+1;
  
  //===========================================================================
  //Jan R=0.2/R=0.4 data pTlead>5GeV by Web Digitizer//need number from him
  float Xjan_chrg[100],Yjan_chrg[100];
  int i=0;
  Xjan_chrg[i]=7.4074074074074066;Yjan_chrg[i]=0.8695278969957081;++i;
  Xjan_chrg[i]=8.37037037037037;Yjan_chrg[i]=0.7596566523605148;++i;
  Xjan_chrg[i]=9.555555555555557;Yjan_chrg[i]=0.7081545064377681;++i;
  Xjan_chrg[i]=10.74074074074074;Yjan_chrg[i]=0.6738197424892705;++i;
  Xjan_chrg[i]=13.11111111111111;Yjan_chrg[i]=0.6738197424892705;++i;
  Xjan_chrg[i]=14.74074074074074;Yjan_chrg[i]=0.7012875536480685;++i;
  Xjan_chrg[i]=16.296296296296294;Yjan_chrg[i]=0.7150214592274677;++i;
  Xjan_chrg[i]=19.77777777777778;Yjan_chrg[i]=0.7768240343347639;++i;
  Xjan_chrg[i]=21.481481481481477;Yjan_chrg[i]=0.8145922746781116;++i;
  Xjan_chrg[i]=26.22222222222222;Yjan_chrg[i]=0.7253218884120172;++i;
  Xjan_chrg[i]=31.925925925925924;Yjan_chrg[i]=0.550214592274678;++i;
  Xjan_chrg[i]=32.07407407407407;Yjan_chrg[i]=0.29957081545064357;++i;
  Xjan_chrg[i]=26.148148148148145;Yjan_chrg[i]=0.5090128755364804;++i;
  Xjan_chrg[i]=21.555555555555554;Yjan_chrg[i]=0.6394849785407724;++i;
  Xjan_chrg[i]=19.703703703703702;Yjan_chrg[i]=0.6120171673819743;++i;
  Xjan_chrg[i]=16.444444444444443;Yjan_chrg[i]=0.5193133047210299;++i;
  Xjan_chrg[i]=14.74074074074074;Yjan_chrg[i]=0.49184549356223184;++i;
  Xjan_chrg[i]=13.11111111111111;Yjan_chrg[i]=0.4952789699570814;++i;
  Xjan_chrg[i]=10.666666666666668;Yjan_chrg[i]=0.5158798283261803;++i;
  Xjan_chrg[i]=9.407407407407405;Yjan_chrg[i]=0.5296137339055793;++i;
  Xjan_chrg[i]=8.518518518518517;Yjan_chrg[i]=0.5433476394849786;++i;
  Xjan_chrg[i]=7.4074074074074066;Yjan_chrg[i]=0.6154506437768239;++i;
  Xjan_chrg[i]=Xjan_chrg[0];Yjan_chrg[i]=Yjan_chrg[0];++i;
  //TGraph *gjan=new TGraph(i,Xjan_chrg,Yjan_chrg);
  
  TString fname=Form("../../../plotting_out/systematics/%s/%s/bining%i/systematics_R0.4_pT%0.lf.root",cent_dir.Data(),evo.Data(),binning,pTlead);
  TFile *fp=new TFile(fname,"read");
  cout<<"opening file: "<<fname.Data()<<endl;
  
  	TGraphAsymmErrors *gjan_unf=fp->Get("RR_unfold_BSL");
	TGraphAsymmErrors *gjan_corr=fp->Get("RR_corr_BSL");
	TGraphAsymmErrors *gjan=fp->Get("RR_BSL");

	float xMax= (system==cent) ? 30 : 25;
	ShrinkGraph(gjan_unf,pTlead,xMax);
	ShrinkGraph(gjan_corr,pTlead,xMax);
	ShrinkGraph(gjan,pTlead,xMax);
	TGraphAsymmErrors *gjan_sys=gcombine_syserr(gjan_unf,gjan_corr,1);
  


  //---------------------------------------------------------------------------
  //Soyez paper
  float Xsoyez_full[10]={12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5};
  float Ysoyez_full[10]={0.297,0.353,0.402,0.426,0.446,0.446,0.452,0.446,0.437,0.399};
  float Esoyez_full[10]={0.20,0.16,0.13,0.13,0.12,0.12,0.11,0.11,0.12,0.13};
  TGraph *gsoyez=graph_band(10,Xsoyez_full,Ysoyez_full,Esoyez_full,ptmax-2);

  //---------------------------------------------------------------------------
  //Vitev
  gROOT->LoadMacro("./util.C");
  TGraphErrors *gppR2R4=gReadData("theory/jet_spectra/Vitev/PPXsecRatio_200GeV_R2overR4.txt",1,0,ptmax-2);
  TGraphErrors *gAAR2R4dn=gReadData("theory/jet_spectra/Vitev/AuAuXsecRatio_smallCNM_200GeV_DN_R2overR4.txt",1,0,ptmax-2);
  TGraphErrors *gAAR2R4up=gReadData("theory/jet_spectra/Vitev/AuAuXsecRatio_smallCNM_200GeV_UP_R2overR4.txt",1,0,ptmax-2);

  TGraphErrors *gppR2R4_dEdx=gReadData("theory/jet_spectra/Vitev/RRpp_02to04_200GeVNLO.dat",1,0,ptmax-2);
  TGraphErrors *gAAR2R4up_dEdx=gReadData("theory/jet_spectra/Vitev/RR_NLO-Au0.2to0.4coldUP.dat",1,0,ptmax-2);
  TGraphErrors *gAAR2R4dn_dEdx=gReadData("theory/jet_spectra/Vitev/RR_NLO-Au0.2to0.4coldDN.dat",1,0,ptmax-2);

  //---------------------------------------------------------------------------
  //Krishna
  TGraphErrors *gKrishnaR2R4=gReadData3(Form("theory/hybrid_results/cent010jetOver02_R04.dat"),5);
  
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
  gROOT->LoadMacro("./UtilityGraph.C");
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
  


  //===========================================================================
  //Draw plot
  //===========================================================================

  gROOT->LoadMacro("./Utility.C");
  gStyle->SetOptStat(0);

  	TLatex *latexL = new TLatex();
	latexL->SetNDC();
	latexL->SetTextSize(G_latex_sz);
	latexL->SetTextColor(G_latex_cl);
  
  TCanvas *cfig1=new TCanvas("cfig1","fig1",10,10,800,600);
  cfig1->SetMargin(0.15,0.03,0.17,0.03);
  cfig1->SetTicks(1);//cfig1->SetLogy(1);
  TH1D *htmp1=new TH1D("htmp1","",100,0,ptmax);
  htmp1->SetXTitle("p_{T, jet}  or  p_{T, jet}^{ch}  (GeV/c)");htmp1->SetYTitle("R=0.2/R=0.4");
  htmp1->SetTitleOffset(1.2,"x");htmp1->SetTitleOffset(1.1,"y");
  htmp1->SetTitleSize(0.06,"x");htmp1->SetTitleSize(0.06,"y");
  htmp1->SetLabelSize(0.05,"x");htmp1->SetLabelSize(0.05,"y");
  //htmp1->SetNdivisions(505,"x");htmp1->SetNdivisions(505,"y");
  htmp1->SetMinimum(0.1);htmp1->SetMaximum(2.0);
  htmp1->Draw();
  float ytxt=-1,dytxt=0.055,tsize=0.038;
  TLine *line=new TLine();line->SetLineStyle(1);line->DrawLine(0,1,ptmax,1);
  TLatex *ltx=new TLatex();ltx->SetNDC(1);ltx->SetTextSize(tsize);
  ltx->SetTextColor(1);ltx->DrawLatex(0.25,0.90,"Charged jets, anti-k_{T}");

  //Pythia and Herwig
  if(showHerwig) plot_h(hptHw,cHw,lsHw,lwHw,"Herwig pp (ch.)",0.55,ytxt,tsize,0,pTlead,ptmax-2);ytxt-=dytxt;
  if(showPythia) plot_h(hptPy,cPy,lsPy,lwPy,"Pythia pp (ch.)",0.55,ytxt,tsize,0,pTlead,ptmax-2);ytxt-=dytxt;

  //Soyez NLO
  gsoyez->SetLineColor(cNLO);
  gsoyez->SetLineWidth(lwNLO);
  gsoyez->SetMarkerColor(cNLO);
  gsoyez->SetLineStyle(lsNLO);
  if(showSoyez) gsoyez->Draw("l");
  //keyLine(0.55,ytxt,"Soyez pp NLO (hadr.)",cpp1,3,tsize);ytxt-=dytxt;

  //Vitev AuAu
  /*
  //gppR2R4_dEdx->SetLineColor(cpp1);gppR2R4_dEdx->SetLineStyle(12);gppR2R4_dEdx->Draw("l");keyLine(0.55,ytxt,"Vitev pp (Model 1)",cpp1,12,tsize);
  if(showVitevAuAu) plot_g2(gAAR2R4up_dEdx,gAAR2R4dn_dEdx,cVit1,lsVit1,lwVit1,"Vitev AuAu (Model 1)",0.08,ytxt,tsize);ytxt-=dytxt;
  gAAR2R4up_dEdx->SetLineColor(cVit1);
  gAAR2R4up_dEdx->SetLineStyle(lsVit1);
  gAAR2R4up_dEdx->SetLineWidth(lwVit1);
  gAAR2R4up_dEdx->SetFillStyle(vit1Fill);
  gAAR2R4up_dEdx->SetFillColor(cVit1);

  //gppR2R4->SetLineColor(cpp1);gppR2R4->SetLineStyle(14);gppR2R4->Draw("l");keyLine(0.55,ytxt,"Vitev pp (Model 2)",cpp1,14,tsize);
  if(showVitevAuAu) plot_g2(gAAR2R4up,gAAR2R4dn,cVit2,lsVit2,lwVit2,"Vitev AuAu (Model 2)",0.08,-1,tsize);ytxt-=dytxt;
  gAAR2R4up->SetLineColor(cVit2);
  gAAR2R4up->SetLineStyle(lsVit2);
  gAAR2R4up->SetLineWidth(lwVit2);
  gAAR2R4up->SetFillStyle(vit2Fill);
  gAAR2R4up->SetFillColor(cVit2);
  */
  TGraphErrors* gvit1=area_graph(gAAR2R4up_dEdx,gAAR2R4dn_dEdx);
  gvit1->SetLineColor(cVit1);
  //gvit1->SetLineStyle(lsVit1);
  //gvit1->SetLineWidth(lwVit1);
  gvit1->SetFillStyle(vit1Fill);
  gvit1->SetFillColor(cVit1);
  if(showVitevAuAu) gvit1->Draw("3");
    
  TGraphErrors* gvit2=area_graph(gAAR2R4up,gAAR2R4dn);
  gvit2->SetLineColor(cVit2);
  //gvit2->SetLineStyle(lsVit2);
  //gvit2->SetLineWidth(lwVit2);
  gvit2->SetFillStyle(vit2Fill);
  gvit2->SetFillColor(cVit2);
  if(showVitevAuAu) gvit2->Draw("3");
  
  //Krishna
  gKrishnaR2R4->SetLineColor(ckrsh);
  gKrishnaR2R4->SetFillStyle(krshFill);
  gKrishnaR2R4->SetFillColor(ckrsh);
  if(showKrishna) gKrishnaR2R4->Draw("3");
  
    //STAR data
  gjan->SetLineWidth(lwSTAR);
  gjan->SetLineColor(cSTAR);
  gjan->SetMarkerColor(cSTAR);
  gjan->Draw("p");
  gjan_sys->SetLineColor(cSTAR);
  gjan_sys->SetFillStyle(0);
  gjan_sys->SetLineWidth(lwSTAR_sys);
  gjan_sys->Draw("2");
  
    //legend
	double posx1=0.6;
	double posx2=0.9;
	double posy1=(system==cent) ? 0.55 : 0.65;  //gppR2R4->SetLineColor(cpp1);gppR2R4->SetLineStyle(14);gppR2R4->Draw("l");keyLine(0.55,ytxt,"Vitev pp (Model 2)",cpp1,14,tsize);
	double posy2=(system==cent) ? 0.95 : 0.85;
	TLegend *legraa = new TLegend(posx1,posy1,posx2,posy2);
	legraa->SetTextSize(0.032);
	legraa->SetFillStyle(0);
	legraa->SetBorderSize(0);
	legraa->AddEntry("","Au+Au #sqrt{s_{NN}}= 200 GeV","");
	legraa->AddEntry("",cent_leg,"");
	legraa->AddEntry(gjan, Form("  STAR, p_{T,lead}^{min}= %.1lf GeV/c",pTlead),"lp");
	//legraa->AddEntry("",Form("    p_{T}^{lead}>%.1lf",pTlead),"");
	if(system==cent)
	if(showVitevAuAu) legraa->AddEntry(gvit1/*gAAR2R4up_dEdx*/, "  NLO pQCD (full jets)","f");
	if(showVitevAuAu) legraa->AddEntry(gvit2/*gAAR2R4up*/, "  SCET (full jets)","f");
	if(showKrishna) legraa->AddEntry(gKrishnaR2R4, "  Hybrid model (full jets)", "f");
	legraa->Draw("same");
	
	    //legend
	posx1=0.35;
	posx2=0.6;
	posy1=0.58;
	posy2=0.88;
	TLegend *legraa = new TLegend(posx1,posy1,posx2,posy2);
	legraa->SetTextSize(0.032);
	legraa->SetFillStyle(0);
	legraa->SetBorderSize(0);
	legraa->AddEntry("","p+p #sqrt{s_{NN}}=200GeV","");
	if(showHerwig) 
	{
		legraa->AddEntry(hptHw, "  Herwig","lp");
		if (biasedPP) legraa->AddEntry("", Form("    p_{T,lead}^{min}= %.1lf GeV/c",pTlead),"");
	}
	//legraa->AddEntry("",Form("    p_{T}^{lead}>%.1lf",pTlead),"");
	if(showPythia) 
	{
		legraa->AddEntry(hptPy,"  PYTHIA","lp");
		if (biasedPP)  legraa->AddEntry("", Form("    p_{T,lead}^{min}= %.1lf GeV/c",pTlead),"");
	}
	//legraa->AddEntry("",Form("    p_{T}^{lead}>%.1lf",pTlead),"");
	if(showSoyez) legraa->AddEntry(gsoyez, "  Soyez NLO (full jets)","lp");
	legraa->Draw("same");
  
	latexL->DrawLatex(0.2,0.25,label);
	
  cfig1->SaveAs(Form("%s/R2R4_pTl%.0lf_bin%i_%s.%s",outDir.Data(),pTlead,binning,cent_dir.Data(),ext.Data()));
}

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
  const int N=n;
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
	const int n=npoints;
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
