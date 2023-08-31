void R2R4(){
  float ptlead=0;
  int Jan_binning=0;
  const int ngrp=13;//used for Jan's binning
  double xbins[ngrp+1]={5,6,7,8,9,10,12,14,16,18,20,25,30,40};

  float ptmax=40;

  gStyle->SetLineStyleString(12,"[40 20]");
  gStyle->SetLineStyleString(13,"[13 13]");
  gStyle->SetLineStyleString(14,"[40 13 13 13]");

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
  TGraph *gjan=new TGraph(i,Xjan_chrg,Yjan_chrg);

  //---------------------------------------------------------------------------
  //Soyez paper
  float Xsoyez_full[10]={12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5};
  float Ysoyez_full[10]={0.297,0.353,0.402,0.426,0.446,0.446,0.452,0.446,0.437,0.399};
  float Esoyez_full[10]={0.20,0.16,0.13,0.13,0.12,0.12,0.11,0.11,0.12,0.13};
  TGraph *gsoyez=graph_band(10,Xsoyez_full,Ysoyez_full,Esoyez_full,ptmax-2);

  //---------------------------------------------------------------------------
  //Vitev
  gROOT->LoadMacro("../jet_spectra/util.C");
  TGraphErrors *gppR2R4=gReadData("../jet_spectra/Vitev/PPXsecRatio_200GeV_R2overR4.txt",1,0,ptmax-2);
  TGraphErrors *gAAR2R4dn=gReadData("../jet_spectra/Vitev/AuAuXsecRatio_smallCNM_200GeV_DN_R2overR4.txt",1,0,ptmax-2);
  TGraphErrors *gAAR2R4up=gReadData("../jet_spectra/Vitev/AuAuXsecRatio_smallCNM_200GeV_UP_R2overR4.txt",1,0,ptmax-2);

  TGraphErrors *gppR2R4_dEdx=gReadData("../jet_spectra/Vitev/RRpp_02to04_200GeVNLO.dat",1,0,ptmax-2);
  TGraphErrors *gAAR2R4up_dEdx=gReadData("../jet_spectra/Vitev/RR_NLO-Au0.2to0.4coldUP.dat",1,0,ptmax-2);
  TGraphErrors *gAAR2R4dn_dEdx=gReadData("../jet_spectra/Vitev/RR_NLO-Au0.2to0.4coldDN.dat",1,0,ptmax-2);

  //===========================================================================
  //Herwig, Pythia
  float binsize=1;
  float scale=1000;//mb->mub
  TH1D *h_ptPy_R02=getPythiaHisto("../jet_spectra/pythia_chrgEta1/histos_pythiajet_R0.2.root",scale*1e-6,binsize,ptlead);//1e6 is nevents
  TH1D *h_ptPy_R04=getPythiaHisto("../jet_spectra/pythia_chrgEta1/histos_pythiajet_R0.4.root",scale*1e-6,binsize,ptlead);
  TH1D *h_ptPy_fullR02=getPythiaHisto("../jet_spectra/pythia_fullEta1/histos_pythiajet_R0.2.root",scale*1e-6,binsize,ptlead);
  TH1D *h_ptPy_fullR04=getPythiaHisto("../jet_spectra/pythia_fullEta1/histos_pythiajet_R0.4.root",scale*1e-6,binsize,ptlead);
  TH1D *h_ptHw_R02=getPythiaHisto("../jet_spectra/herwig_chrgEta1/starjet_herwig_R0.2_chrg.root",scale,binsize,ptlead);//Herwig nevents has already been divided
  TH1D *h_ptHw_R04=getPythiaHisto("../jet_spectra/herwig_chrgEta1/starjet_herwig_R0.4_chrg.root",scale,binsize,ptlead);
  TH1D *h_ptHw_fullR02=getPythiaHisto("../jet_spectra/herwig_fullEta1/starjet_herwig_R0.2_full.root",scale,binsize,ptlead);
  TH1D *h_ptHw_fullR04=getPythiaHisto("../jet_spectra/herwig_fullEta1/starjet_herwig_R0.4_full.root",scale,binsize,ptlead);
  TH1D *h_ptHw_partonR02=getPythiaHisto("../jet_spectra/herwig_partonEta1/starjet_herwig_R0.2_parton.root",scale,binsize,ptlead);
  TH1D *h_ptHw_partonR04=getPythiaHisto("../jet_spectra/herwig_partonEta1/starjet_herwig_R0.4_parton.root",scale,binsize,ptlead);

  //---------------------------------------------------------------------------
  //NLO theory curve
  gROOT->LoadMacro("~/SkyDrive/Research/Resources/rootMacros/UtilityGraph.C");
  TGraphErrors *gNLO2=gReadData("../jet_spectra/Vogelsang/star_ct14_200_05_R02_sc1.dat",scale*1e-9,0);
  TGraphErrors *gNLO4=gReadData("../jet_spectra/Vogelsang/star_ct14_200_05_R04_sc1.dat",scale*1e-9,0);
  TGraphErrors *gNLO=Gdiv(gNLO2,gNLO4);

  //NLO theory curve w/ Herwig parton->full convolution
  TFile *fp=new TFile("herwig_fullEta1partonEta1/NLO_convolute_starjet_herwig_R0.2_parton2full_new2.root","read");
  TH1D *h_NLO_fullR02=(TH1D*)fp->Get("hNLO_conv");h_NLO_fullR02->Scale(scale);
  TFile *fp=new TFile("herwig_fullEta1partonEta1/NLO_convolute_starjet_herwig_R0.4_parton2full_new2.root","read");
  TH1D *h_NLO_fullR04=(TH1D*)fp->Get("hNLO_conv");h_NLO_fullR04->Scale(scale);

  //NLO theory curve w/ Herwig parton->chrg convolution
  TFile *fp=new TFile("herwig_chrgEta1partonEta1/NLO_convolute_starjet_herwig_R0.2_parton2chrg_new2.root","read");
  TH1D *h_NLO_R02=(TH1D*)fp->Get("hNLO_conv");h_NLO_R02->Scale(scale);
  TH1D *h_NLOsc05_R02=(TH1D*)fp->Get("hNLOsc05_conv");h_NLOsc05_R02->Scale(scale);
  TH1D *h_NLOsc2_R02=(TH1D*)fp->Get("hNLOsc2_conv");h_NLOsc2_R02->Scale(scale);
  TFile *fp=new TFile("herwig_chrgEta1partonEta1/NLO_convolute_starjet_herwig_R0.4_parton2chrg_new2.root","read");
  TH1D *h_NLO_R04=(TH1D*)fp->Get("hNLO_conv");h_NLO_R04->Scale(scale);
  TH1D *h_NLOsc05_R04=(TH1D*)fp->Get("hNLOsc05_conv");h_NLOsc05_R04->Scale(scale);
  TH1D *h_NLOsc2_R04=(TH1D*)fp->Get("hNLOsc2_conv");h_NLOsc2_R04->Scale(scale);

  /*
  //NLO Herwig/Pythia convolution
  TFile *fp=new TFile(Form("NLO_Herwig_Pythia_convolute_R%g_new2.root",R),"read");
  TH1D *h_NLO_HwPy=(TH1D*)fp->Get("hNLO_chrg_HWPYconv");h_NLO_HwPy->Scale(scale);
  TH1D *h_NLOsc05_HwPy=(TH1D*)fp->Get("hNLOsc05_chrg_HWPYconv");h_NLOsc05_HwPy->Scale(scale);
  TH1D *h_NLOsc2_HwPy=(TH1D*)fp->Get("hNLOsc2_chrg_HWPYconv");h_NLOsc2_HwPy->Scale(scale);
  */
  //---------------------------------------------------------------------------
  if(Jan_binning){
    h_ptPy_R02->Rebin(ngrp,"hptPy_R02",xbins);
    h_ptHw_R02->Rebin(ngrp,"hptHw_R02",xbins);
    h_NLO_R02->Rebin(ngrp,"hNLO_R02",xbins);
    h_NLOsc05_R02->Rebin(ngrp,"hNLOsc05_R02",xbins);
    h_NLOsc2_R02->Rebin(ngrp,"hNLOsc2_R02",xbins);
    h_ptPy_R04->Rebin(ngrp,"hptPy_R04",xbins);
    h_ptHw_R04->Rebin(ngrp,"hptHw_R04",xbins);
    h_NLO_R04->Rebin(ngrp,"hNLO_R04",xbins);
    h_NLOsc05_R04->Rebin(ngrp,"hNLOsc05_R04",xbins);
    h_NLOsc2_R04->Rebin(ngrp,"hNLOsc2_R04",xbins);

    h_NLO_fullR02->Rebin(ngrp,"hNLO_fullR02",xbins);
    h_NLO_fullR04->Rebin(ngrp,"hNLO_fullR04",xbins);
    h_ptPy_fullR02->Rebin(ngrp,"hptPy_fullR02",xbins);
    h_ptHw_fullR02->Rebin(ngrp,"hptHw_fullR02",xbins);
    h_ptPy_fullR04->Rebin(ngrp,"hptPy_fullR04",xbins);
    h_ptHw_fullR04->Rebin(ngrp,"hptHw_fullR04",xbins);
    h_ptHw_partonR02->Rebin(ngrp,"hptHw_partonR02",xbins);
    h_ptHw_partonR04->Rebin(ngrp,"hptHw_partonR04",xbins);
    /*h_NLO_HwPy->Rebin(ngrp,"hNLO_HwPy",xbins);
    h_NLOsc05_HwPy->Rebin(ngrp,"hNLOsc05_HwPy",xbins);
    h_NLOsc2_HwPy->Rebin(ngrp,"hNLOsc2_HwPy",xbins);*/
    //take care of bin size
    hptPy_R02->Scale(h_ptPy_R02->GetBinWidth(1),"width");
    hptHw_R02->Scale(h_ptHw_R02->GetBinWidth(1),"width");
    hNLO_R02->Scale(h_NLO_R02->GetBinWidth(1),"width");
    hNLOsc05_R02->Scale(h_NLOsc05_R02->GetBinWidth(1),"width");
    hNLOsc2_R02->Scale(h_NLOsc2_R02->GetBinWidth(1),"width");
    hptPy_R04->Scale(h_ptPy_R04->GetBinWidth(1),"width");
    hptHw_R04->Scale(h_ptHw_R04->GetBinWidth(1),"width");
    hNLO_R04->Scale(h_NLO_R04->GetBinWidth(1),"width");
    hNLOsc05_R04->Scale(h_NLOsc05_R04->GetBinWidth(1),"width");
    hNLOsc2_R04->Scale(h_NLOsc2_R04->GetBinWidth(1),"width");

    hNLO_fullR02->Scale(h_NLO_fullR02->GetBinWidth(1),"width");
    hNLO_fullR04->Scale(h_NLO_fullR04->GetBinWidth(1),"width");
    hptPy_fullR02->Scale(h_ptPy_fullR02->GetBinWidth(1),"width");
    hptHw_fullR02->Scale(h_ptHw_fullR02->GetBinWidth(1),"width");
    hptPy_fullR04->Scale(h_ptPy_fullR04->GetBinWidth(1),"width");
    hptHw_fullR04->Scale(h_ptHw_fullR04->GetBinWidth(1),"width");
    hptHw_partonR02->Scale(h_ptHw_partonR02->GetBinWidth(1),"width");
    hptHw_partonR04->Scale(h_ptHw_partonR04->GetBinWidth(1),"width");
    /*hNLO_HwPy->Scale(h_NLO_HwPy->GetBinWidth(1),"width");
    hNLOsc05_HwPy->Scale(h_NLOsc05_HwPy->GetBinWidth(1),"width");
    hNLOsc2_HwPy->Scale(h_NLOsc2_HwPy->GetBinWidth(1),"width");*/
  }else{
    TH1D *hptPy_R02=(TH1D*)h_ptPy_R02->Clone();
    TH1D *hptHw_R02=(TH1D*)h_ptHw_R02->Clone();
    TH1D *hNLO_R02=(TH1D*)h_NLO_R02->Clone();
    TH1D *hNLOsc05_R02=(TH1D*)h_NLOsc05_R02->Clone();
    TH1D *hNLOsc2_R02=(TH1D*)h_NLOsc2_R02->Clone();
    TH1D *hptPy_R04=(TH1D*)h_ptPy_R04->Clone();
    TH1D *hptHw_R04=(TH1D*)h_ptHw_R04->Clone();
    TH1D *hNLO_R04=(TH1D*)h_NLO_R04->Clone();
    TH1D *hNLOsc05_R04=(TH1D*)h_NLOsc05_R04->Clone();
    TH1D *hNLOsc2_R04=(TH1D*)h_NLOsc2_R04->Clone();

    TH1D *hNLO_fullR02=(TH1D*)h_NLO_fullR02->Clone();
    TH1D *hNLO_fullR04=(TH1D*)h_NLO_fullR04->Clone();
    TH1D *hptPy_fullR02=(TH1D*)h_ptPy_fullR02->Clone();
    TH1D *hptHw_fullR02=(TH1D*)h_ptHw_fullR02->Clone();
    TH1D *hptPy_fullR04=(TH1D*)h_ptPy_fullR04->Clone();
    TH1D *hptHw_fullR04=(TH1D*)h_ptHw_fullR04->Clone();
    TH1D *hptHw_partonR02=(TH1D*)h_ptHw_partonR02->Clone();
    TH1D *hptHw_partonR04=(TH1D*)h_ptHw_partonR04->Clone();
    /*TH1D *hNLO_HwPy=(TH1D*)h_NLO_HwPy->Clone();
    TH1D *hNLOsc05_HwPy=(TH1D*)h_NLOsc05_HwPy->Clone();
    TH1D *hNLOsc2_HwPy=(TH1D*)h_NLOsc2_HwPy->Clone();*/
  }
  hptPy_R02->SetName("hptPy_R02");
  hptHw_R02->SetName("hptHw_R02");
  hNLO_R02->SetName("hNLO_R02");
  hNLOsc05_R02->SetName("hNLOsc05_R02");
  hNLOsc2_R02->SetName("hNLOsc2_R02");
  hptPy_R04->SetName("hptPy_R04");
  hptHw_R04->SetName("hptHw_R04");
  hNLO_R04->SetName("hNLO_R04");
  hNLOsc05_R04->SetName("hNLOsc05_R04");
  hNLOsc2_R04->SetName("hNLOsc2_R04");

  hNLO_fullR02->SetName("hNLO_fullR02");
  hNLO_fullR04->SetName("hNLO_fullR04");
  hptPy_fullR02->SetName("hptPy_fullR02");
  hptHw_fullR02->SetName("hptHw_fullR02");
  hptPy_fullR04->SetName("hptPy_fullR04");
  hptHw_fullR04->SetName("hptHw_fullR04");
  hptHw_partonR02->SetName("hptHw_partonR02");
  hptHw_partonR04->SetName("hptHw_partonR04");
  /*hNLO_HwPy->SetName("hNLO_HwPy");
  hNLOsc05_HwPy->SetName("hNLOsc05_HwPy");
  hNLOsc2_HwPy->SetName("hNLOsc2_HwPy");*/

  //ratio R=0.2/R=0.4
  TH1D *hptPy=(TH1D*)hptPy_R02->Clone();hptPy->Divide(hptPy_R04);
  TH1D *hptHw=(TH1D*)hptHw_R02->Clone();hptHw->Divide(hptHw_R04);
  TH1D *hNLO=(TH1D*)hNLO_R02->Clone();hNLO->Divide(hNLO_R04);
  TH1D *hNLOsc05=(TH1D*)hNLOsc05_R02->Clone();hNLOsc05->Divide(hNLOsc05_R04);
  TH1D *hNLOsc2=(TH1D*)hNLOsc2_R02->Clone();hNLOsc2->Divide(hNLOsc2_R04);

  TH1D *hNLO_full=(TH1D*)hNLO_fullR02->Clone();hNLO_full->Divide(hNLO_fullR04);
  TH1D *hptPy_full=(TH1D*)hptPy_fullR02->Clone();hptPy_full->Divide(hptPy_fullR04);
  TH1D *hptHw_full=(TH1D*)hptHw_fullR02->Clone();hptHw_full->Divide(hptHw_fullR04);
  TH1D *hptHw_parton=(TH1D*)hptHw_partonR02->Clone();hptHw_parton->Divide(hptHw_partonR04);

  //===========================================================================
  gROOT->LoadMacro("~/SkyDrive/Research/Resources/rootMacros/Utility.C");
  gStyle->SetOptStat(0);

  TCanvas *cfig1=new TCanvas("cfig1","fig1",10,10,600,500);
  cfig1->SetMargin(0.15,0.03,0.17,0.03);
  cfig1->SetTicks(1);//cfig1->SetLogy(1);
  TH1D *htmp1=new TH1D("htmp1","",100,0,ptmax);
  htmp1->SetXTitle("p_{T, jet}  or  p_{T, jet}^{ch}  [GeV/c]");htmp1->SetYTitle("Jet spectra ratio R=0.2/R=0.4");
  htmp1->SetTitleOffset(1.2,"x");htmp1->SetTitleOffset(1.1,"y");
  htmp1->SetTitleSize(0.06,"x");htmp1->SetTitleSize(0.06,"y");
  htmp1->SetLabelSize(0.05,"x");htmp1->SetLabelSize(0.05,"y");
  //htmp1->SetNdivisions(505,"x");htmp1->SetNdivisions(505,"y");
  htmp1->SetMinimum(0.1);htmp1->SetMaximum(1.5);
  htmp1->Draw();
  float ytxt=0.90,dytxt=0.055,tsize=0.038;
  TLine *line=new TLine();line->SetLineStyle(1);line->DrawLine(0,1,ptmax,1);
  TLatex *ltx=new TLatex();ltx->SetNDC(1);ltx->SetTextSize(tsize);
  ltx->SetTextColor(1);ltx->DrawLatex(0.25,0.90,"200 GeV, Anti-k_{T}");

  //gjan->SetFillStyle(3003);gjan->SetFillColor(1);gjan->Draw("f");
  gjan->SetLineWidth(2);gjan->SetLineColor(2);gjan->Draw("l");
  keyLine(0.08,0.86,"STAR AuAu 0-10%",2,1,tsize,2);
  ltx->SetTextColor(2);ltx->SetTextSize(0.035);ltx->DrawLatex(0.3,0.80,"p_{T}^{lead}>5 GeV/c");

  //plot_h(hNLOsc05,2,3,1," ",0.2,ytxt+0.01,tsize,0,5,ptmax-2);
  //plot_h(hNLO    ,2,1,2,"NLO Herwig parton->chrg",0.2,ytxt,tsize,0,5,ptmax-2);
  //plot_h(hNLOsc2 ,2,3,1," ",0.2,ytxt-0.01,tsize,0,5,ptmax-2);ytxt-=dytxt;
  plot_h(hptHw,4,1,1,"Herwig pp (ch.)",0.55,ytxt,tsize,0,5,ptmax-2);ytxt-=dytxt;
  plot_h(hptPy,4,2,1,"Pythia pp (ch.)",0.55,ytxt,tsize,0,5,ptmax-2);ytxt-=dytxt;

  //full jets
  //plot_h(hNLO_full,2,2,2,"NLO Herwig parton->full",0.2,ytxt,tsize,0,5,ptmax-2);ytxt-=dytxt;
  //plot_h(hptHw_full,1,2,2,"Herwig full",0.2,ytxt,tsize,0,5,ptmax-2);ytxt-=dytxt;
  //plot_h(hptPy_full,4,2,2,"Pythia full",0.2,ytxt,tsize,0,5,ptmax-2);ytxt-=dytxt;
  //gNLO->SetLineColor(2);gNLO->Draw("l");ytxt+=0.04;keyLine(0.2,ytxt,"NLO",2,1);ytxt-=dytxt;
  //plot_h(hptHw_parton,1,1,1,"",0.2,ytxt,tsize);ytxt-=dytxt;keyLine(0.2,ytxt,"Herwig parton",1,1);

  gsoyez->SetLineColor(4);gsoyez->SetMarkerColor(4);gsoyez->SetLineStyle(3);gsoyez->Draw("l");keyLine(0.55,ytxt,"Soyez pp NLO (hadr.)",4,3,tsize);ytxt-=dytxt;

  //float ytxt=0.82,dytxt=0.055,tsize=0.038;
  //ltx->SetTextColor(1);ltx->DrawLatex(xtxt1,ytxt,"NLO/dEdx [PRL104(2010)132001]");ytxt-=dytxt;
  gppR2R4_dEdx->SetLineColor(4);gppR2R4_dEdx->SetLineStyle(12);gppR2R4_dEdx->Draw("l");keyLine(0.55,ytxt,"Vitev pp (Model 1)",4,12,tsize);
  plot_g2(gAAR2R4up_dEdx,gAAR2R4dn_dEdx,2,12,1,"Vitev AuAu (Model 1)",0.08,ytxt,tsize);ytxt-=dytxt;

  //ltx->SetTextColor(2);ltx->DrawLatex(xtxt1,ytxt,"SCET_{G}&med-ind split func[1509.07257], CNM[1509.02936]");ytxt-=dytxt;
  gppR2R4->SetLineColor(4);gppR2R4->SetLineStyle(14);gppR2R4->Draw("l");keyLine(0.55,ytxt,"Vitev pp (Model 2)",4,14,tsize);
  plot_g2(gAAR2R4up,gAAR2R4dn,2,14,1,"Vitev AuAu (Model 2)",0.08,ytxt,tsize);ytxt-=dytxt;

  cfig1->SaveAs("plots/R2R4.gif");
  cfig1->SaveAs("plots/R2R4.pdf"); 
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

TGraph *graph_band(int n,float *x,float *y,float *e,float xmax=10){
  for(int i=0;i<n;++i)e[i]*=y[i];
  const int nn=2*n+1;float xx[nn],yy[nn];
  int N=0;
  for(int i=0;i<n;++i){
    if(x[i]>xmax)continue;
    xx[N]=x[i];yy[N]=y[i]+e[i];
    N++;
  }
  for(int i=n;i<2*n;++i){
    if(x[2*n-1-i]>xmax)continue;
    xx[N]=x[2*n-1-i];yy[N]=y[2*n-1-i]-e[2*n-1-i];
    N++;
  }
  xx[N]=xx[0];yy[N]=yy[0];
  TGraph *g=new TGraphErrors(N+1,xx,yy);
  return(g);
}

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
