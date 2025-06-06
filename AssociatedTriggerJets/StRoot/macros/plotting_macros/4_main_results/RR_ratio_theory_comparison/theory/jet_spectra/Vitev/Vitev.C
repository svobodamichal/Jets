void Vitev(){
  //Vitev_RAA(0);
  //Vitev_RAA(1);
  Vitev_R2R4();

}

void Vitev_R2R4(){
  //the following numbers are as same as in chrg_vs_full/compare_chrg_ratio.C
  //Jan R=0.2/R=0.4 data pTlead>5GeV by Web Digitizer
  float Xjan_chrg[22],Yjan_chrg[22];
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
  TGraph *gjan=new TGraph(22,Xjan_chrg,Yjan_chrg);

  TGraphErrors *gppR2R4=gReadData("PPXsecRatio_200GeV_R2overR4.txt",1,0);
  TGraphErrors *gAAR2R4dn=gReadData("AuAuXsecRatio_smallCNM_200GeV_DN_R2overR4.txt",1,0);
  TGraphErrors *gAAR2R4up=gReadData("AuAuXsecRatio_smallCNM_200GeV_UP_R2overR4.txt",1,0);

  TGraphErrors *gppR2R4_dEdx=gReadData("RRpp_02to04_200GeVNLO.dat",1,0);
  TGraphErrors *gAAR2R4up_dEdx=gReadData("RR_NLO-Au0.2to0.4coldUP.dat",1,0);
  TGraphErrors *gAAR2R4dn_dEdx=gReadData("RR_NLO-Au0.2to0.4coldDN.dat",1,0);

  gROOT->LoadMacro("~/SkyDrive/Research/Resources/rootMacros/Utility.C");
  gStyle->SetOptStat(0);
  TLatex *ltx=new TLatex();ltx->SetNDC(1);ltx->SetTextSize(0.038);
  TCanvas *cfig=new TCanvas("cfig","",10,10,600,500);
  cfig->SetMargin(0.15,0.05,0.15,0.05);
  TH2D *h2=new TH2D("h2","",1,0,45,1,0,1.8);
  h2->SetXTitle("p_{T, jet}  [GeV]");h2->SetYTitle("R=0.2/R=0.4");
  h2->SetTitleOffset(1.0,"x");h2->SetTitleOffset(1.1,"y");
  h2->SetTitleSize(0.06,"x");h2->SetTitleSize(0.06,"y");
  h2->SetLabelSize(0.05,"x");h2->SetLabelSize(0.05,"y");
  //h2->SetNdivisions(505,"x");h2->SetNdivisions(505,"y");
  h2->Draw();

  gjan->SetLineColor(1);gjan->SetFillColor(1);gjan->SetFillStyle(3003);gjan->Draw("f");
  ltx->DrawLatex(0.2,0.88,"Hatched: Jan Au+Au data p_{T}^{lead}>5 GeV");
  //ltx->SetTextColor(1);ltx->DrawLatex(0.2,0.88,"Anti-k_{T}");

  ltx->SetTextColor(1);
  ltx->DrawLatex(0.2,0.82,"NLO; Traditional E-loss [PRL104(2010)132001]");
  gppR2R4_dEdx->SetLineColor(1);gppR2R4_dEdx->SetLineStyle(2);gppR2R4_dEdx->Draw("l");
  keyLine(0.2,0.78,"Vitev pp",1,2,0.038);
  gAAR2R4up_dEdx->SetLineColor(1);gAAR2R4up_dEdx->Draw("l");
  gAAR2R4dn_dEdx->SetLineColor(1);gAAR2R4dn_dEdx->Draw("l");
  keyLine(0.5,0.78,"Vitev AuAu w/ syst.",1,1,0.038);

  ltx->SetTextColor(2);
  ltx->DrawLatex(0.2,0.70,"SCET_{G}&med-ind split func[1509.07257], CNM[1509.02936]");
  gppR2R4->SetLineColor(2);gppR2R4->SetLineStyle(2);gppR2R4->Draw("l");
  keyLine(0.2,0.64,"Vitev pp",2,2,0.038);
  gAAR2R4up->SetLineColor(2);gAAR2R4up->Draw("l");
  gAAR2R4dn->SetLineColor(2);gAAR2R4dn->Draw("l");
  keyLine(0.5,0.64,"Vitev AuAu w/ syst.",2,1,0.038);

  cfig->SaveAs("Vitev_R2R4.gif");
  cfig->SaveAs("Vitev_R2R4.pdf");
}

void Vitev_RAA(int include_data=0){
  //Jan RAA data pTlead>5GeV by Web Digitizer
  float XjanR02[100],YjanR02[100];
  int i=0;
  XjanR02[i]=5.156638488314271;YjanR02[i]=0.2768171154707541;++i;
  XjanR02[i]=6.495524614619593;YjanR02[i]=0.20675332909346605;++i;
  XjanR02[i]=8.191198408751863;YjanR02[i]=0.18055498403712578;++i;
  XjanR02[i]=10.459970164097463;YjanR02[i]=0.19550086453242063;++i;
  XjanR02[i]=12.727498756837392;YjanR02[i]=0.21585499772801198;++i;
  XjanR02[i]=15.502237692690201;YjanR02[i]=0.26284014352915547;++i;
  XjanR02[i]=18.422426653406266;YjanR02[i]=0.3263900770369859;++i;
  XjanR02[i]=21.486822476379913;YjanR02[i]=0.4214752596281636;++i;
  XjanR02[i]=26.537792143212332;YjanR02[i]=0.4943055931633678;++i;
  XjanR02[i]=31.741670810542015;YjanR02[i]=0.5258831742065386;++i;
  XjanR02[i]=31.890850323222278;YjanR02[i]=0.1599492706945583;++i;
  XjanR02[i]=27.378170064644454;YjanR02[i]=0.2920245478582252;++i;
  XjanR02[i]=26.495524614619594;YjanR02[i]=0.3034723850218141;++i;
  XjanR02[i]=21.578816509199402;YjanR02[i]=0.3145425691980937;++i;
  XjanR02[i]=18.363998010939834;YjanR02[i]=0.2582411971088554;++i;
  XjanR02[i]=15.518398806563898;YjanR02[i]=0.20395143691595036;++i;
  XjanR02[i]=12.74241670810542;YjanR02[i]=0.17079352706534928;++i;
  XjanR02[i]=10.47613127797116;YjanR02[i]=0.15169936260240582;++i;
  XjanR02[i]=8.063152660367976;YjanR02[i]=0.13472700826987977;++i;
  XjanR02[i]=7.105917454002982;YjanR02[i]=0.14276047667603506;++i;
  XjanR02[i]=6.210840377921429;YjanR02[i]=0.18032255737777533;++i;
  XjanR02[i]=5.171556439582296;YjanR02[i]=0.21902931134760897;++i;
  TGraph *gjanR02=new TGraph(i,XjanR02,YjanR02);
  float XjanR04[100],YjanR04[100];
  int i=0;
  XjanR04[i]=5.43583535108959;YjanR04[i]=0.09492278690306391;++i;
  XjanR04[i]=6.525423728813559;YjanR04[i]=0.10252251705747872;++i;
  XjanR04[i]=8.631961259079905;YjanR04[i]=0.12194838998440816;++i;
  XjanR04[i]=10.593220338983052;YjanR04[i]=0.1508173775316881;++i;
  XjanR04[i]=12.5544794188862;YjanR04[i]=0.17941099782799716;++i;
  XjanR04[i]=15.605326876513317;YjanR04[i]=0.22172664258211916;++i;
  XjanR04[i]=17.929782082324458;YjanR04[i]=0.253650519226906;++i;
  XjanR04[i]=22.36077481840194;YjanR04[i]=0.3385100934426728;++i;
  XjanR04[i]=26.86440677966102;YjanR04[i]=0.4696395831076281;++i;
  XjanR04[i]=32.021791767554475;YjanR04[i]=0.6903833738878398;++i;
  XjanR04[i]=32.021791767554475;YjanR04[i]=0.27699112322642705;++i;
  XjanR04[i]=26.791767554479417;YjanR04[i]=0.25219875176725975;++i;
  XjanR04[i]=22.36077481840194;YjanR04[i]=0.22075859137851356;++i;
  XjanR04[i]=18.002421307506054;YjanR04[i]=0.19322880265504894;++i;
  XjanR04[i]=15.677966101694917;YjanR04[i]=0.16565905112372373;++i;
  XjanR04[i]=12.627118644067796;YjanR04[i]=0.1340436824155624;++i;
  XjanR04[i]=10.520581113801452;YjanR04[i]=0.11269107955138884;++i;
  XjanR04[i]=8.559322033898304;YjanR04[i]=0.09658966974125825;++i;
  XjanR04[i]=6.525423728813559;YjanR04[i]=0.08607358522758517;++i;
  XjanR04[i]=5.43583535108959;YjanR04[i]=0.08125684054861812;++i;
  TGraph *gjanR04=new TGraph(i,XjanR04,YjanR04);

  TGraphErrors *g02up=gReadData("Raa_NLO-Au-Wmin0.R0.2coldUP.dat",1,0);
  TGraphErrors *g02dn=gReadData("Raa_NLO-Au-Wmin0.R0.2coldDN.dat",1,0);
  TGraphErrors *g04up=gReadData("Raa_NLO-Au-Wmin0.R0.4coldUP.dat",1,0);
  TGraphErrors *g04dn=gReadData("Raa_NLO-Au-Wmin0.R0.4coldDN.dat",1,0);

  gStyle->SetOptStat(0);
  TLatex *ltx=new TLatex();ltx->SetNDC(1);
  TCanvas *cfig=new TCanvas("cfig","",10,10,600,500);
  cfig->SetMargin(0.15,0.05,0.15,0.05);
  TH2D *h2=new TH2D("h2","",1,0,45,1,0,0.7);
  h2->SetXTitle("p_{T, jet}  [GeV]");h2->SetYTitle("R_{AA}");
  h2->SetTitleOffset(1.0,"x");h2->SetTitleOffset(1.1,"y");
  h2->SetTitleSize(0.06,"x");h2->SetTitleSize(0.06,"y");
  h2->SetLabelSize(0.05,"x");h2->SetLabelSize(0.05,"y");
  //h2->SetNdivisions(505,"x");h2->SetNdivisions(505,"y");
  h2->Draw();

  if(include_data){
    gjanR02->SetLineColor(1);gjanR02->SetFillColor(1);gjanR02->SetFillStyle(3004);gjanR02->Draw("f");
    gjanR04->SetLineColor(2);gjanR04->SetFillColor(2);gjanR04->SetFillStyle(3005);gjanR04->Draw("f");
    ltx->SetTextColor(1);ltx->DrawLatex(0.2,0.88,"Hatched: data");
    ltx->SetTextColor(1);ltx->DrawLatex(0.2,0.82,"Curves: Vitev");
  }
  g02up->SetLineColor(1);g02up->Draw("l");
  g02dn->SetLineColor(1);g02dn->Draw("l");
  g04up->SetLineColor(2);g04up->Draw("l");
  g04dn->SetLineColor(2);g04dn->Draw("l");
  ltx->SetTextColor(1);ltx->DrawLatex(0.8,0.88,"R=0.2");
  ltx->SetTextColor(2);ltx->DrawLatex(0.8,0.82,"R=0.4");

  if(include_data){
    cfig->SaveAs("Vitev_RAA_vs_JanData.gif");
    cfig->SaveAs("Vitev_RAA_vs_JanData.pdf");
  }else{
    cfig->SaveAs("Vitev_RAA.gif");
    cfig->SaveAs("Vitev_RAA.pdf");
  }
}

TGraphErrors *gReadData(char *file="star_6m_200_0208_sc1.dat",float scale=1,float err=0){
  FILE *fp=fopen(file,"read");
  float x[1000],y[1000],e[1000];
  int n=0;
  while(fscanf(fp,"%f %f",x+n,y+n)==2)++n;
  //for(int i=0;i<n;++i)printf("%g %f\n",x[i],y[i]);
  for(int i=0;i<n;++i){y[i]*=scale;e[i]=err*y[i];}
  TGraphErrors *g=new TGraphErrors(n,x,y,0,e);
  return(g);
}
