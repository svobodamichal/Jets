void NLO_all(){
  gROOT->LoadMacro("~/SkyDrive/Research/Resources/rootMacros/Utility.C");
  gROOT->LoadMacro("~/SkyDrive/Research/Resources/rootMacros/UtilityGraph.C");
  gROOT->LoadMacro("../util.C");
  /*
  TGraphErrors *g2006=gReadData("NLO_PRL2006_curve.dat");
  TGraphErrors *g6m=gReadData("star_6m_200_0208_sc1.dat");
  //TGraphErrors *gct10=gReadData("star_ct10_200_0208_sc1.dat");
  //TGraphErrors *gct14=gReadData("star_ct14_200_0208_sc1.dat");
  //TGraphErrors *gct14sc05=gReadData("star_ct14_200_0208_sc05.dat");
  //TGraphErrors *gct14sc2=gReadData("star_ct14_200_0208_sc2.dat");
  TGraphErrors *gct10=gReadData(Form("star_ct10_200_05_%s_sc1.dat",R));
  */
  TGraphErrors *gct14R04=gReadData("star_ct14_200_05_R04_sc1.dat");
  TGraphErrors *gct14R04sc05=gReadData("star_ct14_200_05_R04_sc05.dat");
  TGraphErrors *gct14R04sc2=gReadData("star_ct14_200_05_R04_sc2.dat");

  //TGraphErrors *gct14sc05rel=Gdiv(gct14sc05,gct14);
  //TGraphErrors *gct14sc2rel=Gdiv(gct14sc2,gct14);

  TGraphErrors *gct14R02=gReadData("star_ct14_200_05_R02_sc1.dat");
  TGraphErrors *gct14R02sc05=gReadData("star_ct14_200_05_R02_sc05.dat");
  TGraphErrors *gct14R02sc2=gReadData("star_ct14_200_05_R02_sc2.dat");

  TGraphErrors *gct14R06=gReadData("star_ct14_200_05_R06_sc1.dat");
  TGraphErrors *gct14R06sc05=gReadData("star_ct14_200_05_R06_sc05.dat");
  TGraphErrors *gct14R06sc2=gReadData("star_ct14_200_05_R06_sc2.dat");

  TGraphErrors *gct14R02R04=Gdiv(gct14R02,gct14R04);
  TGraphErrors *gct14R06R04=Gdiv(gct14R06,gct14R04);

  gStyle->SetOptStat(0);
  TCanvas *cfig1=new TCanvas("cfig1","fig1",100,10,500,400);
  cfig1->SetMargin(0.12,0.05,0.12,0.05);
  cfig1->SetLogy(1);
  TH1D *h1=new TH1D("h1","",10,0,100);
  h1->SetXTitle("p_{T}  (GeV)");h1->SetYTitle("d#sigma/2#pid#etadp_{T}  (pb/GeV^{-1})");
  h1->SetMinimum(1e-6);h1->SetMaximum(1e8);
  h1->Draw();
  //g2006->SetMarkerColor(1);g2006->SetMarkerStyle(24);g2006->Draw("p");
  //g6m->SetLineColor(209);g6m->Draw("l");
  //gct10->SetLineColor(4);gct10->Draw("l");
  gct14R04->SetLineColor(2);gct14R04->Draw("l");
  keyLine(0.38,0.82,"NLO Vogelsang ct14 R=0.4",2,1);
  //gct14R04sc05->SetLineColor(2);gct14R04sc05->SetLineStyle(2);gct14R04sc05->Draw("l");
  //gct14R04sc2->SetLineColor(2);gct14R04sc2->SetLineStyle(2);gct14R04sc2->Draw("l");

  gct14R02->SetLineColor(1);gct14R02->Draw("l");
  keyLine(0.38,0.76,"NLO Vogelsang ct14 R=0.2",1,1);
  //gct14R02sc05->SetLineColor(1);gct14R02sc05->SetLineStyle(1);gct14R02sc05->Draw("l");
  //gct14R02sc2->SetLineColor(1);gct14R02sc2->SetLineStyle(1);gct14R02sc2->Draw("l");

  gct14R06->SetLineColor(209);gct14R06->Draw("l");
  keyLine(0.38,0.88,"NLO Vogelsang ct14 R=0.6",209,1);
  //gct14R06sc05->SetLineColor(209);gct14R06sc05->SetLineStyle(209);gct14R06sc05->Draw("l");
  //gct14R06sc2->SetLineColor(209);gct14R06sc2->SetLineStyle(209);gct14R06sc2->Draw("l");
  //keySymbol(0.5,0.88,"NLO STAR 2006 PRL",1,24);
  //keyLine(0.38,0.82,"NLO Vogelsang cteq6m",209,1);
  //keyLine(0.38,0.76,"NLO Vogelsang ct10",4,1);
  //keyLine(0.38,0.70,"NLO Vogelsang ct14",2,1);
  //cfig1->SaveAs(Form("NLO_%s.gif",R));
  //cfig1->SaveAs(Form("NLO_%s.pdf",R));

  TCanvas *cfig2=new TCanvas("cfig2","fig2",100,10,500,400);
  cfig2->SetMargin(0.12,0.05,0.12,0.05);
  TH1D *h2=new TH1D("h2","",10,0,100);
  h2->SetXTitle("p_{T}  (GeV)");h2->SetYTitle("Relative scale uncertainty");
  h2->SetMinimum(0.5);h2->SetMaximum(1.5);
  h2->Draw();
  gct14R02R04->SetLineColor(2);gct14R02R04->SetLineStyle(1);gct14R02R04->Draw("l");
  gct14R06R04->SetLineColor(209);gct14R06R04->SetLineStyle(1);gct14R06R04->Draw("l");
  keyLine(0.18,0.88,"NLO Vogelsang ct14 R=0.2/R=0.4",2,1);
  keyLine(0.18,0.82,"NLO Vogelsang ct14 R=0.6/R=0.4",209,1);
  TLine *line=new TLine();
  line->SetLineStyle(3);line->DrawLine(0,1,100,1);
  //cfig2->SaveAs(Form("NLO_%s_error.gif",R));
  //cfig2->SaveAs(Form("NLO_%s_error.pdf",R));
}
