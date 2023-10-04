void make_ratio(float R=0.2, bool central=1)
{
TString scent;
Int_t icent;
if(central){
	scent="0_10";
	icent=3;
}
else{
	scent="60_80";
	icent=2;
}

int jet_area_bin=2;
if(R>0.21) jet_area_bin=4;
else if(R>0.31) jet_area_bin=8;

TFile* inputfile_Delta_pt = new TFile(Form("Merge_Jet_200GeV_M42_I1_V29_R0%.0lf_%s_hjr%i.root",R*10,scent.Data(),icent),"OPEN");


//----------------------------------------------------------------------------------------
// Load Delta pt histograms

const Int_t N_jet_areas = 16;
const Int_t N_Psi_bins = 4;

TString HistName;

cout << "----------- Delta pt histograms -----------" << endl;
TH2D* h_Delta_pt_vs_embed_pt[N_jet_areas][2][N_Psi_bins]; //single particle embedded
TH2D* h_Delta_pt_vs_embed_pt_weight[N_jet_areas][2][N_Psi_bins];// single particle embedded, weighted according to v2

for(Int_t i_orig_smear = 0; i_orig_smear < 1; i_orig_smear++) //only 0 is important here
{
for(Int_t i = 0; i < N_jet_areas; i++)
{
for(Int_t i_Psi = 0; i_Psi < N_Psi_bins; i_Psi++)
{
HistName = "Delta_pt/h_Delta_pt_vs_embed_pt_";
HistName += i;
HistName += "_P_";
HistName += i_Psi;
if(i_orig_smear == 1) HistName += "_smear";
h_Delta_pt_vs_embed_pt[i][i_orig_smear][i_Psi] =(TH2D*)inputfile_Delta_pt->Get(HistName.Data());

HistName = "Delta_pt/h_Delta_pt_vs_embed_pt_weight_";
HistName += i;
HistName += "_P_";
HistName += i_Psi;
if(i_orig_smear == 1) HistName += "_smear";
h_Delta_pt_vs_embed_pt_weight[i][i_orig_smear][i_Psi]= (TH2D*)inputfile_Delta_pt->Get(HistName.Data());


}
}
}

// Integrate over Psi
TH2D* h2D_Delta_pt_vs_embed_pt_use =(TH2D*)h_Delta_pt_vs_embed_pt[jet_area_bin][0][0]->Clone("h2D_Delta_pt_vs_embed_pt_use");
TH2D* h2D_Delta_pt_vs_embed_pt_weight_use =(TH2D*)h_Delta_pt_vs_embed_pt_weight[jet_area_bin][0][0]->Clone("h2D_Delta_pt_vs_embed_pt_weight_use");

for(Int_t i_Psi = 1; i_Psi < N_Psi_bins; i_Psi++)
{
h2D_Delta_pt_vs_embed_pt_use->Add(h_Delta_pt_vs_embed_pt[jet_area_bin][0][i_Psi]);
h2D_Delta_pt_vs_embed_pt_weight_use->Add(h_Delta_pt_vs_embed_pt_weight[jet_area_bin][0][i_Psi]);
}

TH2D* hdiv=(TH2D*)h2D_Delta_pt_vs_embed_pt_weight_use->Clone("hdiv");
hdiv->Divide(h2D_Delta_pt_vs_embed_pt_use);

TCanvas *c1 =new TCanvas("c1","c1",10,10,1200,800);
  c1->cd();
  c1->SetGrid();
  h2D_Delta_pt_vs_embed_pt_use->SetTitle("uncorrected");
h2D_Delta_pt_vs_embed_pt_use->Draw("COLZ");

TCanvas *c2 =new TCanvas("c2","c2",10,10,1200,800);
  c2->cd();
  c2->SetGrid();
  h2D_Delta_pt_vs_embed_pt_weight_use->SetTitle("weighted");
h2D_Delta_pt_vs_embed_pt_weight_use->Draw("COLZ");

TCanvas *c3 =new TCanvas("c3","c3",10,10,1200,800);
  c3->cd();
  c3->SetGrid();
  hdiv->SetTitle("weighted/uncorrected");
hdiv->Draw("COLZ");

TFile *f = new TFile("v2corr.root","UPDATE");
f->cd();
TString name;
if(central) name=Form("dpt_v2corrVSuncorr_R0%.0lf_central",R*10);
else name=Form("dpt_v2corrVSuncorr_R0%.0lf_peripheral",R*10);
hdiv->Write(name);
f->Close();
}
//=========================================================

void compare_jan_alex(float R=0.2, bool central=1, TString suffix="")
{

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

	Color_t colorList[]={kRed,kGreen+3,kMagenta+2,kBlack,kBlue+3};
	Int_t marker1[]={20,21,22};
	Int_t marker2[]={24,25,26};
//----------------------------------
// Load Alex's Delta pt histograms
TString scent;
Int_t icent;
if(central){
	scent="0_10";
	icent=3;
}
else{
	scent="60_80";
	icent=2;
}

int jet_area_bin=2;
if(R>0.21) jet_area_bin=4;
else if(R>0.31) jet_area_bin=8;

TFile* inputfile_Delta_pt = new TFile(Form("Merge_Jet_200GeV_M42_I1_V29_R0%.0lf_%s_hjr%i.root",R*10,scent.Data(),icent),"OPEN");

const Int_t N_jet_areas = 16;
const Int_t N_Psi_bins = 4;

TString HistName;

//cout << "----------- Delta pt histograms -----------" << endl;
TH2D* h_Delta_pt_vs_embed_pt[N_jet_areas][2][N_Psi_bins]; //single particle embedded
TH2D* h_Delta_pt_vs_embed_pt_weight[N_jet_areas][2][N_Psi_bins];// single particle embedded, weighted according to v2

for(Int_t i_orig_smear = 0; i_orig_smear < 1; i_orig_smear++) //only 0 is important here
{
for(Int_t i = 0; i < N_jet_areas; i++)
{
for(Int_t i_Psi = 0; i_Psi < N_Psi_bins; i_Psi++)
{
HistName = "Delta_pt/h_Delta_pt_vs_embed_pt_";
HistName += i;
HistName += "_P_";
HistName += i_Psi;
if(i_orig_smear == 1) HistName += "_smear";
h_Delta_pt_vs_embed_pt[i][i_orig_smear][i_Psi] =(TH2D*)inputfile_Delta_pt->Get(HistName.Data());

HistName = "Delta_pt/h_Delta_pt_vs_embed_pt_weight_";
HistName += i;
HistName += "_P_";
HistName += i_Psi;
if(i_orig_smear == 1) HistName += "_smear";
h_Delta_pt_vs_embed_pt_weight[i][i_orig_smear][i_Psi]= (TH2D*)inputfile_Delta_pt->Get(HistName.Data());


}
}
}

// Integrate over Psi
TH2D* h2D_Delta_pt_vs_embed_pt_use =(TH2D*)h_Delta_pt_vs_embed_pt[jet_area_bin][0][0]->Clone("h2D_Delta_pt_vs_embed_pt_use");
TH2D* h2D_Delta_pt_vs_embed_pt_weight_use =(TH2D*)h_Delta_pt_vs_embed_pt_weight[jet_area_bin][0][0]->Clone("h2D_Delta_pt_vs_embed_pt_weight_use");

for(Int_t i_Psi = 1; i_Psi < N_Psi_bins; i_Psi++)
{
h2D_Delta_pt_vs_embed_pt_use->Add(h_Delta_pt_vs_embed_pt[jet_area_bin][0][i_Psi]);
h2D_Delta_pt_vs_embed_pt_weight_use->Add(h_Delta_pt_vs_embed_pt_weight[jet_area_bin][0][i_Psi]);
}

//------------------------------
//Load Jan's Delta pt histograms

  const Int_t nprobes=2;
  const Float_t embPt[nprobes]={1.0, 10.0};
  TString dir[2]={"peripheral","central"};
  TString input=Form("%s/histos_embeddedjet_R%.1lf%s.root",dir[central].Data(),R,suffix.Data());
  TFile *f =new TFile(input);
  f->cd(); 

  TH2D* hdeltapT_pt=(TH2D*)f->Get("delta_pt_BG_sp_0");

//-----------------------------
  
  TH1D* hdpt_alex[nprobes];
  TH1D* hdpt_jan[nprobes];
  TH1D* hdpt_ratio[nprobes];

  for(int i=0; i<nprobes; i++)
  {
	  Float_t pTemb=embPt[i];
	  Int_t xbin=hdeltapT_pt->GetXaxis()->FindBin(pTemb);
	  TString hname=Form("hdpt_jan_%i",i);
	  hdpt_jan[i]=(TH1D*) hdeltapT_pt->ProjectionY(hname,xbin,xbin);
	  hdpt_jan[i]->SetLineColor(colorList[i]);
	  hdpt_jan[i]->SetMarkerColor(colorList[i]);
	  hdpt_jan[i]->SetMarkerStyle(marker1[i]);
	  hdpt_jan[i]->Scale(1./hdpt_jan[i]->Integral("width"));
	  
	  xbin=h2D_Delta_pt_vs_embed_pt_use->GetXaxis()->FindBin(pTemb);
	  hname=Form("hdpt_alex_%i",i);
	  hdpt_alex[i]=(TH1D*) h2D_Delta_pt_vs_embed_pt_use->ProjectionY(hname,xbin,xbin);
	  hdpt_alex[i]->SetLineColor(colorList[i]);
	  hdpt_alex[i]->SetMarkerColor(colorList[i]);
	  hdpt_alex[i]->SetMarkerStyle(marker2[i]);
	  hdpt_alex[i]->Scale(1./hdpt_alex[i]->Integral("width"));

	  hname=Form("hdpt_ratio_%i",i);
	  hdpt_ratio[i]=(TH1D*)hdpt_alex[i]->Clone(hname);
	  hdpt_ratio[i]->Reset("MICE");
	  for(int bin=1; bin<=hdpt_ratio[i]->GetNbinsX(); bin++)
	  {
				 Double_t pT=hdpt_ratio[i]->GetBinCenter(bin);
				 Int_t bin_jan=hdpt_jan[i]->GetXaxis()->FindBin(pT);
				 Double_t yalex=hdpt_alex[i]->GetBinContent(bin);
				 Double_t yjan=hdpt_jan[i]->GetBinContent(bin_jan);
				 Double_t ratio=0;
				 if(yalex>0)ratio=yjan/yalex;
				 hdpt_ratio[i]->SetBinContent(bin,ratio);
	  }
  }
  
  Float_t spectraXmin=-15;
  if(!central)spectraXmin=-5;
  Float_t spectraXmax=20;
  Float_t spectraYmin=1E-4;
  Float_t spectraYmax=1;
  //if(!central)spectraYmax=0.45;
  if(!central)spectraYmax=5;
  
  TH1 *frame = new TH1I("frame", "", 1000, -100, +100);
  
  TCanvas *c1 =new TCanvas("c1","c1",10,10,1200,800);
  c1->cd();
  c1->SetGrid();
  Double_t eps=0.02;
  TPad* p1 = new TPad("p1","p1",0,0.35-eps,1,1,0); p1->Draw();
  p1->SetBottomMargin(eps);
  p1->SetGrid();
  p1->SetLogy();
  
  TPad* p2 = new TPad("p2","p2",0,0,1,0.35*(1.-eps),0); p2->Draw();
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0.25);
  p2->SetGrid();
  p2->SetFillColor(0);
  p2->SetFillStyle(0);
  
  p1->cd();

  frame->SetTitle(Form("R=%.1lf, %s collisions",R,dir[central].Data()));
  frame->GetXaxis()->SetTitle("#delta p_{T} (GeV/c)");
  frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
  frame->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
  frame->DrawClone();
   for(int i=0; i<nprobes; i++)
  {
	hdpt_jan[i]->Draw("same");
	hdpt_alex[i]->Draw("same");
  }
  
   //Legend
  TLegend *legspectra = new TLegend(0.6467, 0.60, 0.89, 0.90);
  legspectra->SetTextSize(0.03);
  legspectra->SetFillStyle(0);
  legspectra->SetBorderSize(0);
    for(int i=0; i<nprobes; i++)
  {
  legspectra->AddEntry(hdpt_jan[i], Form("Jan: #delta p_{T}, p{T,emb}=%.1lf GeV/c",embPt[i]), "lp");
  legspectra->AddEntry(hdpt_alex[i], Form("Alex: #delta p_{T}, p{T,emb}=%.1lf GeV/c",embPt[i]), "lp");
  }
	legspectra->Draw("same");

	p2->cd();
 	frame->GetYaxis()->SetTitle("Jan/Alex");
	frame->SetTitle("");
   frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
   frame->GetYaxis()->SetRangeUser(0.2, 1.5);
	TAxis *yaxis = frame->GetYaxis();
   yaxis->SetLabelSize(0.08);
	yaxis->SetTitleSize(0.1);
	yaxis->SetTitleOffset(0.3);
   TAxis *xaxis = frame->GetXaxis();
   xaxis->SetLabelSize(0.08);
	xaxis->SetTitleSize(0.1);
	frame->DrawClone();
	
	TLine *one = new TLine(spectraXmin, 1, spectraXmax, 1);
  one->SetLineWidth(2);
  one->SetLineStyle(2);
  one->SetLineColor(kBlack);
  one->DrawClone("same");

	for(int i=0; i<nprobes; i++)
	  {
		  hdpt_ratio[i]->Draw("same");
	  }
	  c1->SaveAs(Form("jan_vs_alex_deltapT_%s_R%.1lf.gif",dir[central].Data(),R));
}
