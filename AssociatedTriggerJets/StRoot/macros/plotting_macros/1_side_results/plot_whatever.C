//available plotting functions:
//plot_measured_R() - plot measured spectra for different R values
//plot_measured_pTlead() - plot measured spectra for different pTleading values
//plot_measured_centrality() - plot measured spectra for different centralities
//plot_measured_vs_toymodel() - compare star data with toymodel
//plot_measured_vs_pythia() - compare star data with PYTHIA
//plot_epsilon() - plot jet reconstruction efficiency
//plot_dpT() - plot deltapT distributions
//plot_area_and_rho() - plot jet area and rho distribution
//plot_RM() - plot response matrix
//show_bbc_ntr() - show <N> tracks vs BBC rate (=pileup removal proof)

//include tracking efficiency function for comparison with jet reconstruction efficiency
#include "../3_additional_studies/tracking_efficiency/efficiency_jan/include/eff_functions.C"

TString G_system_info[]={"Central (0-10%)", "Peripheral (60-80%)", "pp"};
enum systems { cent, peri, pp };

double new_bins1[48]={-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,35,40};
double new_bins2[33]={-20,-15,-12,-10,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,15,17,20,25,30,35,40};

//1D histograms visual settings
//lines
const float G_linewidth=2.0;
Color_t G_colorList[]={kRed,kGray+3,kBlue,kBlack};
int G_lineStyle[]={1,3,2,4};

//fills
Color_t G_colorFill[]={kGray,kGray,kGray,kGray};
int G_styleFill[]={0,3004,0,0};
//markers
const float G_marksize=1.0;
int G_markerList[]={20,21,22,23,29,33,34};

//draw options
TString G_histoDraw[]={"histo same","histo same","histo same","histo same"};
TString G_legDraw[]={"l","f","l","l"};

float G_latex_sz=0.032;
Color_t G_latex_cl=kGray+1;
//TString G_label="THIS THESIS";
/*
 *	TLatex *latex = new TLatex();
 *  latex->SetNDC();
 *  latex->SetTextSize(0.035);
 *	latex->DrawLatex(0.4, 0.3,G_label);*/


//=================================================================
//plot measured spectra for different Rs
//=================================================================
plot_measured_R(float pTlead=5.0, systems system=cent, TString version="GPC2" ,TString label="THIS THESIS", TString ext = "pdf")
{ 
	//gStyle->SetPadTickY(0);
    TString outdir = Form("../../plotting_out/obr/%s/uncorrected",version.Data()); //figure output directory
	
	
	const float marker_size=1.5;
	const float line_width=2;
	
	const Color_t colorList[]={kBlack,kRed,kGreen+3,kMagenta+2,kBlue,kOrange+2,kYellow+2,kBlue-4,kOrange-2,kGreen,kRed+2,7,kGreen+4,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,28};
	const int marker_list[]={21,20,22,30};
	const int marker[]={20,21,22,23};
	
	//canvas size
	const int can_x=1000; //1600
	const int can_y=750; //900
	
	TString trigger="MB";
	if(system==peri) trigger+="_peripheral";
	else if(system==cent) trigger+="_central";
	else if(system==pp) trigger+="HT_pp";
	TString cent_suff[]={"","_peri",""};
	
	float spectraXmin=-25;
	float spectraXmax=45;
	float spectraYmin=1E-9;
	float spectraYmax=1E3;
	
	if(system==pp)spectraYmin=1E-10;
	if(system==peri)
	{
		spectraYmin=1E-9;
		spectraYmax=1E3;
		spectraXmin=-25;
		spectraXmax=45;
	}
	
	int binmerge=2; //rebin
	
	const int nhistos=3;
	const double R[nhistos]={0.2,0.3,0.4};
	const double Acut[nhistos]={0.07,0.3,0.4};
	TH1D* hmeasured[nhistos];
	TH1D* hplot[nhistos];
	
	int mevents=0;
	
	for(int i=0; i<nhistos; i++){
		TString str = Form("../../plotting_out/root/%s/%s/histos_inclusivejet_normal.root",trigger.Data(),version.Data());
		TFile *f = new TFile(str.Data(), "OPEN");
		TH1I* hevents;
		hevents= (TH1I*) f->Get("hevents");
		int nevents=hevents->GetBinContent(2);
		if(i==0)mevents=(int)nevents/1000000;
		double hole=0;//(1.0/12.)+((R[i]*2.)/(2.*TMath::Pi()));
		double scale_jets = 1./(2*(1.-R[i])*2.*TMath::Pi()*nevents*(1-hole));
		if(system!=pp)
		{
			TH2D* histo2d=(TH2D*) f->Get(Form("hpT_pTlead_R0%.0lf",R[i]*10));
			hmeasured[i]=histo2d->ProjectionX(Form("histo1d_%i",i),histo2d->GetYaxis()->FindBin(pTlead),histo2d->GetYaxis()->GetNbins());
		}
		else
		{
			hmeasured[i]=(TH1D*) f->Get(Form("hpT_pTl%0.lf_R0%.0lf",pTlead,R[i]*10));
		}
		//hmeasured[i]->SetMarkerStyle(marker_list[i+1]);
		hmeasured[i]->SetMarkerStyle(marker[i]);
		hmeasured[i]->SetLineColor(colorList[i]);
		hmeasured[i]->SetLineWidth(line_width);
		hmeasured[i]->SetMarkerColor(colorList[i]);
		hmeasured[i]->SetMarkerSize(marker_size);
		hmeasured[i]->Sumw2();
		hplot[i]=(TH1D*) hmeasured[i]->Rebin(47,Form("hplot_%i",i),new_bins1);
		hplot[i]->Scale(scale_jets,"width");
	}
	
	TH1 *frame = new TH1I("frame", "", 1000, -100, +100);
	
	TCanvas *spectra = new TCanvas("cspectra","cspectra",10,10,can_x,can_y);
	cspectra->cd();
	//cspectra->SetGrid();
	cspectra->SetLogy();
	frame->GetXaxis()->SetTitle("p_{T, jet}^{reco, ch} (GeV/#it{c})");
	frame->GetYaxis()->SetTitle("1/N_{events} 1/2#pi d^{2}N/dp_{T}^{ch}d#eta (GeV/#it{c})^{-1}");
	frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
	frame->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
	//frame->GetYaxis()->SetRangeUser(1E-3, 1E5);
	str="";
	frame->SetTitle(str);
	frame->DrawCopy("");
	for(int i=0; i<nhistos; i++){
		hplot[i]->DrawCopy("E same");
	}
	TLegend *model_info = new TLegend(0.20, 0.65, 0.35, 0.900);
	model_info->SetTextSize(0.047);
	model_info->SetFillStyle(0);
	model_info->SetBorderSize(0);
	model_info->SetMargin(0.05);
	//model_info->SetHeader("Run 11 Au+Au #sqrt{s_{NN}} = 200 GeV");
	model_info->SetHeader("Au+Au #sqrt{s_{NN}} = 200 GeV");
	model_info->AddEntry("", G_system_info[system], "");
	model_info->AddEntry("", "anti-k_{T}", "");
	//model_info->AddEntry("", Form("%iM events",mevents), "");
	if(pTlead>0) model_info->AddEntry("",Form("p_{T,lead}^{min} = %.0lf GeV/#it{c}",pTlead),"");
	model_info->DrawClone("same");
	
	TLatex *latex = new TLatex();
	latex->SetNDC();
	latex->SetTextSize(G_latex_sz);
	latex->SetTextColor(G_latex_cl);
	latex->DrawLatex(0.4, 0.2,label);
	
	TLegend *legspectra = new TLegend(0.54, 0.60, 0.85, 0.90);
	legspectra->SetTextSize(0.045);
	legspectra->SetFillStyle(0);
	legspectra->SetBorderSize(0);
	for(int i=0; i<nhistos; i++){
		TString sleg=Form("R = %.1lf, A > %.1lf Sr",R[i],Acut[i]);
		if(i==0)sleg=Form("R = %.1lf, A > %.2lf Sr",R[i],Acut[i]);
		legspectra->AddEntry(hmeasured[i],sleg , "lp");
	}
	legspectra->DrawClone("same");
	
	str = Form("%s/pTreco_R-dep_pTl%.0lf%s.%s", outdir.Data(),pTlead,cent_suff[system].Data(),ext.Data());
	cspectra->SaveAs(str.Data());
	
}


//=================================================================
//plot measured spectra for different pTleadings
//=================================================================
plot_measured_pTlead(float R=0.3, systems system=cent, TString version="GPC2", TString label="THIS THESIS", TString ext = "pdf")
{  
	TString outdir = Form("../../plotting_out/obr/%s/uncorrected",version.Data()); //figure output directory
	
	float binmerge=4; //how many bins to merge
	
	const float line_width=2;
	const float marker_size=1.0;
	const Color_t colorList[]={kBlack,kMagenta+2,kBlue,kRed,kGreen+3,kOrange+2,kYellow+2,kBlue-4,kOrange-2,kGreen,kRed+2,7,kGreen+4,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,28};
	const int marker_list[]={20,21,22,23};
	
	//canvas size
	const int can_x=1000; //1600
	const int can_y=750; //900
	
	
	float spectraXmin=-30;
	float spectraXmax=40;
	float spectraYmin=1E-9;
	if(system==pp)spectraYmin=1E-10;
	float spectraYmax=1E2;
	
	const int nhistos=4;
	const float pTlead[nhistos]={0, 3.,5.,7.};
	TH1D* hmeasured[nhistos];
	TH1D* hplot[nhistos];
	TH1D* hdiv[nhistos];
	int mevents=0;
	
	int ridx=calc_R_index(R);
	
	const int nPanels=4; //maximum number of columns
	const int nRows=3; //maximum number of rows
	//in which panels do we want to show info, legend and labels?
	bool showInfo[nRows][nPanels]={/*row1*/{1,1,1,0},
	/*row2*/{1,1,1,0},
	/*row3*/{1,1,1,0}};
	bool showLeg[nRows][nPanels]={/*row1*/{1,1,1,0},
	/*row2*/{1,1,1,0},
	/*row3*/{1,1,1,0}};
	bool showLabel[nRows][nPanels]={/*row1*/{1,1,1,0},
	/*row2*/{1,1,1,0},
	/*row3*/{1,1,1,0}};
	bool showWatermark[nRows][nPanels]={/*row1*/{1,1,1,0},
	/*row2*/{1,1,1,0},
	/*row3*/{1,1,1,0}};
	short panelType=0; //what does show different panels? 0: different Rs | 1: different centralities
	int pidx=(panelType==0)? ridx : system; //panel index
	int rowidx=(panelType==1)? ridx : system; //row index
	
	//set description of panels and descriptor in information legend based on the panel type
	TString descR=Form("R = %.1lf",R); 
	TString descC=G_system_info[system];
	/*TString descPanel=descR;//(panelType==0)? descR : descC;
	 *	TString descInfo=descC;//(panelType==1)? descR : descC;
	 *	float panx[]={0.75,0.25,0.25,0.25}; //x-position of the panel desc.
	 *	float pany=0.8;//y-position of the panel desc.*/
	
	
	TString trigger="MB";
	if(system==peri) trigger+="_peripheral";
	else if(system==cent) trigger+="_central";
	else if(system==pp) trigger+="HT_pp";
	TString cent_suff[]={"","_peri",""};
	
	
	
	for(int i=0; i<nhistos; i++){
		TString str = Form("../../plotting_out/root/%s/%s/histos_inclusivejet_normal.root",trigger.Data(),version.Data());
		TFile *f = new TFile(str.Data(), "OPEN");
		TH1I* hevents;
		hevents= (TH1I*) f->Get("hevents");
		int nevents=hevents->GetBinContent(2);
		if(i==0)mevents=(int)nevents/1000000;
		double hole=0;
		double scale_jets = 1./(2*(1.-R)*2.*TMath::Pi()*nevents*(1-hole));
		if(system!=pp){
			TH2D* histo2d=(TH2D*) f->Get(Form("hpT_pTlead_R0%.0lf",R*10));
			hmeasured[i]=histo2d->ProjectionX(Form("histo1d_%i",i),histo2d->GetYaxis()->FindBin(pTlead[i]),histo2d->GetYaxis()->GetNbins());
		}
		else
		{
			hmeasured[i]=(TH1D*) f->Get(Form("hpT_pTl%0.lf_R0%.0lf",pTlead[i],R*10));
		}
		hmeasured[i]->SetMarkerStyle(marker_list[i]);
		//hmeasured[i]->SetMarkerStyle(30);
		hmeasured[i]->SetLineColor(colorList[i]);
		hmeasured[i]->SetLineWidth(line_width);
		hmeasured[i]->SetMarkerColor(colorList[i]);
		hmeasured[i]->SetMarkerSize(marker_size);
		hmeasured[i]->Sumw2();
		hplot[i]=(TH1D*) hmeasured[i]->Rebin(47,Form("hplot_%i",i),new_bins1);
		hplot[i]->Scale(scale_jets,"width");
		if(i>0)
		{
			hdiv[i]=(TH1D*)hplot[i]->Clone(Form("hdiv_%i",i));
			hdiv[i]->Divide(hplot[0]);
		}
	}
	
	TH1 *frame = new TH1I("frame", "", 1000, -100, +100);
	
	TCanvas *spectra = new TCanvas("cspectra","cspectra",10,10,can_x,can_y);
	cspectra->SetLogy();
	cspectra->cd();
	
	 	/*double eps=0.02;
	 	TPad* p1 = new TPad("p1","p1",0,0.35-eps,1,1,0); p1->Draw();
	 	p1->SetBottomMargin(eps);
	    p1->SetLogy();
	  
	    TPad* p2 = new TPad("p2","p2",0,0,1,0.35*(1.-eps),0); p2->Draw(); 
	    p2->SetTopMargin(0);
	    p2->SetBottomMargin(0.25);
	    p2->SetFillColor(0);
	    p2->SetFillStyle(0);
	    p2->SetLogy();
	    
	 	p1->cd();*/
	frame->GetXaxis()->SetTitle("p_{T, jet}^{reco, ch} (GeV/#it{c})");
	frame->GetYaxis()->SetTitle("1/N_{events} 1/2#pi d^{2}N/dp_{T}^{ch}d#eta (GeV/#it{c})^{-1}");
	frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
	frame->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
	//frame->GetYaxis()->SetRangeUser(1E-3, 1E5);
	str="";
	frame->SetTitle(str);
	frame->DrawCopy("");
	for(int i=0; i<nhistos; i++){
		hplot[i]->DrawCopy("E same");
		/*int pbin=hplot[i]->FindBin(31);
		double counts=hplot[i]->GetBinContent(pbin);
		cout<<"pTlead:"<<pTlead[i]<<" bin counts in bin 30-35:"<<counts<<endl;*/
	}
	TLegend *model_info = new TLegend(0.20, 0.6, 0.35, 0.900);
	model_info->SetTextSize(0.045);
	model_info->SetFillStyle(0);
	model_info->SetBorderSize(0);
	model_info->SetMargin(0.05);
	if(showInfo[rowidx][pidx])
	{
		if(system!=pp) model_info->SetHeader("Au+Au #sqrt{s_{NN}} = 200 GeV");
		else model_info->SetHeader("Run 12 p+p #sqrt{s_{NN}} = 200 GeV");
	}
	else model_info->SetHeader(" ");
	model_info->AddEntry("", descC, "");
	model_info->AddEntry("", Form("anti-k_{T}, %s",descR.Data()), "");
	//model_info->AddEntry("", Form("%iM events",mevents), "");
	if(showLabel[rowidx][pidx]) model_info->DrawClone("same");
	
	TLatex *latex = new TLatex();
	latex->SetNDC();
	latex->SetTextSize(G_latex_sz);
	latex->SetTextColor(G_latex_cl);
	if(showWatermark[rowidx][pidx]) latex->DrawLatex(0.5, 0.2,label);
	
	TLatex *latexP = new TLatex();
	latexP->SetNDC();
	latexP->SetTextSize(0.04);
	//if(showLabel[rowidx][pidx])latexP->DrawLatex(panx[ridx],pany,descPanel);
	
	TLegend *legspectra = new TLegend(0.65, 0.55, 0.9, 0.90);
	legspectra->SetTextSize(0.042);
	legspectra->SetFillStyle(0);
	legspectra->SetBorderSize(0);
	legspectra->SetHeader("p_{T,lead}^{min} = ");
	for(int i=0; i<nhistos; i++){
	  // legspectra->AddEntry(hmeasured[i], Form("%.1lf GeV/#it{c}",pTlead[i]), "lp"); //Jana
	  legspectra->AddEntry(hmeasured[i], Form("%.0lf GeV/#it{c}",pTlead[i]), "lp"); //Jana
	}
	if(showLeg[rowidx][pidx]) legspectra->DrawClone("same");
	/*
	 	p2->cd();
	 	frame->GetXaxis()->SetTitle("p_{T, reco}^{charged} (GeV/#it{c})");
	 	frame->GetYaxis()->SetTitle("ratio");
	 	frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
	 	frame->GetYaxis()->SetRangeUser(1E-2, 1.5);
	 	frame->SetTitleSize(0.1,"Y");
	 	frame->SetTitleOffset(0.3,"Y");
	 	frame->SetTitleSize(0.1,"X");
	 	frame->SetTitleOffset(0.95,"X");
	 	frame->SetLabelSize(0.06,"X");
	 	frame->SetLabelSize(0.06,"Y");
	   //frame->GetYaxis()->SetRangeUser(1E-3, 1E5);
	 	frame->SetTitle("");
	 	frame->DrawCopy("");
	 	for(int i=1; i<nhistos; i++)
	 	{
	 		hdiv[i]->DrawCopy("esame");
        }
TLegend *legspectra = new TLegend(0.2, 0.4, 0.6, 0.90);
legspectra->SetTextSize(0.08);
legspectra->SetFillStyle(0);
legspectra->SetBorderSize(0);
for(int i=1; i<nhistos; i++)
{
legspectra->AddEntry(hdiv[i], Form("p_{T}^{lead}>%.1lf/p_{T}^{lead}>0",pTlead[i]), "lp");
}
legspectra->DrawClone("same");

TLine *one = new TLine(0, 1, spectraXmax, 1);
one->SetLineWidth(2);
one->SetLineStyle(2);
one->SetLineColor(kGray+1);
one->DrawClone("same");
*/

	str = Form("%s/pTreco_pTlead-dep_R0%.0lf%s.%s", outdir.Data(),R*10,cent_suff[system].Data(),ext.Data());
	cspectra->SaveAs(str.Data());
	
}

//=================================================================
//plot measured spectra for different centralities
//=================================================================
plot_measured_centrality(float R=0.3,float pTlead=0.0, TString version="GPC2")
{ 
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	gStyle->SetOptDate(0);
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
	
	/*
	 *  const int markerU=21; //marker style for unfolding
	 *  const int markerP=20; //marker style for prior
	 *  const int markerM=29; //marker style for measured
	 *  const int markerB=22; //marker style for backfolded
	 *  const int markerU2=25; //marker style for unfolding2
	 *  const int markerP2=24; //marker style for prior2
	 *  const int markerM2=30; //marker style for measured2
	 *  const int markerB2=26; //marker style for backfolded2
	 */
	
	TString outdir = "../../plotting_out/obr/intersteps/reco"; //figure output directory
	TString ext = "gif";
	
	const float marker_size=1.0;
	const float line_width=2;
	
	const Color_t colorList[]={kBlack,kGreen+3,kBlue,kRed,kMagenta+2,kOrange+2,kYellow+2,kBlue-4,kOrange-2,kGreen,kRed+2,7,kGreen+4,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,28};
	const int marker_list[]={21,20,22,30};
	
	//canvas size
	const int can_x=900; //1600
	const int can_y=460; //900
	
	const  float spectraXmin=-30;
	const float spectraXmax=40;
	const float spectraYmin=1E-8;
	const float spectraYmax=1E1;
	
	const int nhistos=2;
	const float TAA[nhistos]={22.75,0.49};
	const TString trigger[nhistos]={"MB_central","MB_peripheral"};
	const TString centrality[nhistos]={"central","peripheral"};
	const TString suffix[nhistos]={"",""};
	TH1D* hmeasured[nhistos];
	
	int mevents=0;
	
	for(int i=0; i<nhistos; i++){
		TString str = Form("../../plotting_out/root/%s/%s/histos_inclusivejet_R%.1lf%s.root",trigger[i].Data(),version.Data(),R,suffix[i].Data());
		TFile *f = new TFile(str.Data(), "OPEN");
		TH1I* hevents;
		hevents= (TH1I*) f->Get("hevents");
		int nevents=hevents->GetBinContent(2);
		if(i==0)mevents=(int)nevents/1000000;
		double hole=0;//(1.0/12.)+((R[i]*2.)/(2.*TMath::Pi()));
		double scale_jets = 1./(2*(1.-R)*2.*TMath::Pi()*nevents*(1-hole));
		TH2D* histo2d=(TH2D*) f->Get("hpT_pTlead");
		hmeasured[i]=histo2d->ProjectionX(Form("histo1d_%i",i),histo2d->GetYaxis()->FindBin(pTlead),histo2d->GetYaxis()->GetNbins());
		//hmeasured[i]->SetMarkerStyle(marker_list[i+1]);
		hmeasured[i]->SetMarkerStyle(30);
		hmeasured[i]->SetLineColor(colorList[i]);
		hmeasured[i]->SetLineWidth(line_width);
		hmeasured[i]->SetMarkerColor(colorList[i]);
		hmeasured[i]->SetMarkerSize(marker_size);
		hmeasured[i]->Rebin(2);
		hmeasured[i]->Scale(scale_jets/TAA[i],"width");
	}
	
	TH1 *frame = new TH1I("frame", "", 1000, -100, +100);
	
	TCanvas *spectra = new TCanvas("cspectra","cspectra",10,10,can_x,1.5*can_y);
	
	double eps=0.02;
	TPad* p1 = new TPad("p1","p1",0,0.35-eps,1,1,0); 
	p1->Draw();
	p1->SetBottomMargin(eps);
	p1->SetGrid();
	p1->SetLogy();
	
	TPad* p2 = new TPad("p2","p2",0,0,1,0.35*(1.-eps),0); 
	p2->Draw(); 
	p2->SetTopMargin(0);
	p2->SetBottomMargin(0.25);
	p2->SetGrid();
	p2->SetLogy();
	p2->SetFillColor(0);
	p2->SetFillStyle(0);
	
	p1->cd();
	
	frame->GetXaxis()->SetTitle("p_{T, reco}^{charged} (GeV/#it{c})");
	frame->GetYaxis()->SetTitle("1/T_{AA} 1/N_{events} 1/2#pi d^{2}N/dp_{T}^{ch}d#eta (GeV/#it{c})^{-1}");
	frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
	frame->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
	//frame->GetYaxis()->SetRangeUser(1E-3, 1E5);
	str="";
	frame->SetTitle(str);
	frame->DrawCopy("");
	for(int i=0; i<nhistos; i++){
		hmeasured[i]->DrawCopy("E same");
	}
	TLegend *model_info = new TLegend(0.20, 0.48, 0.35, 0.900);
	model_info->SetTextSize(0.035);
	model_info->SetFillStyle(0);
	model_info->SetBorderSize(0);
	model_info->SetMargin(0.05);
	model_info->SetHeader("Run 11 Au+Au #sqrt{s_{NN}} = 200 GeV");
	model_info->AddEntry("", "anti-k_{T}", "");
	model_info->AddEntry("", Form("%iM events",mevents), "");
	/*model_info->AddEntry("", "p_{T}^{const} > 0.2 GeV/#it{c}", "");
	 *  model_info->AddEntry("", Form("p_{T}^{leading} > %.1lf GeV/#it{c}",pTthresh), "");
	 *   
	 *  if(R>0.35) model_info->AddEntry("", "A_{reco jet} > 0.4 sr", "");
	 *  else if(R>0.25) model_info->AddEntry("", "A_{reco jet} > 0.2 sr", "");
	 *  else model_info->AddEntry("", "A_{reco jet} > 0.09 sr", "");*/
	model_info->DrawClone("same");
	
	TLatex *latex = new TLatex();
	latex->SetNDC();
	latex->SetTextSize(0.035);
	//latex->DrawLatex(0.4, 0.3,"STAR Preliminary");
	
	TLegend *legspectra = new TLegend(0.6467, 0.60, 0.89, 0.90);
	legspectra->SetTextSize(0.04);
	legspectra->SetFillStyle(0);
	legspectra->SetBorderSize(0);
	for(int i=0; i<nhistos; i++){
		legspectra->AddEntry(hmeasured[i], centrality[i].Data(), "lp");
	}
	legspectra->DrawClone("same");
	
	
	p2->cd();
	
	TH1D *hdiv=(TH1D*)hmeasured[0]->Clone("hdiv");
	hdiv->Divide(hmeasured[1]);
	
	frame->GetXaxis()->SetTitle("p_{T, reco}^{charged} (GeV/#it{c})");
	frame->GetYaxis()->SetTitle("central/peripheral");
	frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
	frame->GetYaxis()->SetRangeUser(1E-1, 1E2);
	frame->SetTitleSize(0.1,"Y");
	frame->SetTitleOffset(0.3,"Y");
	frame->SetTitleSize(0.1,"X");
	frame->SetTitleOffset(0.95,"X");
	frame->SetLabelSize(0.06,"X");
	frame->SetLabelSize(0.06,"Y");
	//frame->GetYaxis()->SetRangeUser(1E-3, 1E5);
	str="";
	frame->SetTitle(str);
	frame->DrawCopy("");
	hdiv->DrawCopy("esame");
	
	TLine *one = new TLine(spectraXmin, 1, spectraXmax, 1);
	one->SetLineWidth(2);
	one->SetLineStyle(2);
	one->SetLineColor(kBlack);
	one->DrawClone("same");
	
	str = Form("%s/spectra_pTcorr_cent-dep.%s", outdir.Data(),ext.Data());
	cspectra->SaveAs(str.Data());
	
}



//========================================================================
//STAR vs Toymodel
//========================================================================

plot_measured_vs_toymodel(float pTlead_toy=0.0, systems system=cent,TString version="GPC2",TString version_toy="GPC2"/*"newBinning2"*/, TString ext = "pdf")
{ 
	
	//Graphics settings
	//using variables from rootlogon.C
	//--------------------------------------------------------
	int graphics=0; //gif, png
	if(ext=="pdf" || ext=="ps" || ext=="eps")  graphics=1;
	
	gStyle->SetGridColor(gridColorG[graphics]);
	gStyle->SetHatchesLineWidth(hatchesLineWidthG[graphics]);
	gStyle->SetGridWidth(gridLineWidthG[graphics]);
	gStyle->SetGridColor(gridColorG[graphics]);
	gStyle->SetGridStyle(gridStyleG[graphics]);
	const int lineLineWidth=lineLineWidthG[graphics];
	const int graphLineWidth=graphLineWidthG[graphics];
	const int fillStyle=fillStyleG[graphics];
	bool showGrid=0;
	//---------------------------------------------------------
	
	const float marker_size=1.0;
	const Color_t colorList[]={kBlack,kRed,kGreen+3,kMagenta+2,kBlue,kOrange+2,kGreen+3,kMagenta-2,kOrange-2,kGreen,kRed+2,7,kGreen+4,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,28};
	const int marker_list[]={30,24,25,26,27};
	
	//canvas size
	const int can_x=700; //1600
	const int can_y=650; //900
	
	bool showHardJet=0; //show also reconstructed hard jet distribution
	
	double RAAtoy=0.5;
	double TAA=22.75;
	TString trigger="MB";
	TString cent_suff="_central";
	if(system==peri)
	{
		cent_suff="_peripheral";
		RAAtoy=0.7;
		TAA=0.5;
	}
	trigger+=cent_suff;
	
	//Which R parameters to show
	const int nR=2;
	float Rs[nR]={0.2,0.4};
	
	const  float spectraXmin=-25;
	const float spectraXmax=42;
	const float spectraYmin=1E-8;
	const float spectraYmax=1E2;
	bool doRebin=1;
	const int rebin_bins=2;
	TString ytitle="1/N_{events} 1/2#pi d^{2}N_{jet}/dp_{T}^{ch}d#eta (GeV/#it{c})^{-1}";
	
	TString outdir = Form("../../plotting_out/obr/%s/toy",version.Data()); //figure output directory
	
	const int nhistos=4;
	const int nhistos_toy=nhistos; 
	const float pTlead[nhistos]={0.0,3.0,5.0,7.0}; 	
	bool compare[nhistos]={1,0,1,0}; //for which pTlead cuts do we want to show the ratio of toymodel/STAR
	TH1D* hmeasured[nhistos][nR];
	TH1D* hplot[nhistos][nR];
	TH1D* hmeasured_toy[nhistos_toy][nR];
	TH1D* hplot_toy[nhistos_toy][nR];
	TH1D* hmeasured_htoy[nR];
	TFile *f[nR];
	TFile *ftoy[nR];
	TFile *ftoy2[nR];
	TFile *fhtoy[nR];
	TH1D* htoydatarat[nhistos][nR];
	double scale_jets[nR];
	double scale_jets_toy_rec[nR];
	double scale_jets_toy_gen[nR];
	double scale_jets_toy_hard[nR];
	float mevents_toy;
	float mevents_htoy;
	float mevents;
	
	int icompare=0; //number of STAR data histogram to compare with toymodel histogram 
	for(int i=0; i<nhistos; i++)
	{
		float epsilon=0.001;
		if((pTlead[i]+epsilon)>pTlead_toy && (pTlead[i]-epsilon)<pTlead_toy) 
		{
			icompare=i;
			break;
		}
	}
	
	//int mevents=0;
	for (int r=0;r<nR;r++)
	{
			
		if(system!=cent)TAA=0.5;
		
		TString str = Form("../../plotting_out/root/%s/%s/histos_inclusivejet_normal.root",trigger.Data(),version.Data());
		f[r] = new TFile(str.Data(), "OPEN");
		TH1I* hevents;
		hevents= (TH1I*) f[r]->Get("hevents");
		int nevents=hevents->GetBinContent(2);
		if(r==0)mevents=(float)nevents/1000000;
		double hole=0;//1/(1.0/12.)+((Rs[r]*2.)/(2.*TMath::Pi()));
		scale_jets[r] = 1./(2*(1.-Rs[r])*2.*TMath::Pi()*nevents*(1-hole));
		//TH2D* histo2d=(TH2D*) f->Get("hpT_pTlead_nobadsecedge1R");
		TH2D* histo2d=(TH2D*) f[r]->Get(Form("hpT_pTlead_R0%.0lf",Rs[r]*10));
		for(int i=0; i<nhistos; i++){
			hmeasured[i][r]=histo2d->ProjectionX(Form("histo1d_%i_%i",i,r),histo2d->GetYaxis()->FindBin(pTlead[i]),histo2d->GetYaxis()->GetNbins());
			hmeasured[i][r]->SetMarkerStyle(marker_list[i]);
			//hmeasured[i][r]->SetMarkerStyle(30);
			hmeasured[i][r]->SetLineColor(colorList[i]);
			hmeasured[i][r]->SetLineWidth(graphLineWidth);
			hmeasured[i][r]->SetMarkerColor(colorList[i]);
			hmeasured[i][r]->SetMarkerSize(marker_size);
			hplot[i][r]=(TH1D*) hmeasured[i][r]->Rebin(32,Form("hplot_%i_%i",i,r),new_bins2);
			hplot[i][r]->Sumw2();
			hplot[i][r]->Scale(scale_jets[r],"width");
		}
		//f->Close();
		//delete f;
		
		TString str = Form("../../plotting_out/root/toymodel%s/%s/root_RAA%.1lf_%s/histos_jets_R%.1lf_pTcut0.2.root",cent_suff.Data(),version_toy.Data(),RAAtoy,version_toy.Data(),Rs[r]);
		cout<<str.Data()<<endl;
		ftoy[r] = new TFile(str.Data(), "OPEN");
		TH1I* hevents_toy;
		hevents_toy= (TH1I*) ftoy[r]->Get("hevents");
		int nevents_toy=hevents_toy->GetEntries();
		if(r==0) mevents_toy=(float)nevents_toy/1000000;
		scale_jets_toy_rec[r] = 1./(2*(1.-Rs[r])*2.*TMath::Pi()*nevents_toy);
		scale_jets_toy_gen[r] = 1./(2*2.*TMath::Pi()*nevents_toy);
		
		TH2D* histo2d=(TH2D*) ftoy[r]->Get("fhDSpTleading");
		for(int i=0; i<nhistos_toy; i++)
		{
			hmeasured_toy[i][r]=histo2d->ProjectionY(Form("histo1d_toy_%i_%i",i,r),histo2d->GetXaxis()->FindBin(pTlead[i]),histo2d->GetXaxis()->GetNbins());
			
			//hmeasured[i][r]->SetMarkerStyle(marker_list[i+1]);
			hmeasured_toy[i][r]->SetMarkerStyle(31);
			hmeasured_toy[i][r]->SetLineColor(colorList[i+nhistos]);
			hmeasured_toy[i][r]->SetLineWidth(lineLineWidth);
			hmeasured_toy[i][r]->SetMarkerColor(colorList[i+nhistos]);
			hmeasured_toy[i][r]->SetMarkerSize(marker_size);
			hplot_toy[i][r]=(TH1D*) hmeasured_toy[i][r]->Rebin(32,Form("hplottoy_%i_%i",i,r),new_bins2);
			hplot_toy[i][r]->Sumw2();
			hplot_toy[i][r]->Scale(scale_jets_toy_rec[r],"width");
			/* cout<<"integral1: "<<hmeasured_toy[i][r]->Integral("width")<<endl;
			*  hmeasured_toy[i][r]->Scale(hmeasured[i][r]->Integral("width")/hmeasured_toy[i][r]->Integral("width"));
			*  cout<<"integral2: "<<hmeasured_toy[i][r]->Integral("width")<<endl;
			*  cout<<"integral3: "<<hmeasured[i][r]->Integral("width")<<endl;
			*/
		}
		//hard jet only distribution
		if(showHardJet)
		{
			str = Form("../../plotting_out/root/toymodel%s/%s/root_RAA%.1lf_%s/jetonly_R%.1lf_pTcut0.2.root",cent_suff.Data(),version_toy.Data(),RAAtoy,version_toy.Data(),Rs[r]);
			cout<<str.Data()<<endl;
			fhtoy[r] = new TFile(str.Data(), "OPEN");
			TH1I* hevents_htoy;
			hevents_htoy= (TH1I*) fhtoy[r]->Get("hevents");
			int nevents_htoy=hevents_htoy->GetEntries();
			if(r==0) mevents_htoy=(float)nevents_htoy/1000000;
			scale_jets_toy_hard[r] = 1./(2*(1.-Rs[r])*2.*TMath::Pi()*nevents_htoy);
			TH2D* histo2dh=(TH2D*) fhtoy[r]->Get("fhPtRecpTleading");
		
			hmeasured_htoy[r]=(TH1D*) histo2dh->ProjectionY(Form("histo1d_toy_%i",i),histo2dh->GetXaxis()->FindBin(pTlead_toy),histo2dh->GetXaxis()->GetNbins());
			hmeasured_htoy[r]->SetMarkerStyle(21);
			hmeasured_htoy[r]->SetLineColor(colorList[i]+nhistos_toy+1);
			hmeasured_htoy[r]->SetLineWidth(graphLineWidth);
			hmeasured_htoy[r]->SetMarkerColor(colorList[i]+nhistos_toy+1);
			hmeasured_htoy[r]->SetMarkerSize(marker_size);
			hmeasured_htoy[r]->Rebin(10);
			hmeasured_htoy[r]->Scale(scale_jets_toy_hard[r],"width");
		}
		
		//ratio toymodel/STAR
		for (int i=0; i<nhistos; i++)
		{
			htoydatarat[i][r]=(TH1D*)hplot_toy[i][r]->Clone(Form("htoydatarat_%i",i));
			htoydatarat[i][r]->Divide(hplot[i][r]);
			htoydatarat[i][r]->SetMarkerStyle(marker_list[i]);
			htoydatarat[i][r]->SetMarkerColor(colorList[i+nhistos]);
		}
		
	}//R loop
	
		//ftoy->Close();
		//delete ftoy;
		/*
		//particle level toymodel hard jet spectrum 
		TString str = Form("../../plotting_out/root/toymodel%s/%s/root_RAA%.1lf_%s/control_histos.root",cent_suff.Data(),version_toy.Data(),RAAtoy,version_toy.Data());
		ftoy2[r] = new TFile(str.Data(), "OPEN");
		TH2D* hgen2d=(TH2D*) ftoy2[r]->Get("hpTpTleadGen");
		TH1D* hgen=(TH1D*)hgen2d->ProjectionY("hgen",hgen2d->GetXaxis()->FindBin(pTlead_toy),hgen2d->GetXaxis()->GetNbins());
		//TH1D* hgen=(TH1D*) ftoy2[r]->Get("hjetgenpT");
		TH2D* hreq2d=(TH2D*) ftoy2[r]->Get("hpTpTleadReq");
		TH1D* hreq=(TH1D*)hreq2d->ProjectionY("hreq",hreq2d->GetXaxis()->FindBin(pTlead_toy),hreq2d->GetXaxis()->GetNbins());
		//TH1D* hreq=(TH1D*) ftoy2[r]->Get("hjetreqpT");
		hgen->Scale(scale_jets_toy_gen,"width");
		hreq->Scale(scale_jets_toy_gen,"width");
		hgen->SetMarkerStyle(31);
		hgen->SetLineWidth(graphLineWidth);
		hgen->SetMarkerSize(marker_size);
		hgen->SetLineColor(kBlue);
		hgen->SetMarkerColor(kBlue);
		hreq->SetMarkerStyle(32);
		hreq->SetLineWidth(graphLineWidth);
		hreq->SetMarkerSize(marker_size);
		hreq->SetLineColor(kGreen+2);
		hreq->SetMarkerColor(kGreen+2);
		*/
		
		double jet_norm=TAA;
		double B[3] = {2.41349e+01,1.75255e+01,1.27430e+01}; //R=0.2,0.3,0.4
		double T[3] = {2.69190e+00,2.53251e+00,2.40521e+00};
		double n[3] = {1.42495e+02,1.03738e+02,8.14043e+01};
		double m0[3] = {-3.00000e+00,-7.83741e+00,-8.40985e+00};
		double mu[3] = {-2.97634e+01,-2.01221e+01,-1.62948e+01};
		double A[3] = {1.12732e+00,5.45759e-01,3.25500e-01};
		double pwr[3] = {4.82216e+00,4.37727e+00,4.08808e+00};
				
		//int ridx=R*10-2; //0:R=0.2 1:R=0.3 2:R=0.4
		
		TH1 *frame = new TH1I("frame", "", 1000, -100, +100);
		/*
		//--------------------------------
		//STAR and Toymodel in two panels
		//--------------------------------
		
		TCanvas *spectra = new TCanvas("cspectra","cspectra",10,10,2*can_x,can_y);
		cspectra->Divide(2,1);
		cspectra->cd(1);
		gPad->SetLogy();
		gPad->SetGrid();
		
		frame->GetXaxis()->SetTitle("p_{T, raw}^{charged} (GeV/#it{c})");
		frame->GetYaxis()->SetTitle(ytitle);
		frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
		frame->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
		//frame->GetYaxis()->SetRangeUser(1E-3, 1E5);
		str="";
		frame->SetTitle(str);
		frame->DrawCopy("");
		for(int i=0; i<nhistos; i++){
			hplot[i][r]->DrawCopy("E same");
		}
		
		TLegend *model_info = new TLegend(0.20, 0.48, 0.35, 0.900);
		model_info->SetTextSize(0.04);
		model_info->SetFillStyle(0);
		model_info->SetBorderSize(0);
		model_info->SetMargin(0.05);
		model_info->SetHeader("STAR Au+Au #sqrt{s_{NN}} = 200 GeV");
		model_info->AddEntry("",G_system_info[system].Data(),"");
		//model_info->AddEntry("", Form("%.1lfM events",mevents), "");
		//model_info->AddEntry("", Form("Run11, %s",trigger.Data()), "");
		model_info->AddEntry("", "Charged jets", "");
		model_info->AddEntry("", Form("anti-k_{T}, #bf{R = %.1lf}",R), "");
		model_info->AddEntry("","p_{T}^{const}>0.2GeV/#it{c}","");
		model_info->DrawClone("same");
		
		TLatex *latex = new TLatex();
		latex->SetNDC();
		latex->SetTextSize(0.035);
		//latex->DrawLatex(0.4, 0.3,"STAR Preliminary");
		
		TLegend *legspectra = new TLegend(0.6, 0.60, 0.85, 0.90);
		legspectra->SetTextSize(0.04);
		legspectra->SetFillStyle(0);
		legspectra->SetBorderSize(0);
		for(int i=0; i<nhistos; i++){
			legspectra->AddEntry(hplot[i][r], Form("p_{T}^{lead}>%.0lf GeV/#it{c}",pTlead[i]), "lp");
		}
		
		legspectra->DrawClone("same");
		
		cspectra->cd(2);
		gPad->SetLogy();
		gPad->SetGrid();
		
		frame->GetXaxis()->SetTitle("p_{T, raw}^{charged} (GeV/#it{c})");
		frame->GetYaxis()->SetTitle(ytitle);
		frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
		frame->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
		//frame->GetYaxis()->SetRangeUser(1E-3, 1E5);
		str="";
		frame->SetTitle(str);
		frame->DrawCopy("");
		hplot[icompare]->DrawCopy("E same");
		//for(int i=0; i<nhistos_toy; i++){
		hplot_toy[icompare]->DrawCopy("E same");
		//}
		if(showHardJet)hmeasured_htoy->DrawCopy("E same");
		//hreq->DrawCopy("E same");
		//hgen->DrawCopy("E same");
		//fjettrue->DrawCopy("same");
		
		
		
		TLegend *model_info2 = new TLegend(0.14, 0.48, 0.35, 0.900);
		model_info2->SetTextSize(0.035);
		model_info2->SetFillStyle(0);
		model_info2->SetBorderSize(0);
		model_info2->SetMargin(0.05);
		model_info2->SetHeader("STAR vs. TOYMODEL");
		//model_info2->AddEntry("", "0-10% Central Collisions", "");
		model_info2->AddEntry("", Form("anti-k_{T}, R = %.1lf",R), "");
		model_info2->AddEntry("", Form("p_{T,lead}^{min} = %.0lf GeV/#it{c}",pTlead_toy), "");
		model_info2->AddEntry("", "STAR Data:", "");
		model_info2->AddEntry("", Form("   %.1lfM events",mevents), "");
		model_info2->AddEntry("", "TOYMODEL:", "");
		model_info2->AddEntry("", Form("   %.1lfM events",mevents_toy), "");

		model_info2->DrawClone("same");
		
		TLegend *legspectra = new TLegend(0.6, 0.70, 0.89, 0.85);
		legspectra->SetTextSize(0.04);
		legspectra->SetFillStyle(0);
		legspectra->SetBorderSize(0);
		legspectra->AddEntry(hplot[icompare], Form("STAR Data"), "lp");
		//for(int i=0; i<nhistos_toy; i++){
		legspectra->AddEntry(hplot_toy[icompare], Form("TOYMODEL"), "lp");
		//}
		if(showHardJet)legspectra->AddEntry(hmeasured_htoy, Form("TOYMODEL - hard jets"), "lp");
		//legspectra->AddEntry(fjettrue, Form("TM hard jet generating function"), "lp");
		//legspectra->AddEntry(hreq, "TM parton level spec.", "lp");
		//legspectra->AddEntry(hgen, "TM detector level spec.", "lp");
		legspectra->DrawClone("same");
		
		str = Form("%s/spectra_pTreco_STAR_vs_toymodel%s_pTlead%.0lf_R0%.0lf.%s", outdir.Data(),cent_suff.Data(),pTlead_toy,R*10,ext.Data());
		cspectra->SaveAs(str.Data());*/
		
		//--------------------------------
		//STAR and Toymodel in one panel
		//--------------------------------
		
		//pTleading dependence + toymodel in one plot
		TCanvas *cspectradt = new TCanvas("cspectradt","cspectradt",10,10,2*can_x,1.2*can_y);
		double eps=0.02;
		TPad* p1[2];
		TPad* p2[2];
		
		TLegend *model_info = new TLegend(0.20, 0.48, 0.35, 0.900);
		model_info->SetTextSize(0.04);
		model_info->SetFillStyle(0);
		model_info->SetBorderSize(0);
		model_info->SetMargin(0.05);
		model_info->SetHeader("Au+Au #sqrt{s_{NN}} = 200 GeV");
		
		//coordinates of the text label (R=...)
		double labx[]={0.2,0.2,0.2};
		double laby[]={0.5,0.8,0.8};
		
		//x-coordinates of the pads
		double pxmin[]={0,0.33,0.66};
		double pxmax[]={0.33-eps,0.66-eps,1};
		if(nR==2)
		{
			pxmin[1]=0.52;
			pxmax[0]=0.52;
			pxmax[1]=1;
		}
		
		for(int r=0;r<nR;r++)
		{
			//top pad - spectra
			TString pname=Form("p%i_top",r);
			p1[r] = new TPad(pname,pname,pxmin[r],0.35-eps,pxmax[r],1,0); p1[r]->Draw();
			p1[r]->SetBottomMargin(eps);
			if(r<nR-1) p1[r]->SetRightMargin(0);
			if(r>0) p1[r]->SetLeftMargin(0);
			if(showGrid)p1[0]->SetGrid();
			p1[r]->SetLogy();
			
			//bottom pad - ratios
			pname=Form("p%i_bottom",r);
			p2[r] = new TPad(pname,pname,pxmin[r],0,pxmax[r],0.35*(1.-eps),0); p2[r]->Draw(); 
			p2[r]->SetTopMargin(0);
			p2[r]->SetBottomMargin(0.3);
			if(r<nR-1) p2[r]->SetRightMargin(0);
			if(r>0) p2[r]->SetLeftMargin(0);
			if(showGrid)p2[r]->SetGrid();
			p2[r]->SetFillColor(0);
			p2[r]->SetFillStyle(0);
			//p2[r]->SetLogy();
		}
		
		for(int r=0;r<nR;r++)
		{
			p1[r]->cd();
		
			frame->GetYaxis()->SetTitle(ytitle);
			frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
			frame->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
			//frame->GetYaxis()->SetRangeUser(1E-3, 1E5);
			str="";
			frame->SetTitle(str);
			frame->DrawCopy("");
			for(int i=0; i<nhistos; i++){
				hplot[i][r]->DrawCopy("E same");
			}
			for(int i=0; i<nhistos_toy; i++){
				hplot_toy[i][r]->DrawCopy("c hist same");
			}
			
			if(r==0)
			{
				model_info->AddEntry("",G_system_info[system].Data(),"");
				//model_info->AddEntry("", Form("%.1lfM events",mevents), "");
				//model_info->AddEntry("", Form("Run11, %s",trigger.Data()), "");
				model_info->AddEntry("", "Charged jets", "");
				model_info->AddEntry("", Form("anti-k_{T}, R = %.1lf",Rs[r]), "");
				//model_info->AddEntry("", "anti-k_{T}", "");
				model_info->AddEntry("","p_{T}^{const} > 0.2 GeV/#it{c}","");
			}
			if(r==0) model_info->DrawClone("same");  
			
			TLatex *latex = new TLatex();
			latex->SetNDC();
			latex->SetTextSize(0.045);
			if(r>0)latex->DrawLatex(labx[r], laby[r],Form("R = %.1lf",Rs[r]));
			
			TLegend *legspectra = new TLegend(0.65, 0.5, 0.9, 0.90);
			legspectra->SetTextSize(0.045);
			legspectra->SetFillStyle(0);
			legspectra->SetBorderSize(0);
			if(r==0)
			{
				legspectra->SetHeader("");
				latex->DrawLatex(0.65, 0.85,"STAR, p_{T,lead}^{min}:");
				for(int i=0; i<nhistos; i++)
				{
					legspectra->AddEntry(hplot[i][r], Form("%.0lf GeV/#it{c}",pTlead[i]), "p");
				}
			}
			else if(r==1)
			{
				legspectra->SetHeader("");
				latex->DrawLatex(0.55, 0.85,"PM, p_{T,lead}^{min}:");
				for(int i=0; i<nhistos_toy; i++)
				{
					legspectra->AddEntry(hplot_toy[i][r],Form("%.0lf GeV/#it{c}",pTlead[i]),"l");
				}
			}
			
			legspectra->DrawClone("same");
			
			p2[r]->cd();
			p2[r]->SetLogy();
			TString xttl="";
			if(r==nR-1) xttl="p_{T, jet}^{reco, ch} (GeV/#it{c})";
			frame->GetXaxis()->SetTitle(xttl);
			if(nR==3) frame->GetXaxis()->CenterTitle();
			frame->GetYaxis()->SetTitle("PM/STAR    ");
			frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
			frame->GetYaxis()->SetRangeUser(0.1, 9);
			
			frame->GetYaxis()->SetLabelSize(0.08);
			frame->GetXaxis()->SetLabelSize(0.08);
			frame->GetXaxis()->SetTitleSize(0.12);
			frame->GetYaxis()->SetTitleSize(0.09);
			frame->GetYaxis()->SetTitleOffset(0.5);
			
			
			frame->SetTitle("");
			frame->DrawCopy("");
			
			for(int i=0;i<nhistos;i++)
			{
				if(compare[i]) htoydatarat[i][r]->Draw("E same");
			}
			
			TLegend *legtoy = new TLegend(0.20,0.7,0.65,0.95);
			legtoy->SetTextSize(0.085);
			legtoy->SetFillStyle(0);
			legtoy->SetBorderSize(0);
			legtoy->SetHeader("p_{T,lead}^{min}:");
			for(int i=0;i<nhistos;i++)
			{
				if(compare[i]) legtoy->AddEntry(htoydatarat[i][r],Form("%.0lf GeV/#it{c}",pTlead[i]),"lp");
			}
			if(r==0) legtoy->DrawClone("same");
			
			
			TLine *one = new TLine(spectraXmin, 1, spectraXmax, 1);
			one->SetLineWidth(lineLineWidth);
			one->SetLineStyle(2);
			one->SetLineColor(kBlack);
			one->DrawClone("same");
		}//R loop
		str = Form("%s/spectra1p_pTreco_STAR_vs_toymodel%s_pTlead%.0lf.%s", outdir.Data(),cent_suff.Data(),pTlead_toy,ext.Data());
		cspectradt->SaveAs(str.Data());
		
		
}

//========================================================================
//STAR vs Toymodel (OLD)
//========================================================================
/*
plot_measured_vs_toymodel(float R=0.3,float pTlead=5.0, systems system=cent, TString version="GPC2", TString version_toy="newBinning2", TString ext = "gif")
{ 
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	gStyle->SetOptDate(0);
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
	
	
	bool showHardJet=0; //show also reconstructed hard jet distribution
	
	TString outdir = "../../plotting_out/obr/intersteps/toy_vs_star"; //figure output directory
	
	double RAAtoy=0.5;
	if(system==peri) RAAtoy=0.7;
	
	const float marker_size=1.0;
	const float line_width=2;
	
	const Color_t colorList[]={kBlack,kRed,kGreen+3,kMagenta+2,kBlue,kOrange+2,kYellow+2,kBlue-4,kOrange-2,kGreen,kRed+2,7,kGreen+4,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,28};
	const int marker_list[]={21,20,22,30};
	
	//canvas size
	const int can_x=1200; //1600
	const int can_y=680; //900
	
	TString trigger="MB";
	TString cent_suff="";
	if(system==peri)
	{
		trigger="MB_peripheral";
		cent_suff="_peripheral";
	}
	const  float spectraXmin=-30;
	const float spectraXmax=50;
	const float spectraYmin=1E-8;
	const float spectraYmax=1E1;
	
	const int nhistos=1;
	//const float pTlead[nhistos]={0, 4.,5.,7.};
	TH1D* hmeasured[nhistos];
	
	//int mevents=0;
	
	double TAA=22.75;
	if(system==peri)TAA=0.5;
	
	TString str = Form("../../plotting_out/root/%s/%s/histos_inclusivejet_normal.root",trigger.Data(),version.Data());
	TFile *f = new TFile(str.Data(), "OPEN");
	TH1I* hevents;
	hevents= (TH1I*) f->Get("hevents");
	int nevents=hevents->GetBinContent(2);
	float mevents=(float)nevents/1000000;
	double hole=0;//1/(1.0/12.)+((R*2.)/(2.*TMath::Pi()));
	double scale_jets = 1./(2*(1.-R)*2.*TMath::Pi()*nevents*(1-hole));
	//TH2D* histo2d=(TH2D*) f->Get("hpT_pTlead_nobadsecedge1R");
	TH2D* histo2d=(TH2D*) f->Get(Form("hpT_pTlead_R0%.0lf",R*10));
	
	for(int i=0; i<nhistos; i++){
		hmeasured[i]=histo2d->ProjectionX(Form("histo1d_%i",i),histo2d->GetYaxis()->FindBin(pTlead),histo2d->GetYaxis()->GetNbins());
		//hmeasured[i]->SetMarkerStyle(marker_list[i+1]);
		hmeasured[i]->SetMarkerStyle(30);
		hmeasured[i]->SetLineColor(colorList[i]);
		hmeasured[i]->SetLineWidth(line_width);
		hmeasured[i]->SetMarkerColor(colorList[i]);
		hmeasured[i]->SetMarkerSize(marker_size);
		hmeasured[i]->Scale(scale_jets,"width");
	}
	//f->Close();
	//delete f;
	
	TString str = Form("../../plotting_out/root/toymodel%s/%s/root_RAA%.1lf_%s/histos_jets_R%.1lf_pTcut0.2.root",cent_suff.Data(),version_toy.Data(),RAAtoy,version_toy.Data(),R);
	cout<<str.Data()<<endl;
	TFile *ftoy = new TFile(str.Data(), "OPEN");
	TH1I* hevents_toy;
	hevents_toy= (TH1I*) ftoy->Get("hevents");
	int nevents_toy=hevents_toy->GetEntries();
	float mevents_toy=(float)nevents_toy/1000000;
	double scale_jets_toy_rec = 1./(2*(1.-R)*2.*TMath::Pi()*nevents_toy);
	double scale_jets_toy_gen = 1./(2*2.*TMath::Pi()*nevents_toy);
	
	TH2D* histo2d=(TH2D*) ftoy->Get("fhDSpTleading");
	
	for(int i=0; i<nhistos; i++){
		hmeasured_toy[i][r]=histo2d->ProjectionY(Form("histo1d_toy_%i",i),histo2d->GetXaxis()->FindBin(pTlead),histo2d->GetXaxis()->GetNbins());
		
		//hmeasured[i]->SetMarkerStyle(marker_list[i+1]);
		hmeasured_toy[i][r]->SetMarkerStyle(31);
		hmeasured_toy[i][r]->SetLineColor(colorList[i]+nhistos);
		hmeasured_toy[i][r]->SetLineWidth(line_width);
		hmeasured_toy[i][r]->SetMarkerColor(colorList[i]+nhistos);
		hmeasured_toy[i][r]->SetMarkerSize(marker_size);
		hmeasured_toy[i][r]->Rebin(10);
		hmeasured_toy[i][r]->Scale(scale_jets_toy_rec,"width");
	}
	
	//hard jet only distribution
	if(showHardJet){
		str = Form("../../plotting_out/root/toymodel%s/%s/root_RAA%.1lf_%s/jetonly.root",cent_suff.Data(),version_toy.Data(),RAAtoy,version_toy.Data());
		cout<<str.Data()<<endl;
		TFile *fhtoy = new TFile(str.Data(), "OPEN");
		TH1I* hevents_htoy;
		hevents_htoy= (TH1I*) fhtoy->Get("hevents");
		int nevents_htoy=hevents_htoy->GetEntries();
		float mevents_htoy=(float)nevents_htoy/1000000;
		double scale_jets_toy_hard = 1./(2*(1.-R)*2.*TMath::Pi()*nevents_htoy);
		TH2D* histo2dh=(TH2D*) fhtoy->Get("fhPtRecpTleading");}
		
		TH1D* hmeasured_htoy;
		if(showHardJet){
			hmeasured_htoy=(TH1D*) histo2dh->ProjectionY(Form("histo1d_toy_%i",i),histo2dh->GetXaxis()->FindBin(pTlead),histo2dh->GetXaxis()->GetNbins());
			hmeasured_htoy->SetMarkerStyle(21);
			hmeasured_htoy->SetLineColor(colorList[i]+nhistos+1);
			hmeasured_htoy->SetLineWidth(line_width);
			hmeasured_htoy->SetMarkerColor(colorList[i]+nhistos+1);
			hmeasured_htoy->SetMarkerSize(marker_size);
			hmeasured_htoy->Rebin(10);
			hmeasured_htoy->Scale(scale_jets_toy_hard,"width");
		}
		
		
		//ftoy->Close();
		//delete ftoy;
		
		//particle level toymodel hard jet spectrum 
		TString str = Form("../../plotting_out/root/toymodel%s/%s/root_RAA%.1lf_%s/control_histos.root",cent_suff.Data(),version_toy.Data(),RAAtoy,version_toy.Data());
		TFile *ftoy2 = new TFile(str.Data(), "OPEN");
		TH2D* hgen2d=(TH2D*) ftoy2->Get("hpTpTleadGen");
		TH1D* hgen=(TH1D*)hgen2d->ProjectionY("hgen",hgen2d->GetXaxis()->FindBin(pTlead),hgen2d->GetXaxis()->GetNbins());
		//TH1D* hgen=(TH1D*) ftoy2->Get("hjetgenpT");
		TH2D* hreq2d=(TH2D*) ftoy2->Get("hpTpTleadReq");
		TH1D* hreq=(TH1D*)hreq2d->ProjectionY("hreq",hreq2d->GetXaxis()->FindBin(pTlead),hreq2d->GetXaxis()->GetNbins());
		//TH1D* hreq=(TH1D*) ftoy2->Get("hjetreqpT");
		hgen->Scale(scale_jets_toy_gen,"width");
		hreq->Scale(scale_jets_toy_gen,"width");
		hgen->SetMarkerStyle(31);
		hgen->SetLineWidth(line_width);
		hgen->SetMarkerSize(marker_size);
		hgen->SetLineColor(kBlue);
		hgen->SetMarkerColor(kBlue);
		hreq->SetMarkerStyle(32);
		hreq->SetLineWidth(line_width);
		hreq->SetMarkerSize(marker_size);
		hreq->SetLineColor(kGreen+2);
		hreq->SetMarkerColor(kGreen+2);
		
		
		double jet_norm=TAA;
		double B[3] = {2.41349e+01,1.75255e+01,1.27430e+01}; //R=0.2,0.3,0.4
		double T[3] = {2.69190e+00,2.53251e+00,2.40521e+00};
		double n[3] = {1.42495e+02,1.03738e+02,8.14043e+01};
		double m0[3] = {-3.00000e+00,-7.83741e+00,-8.40985e+00};
		double mu[3] = {-2.97634e+01,-2.01221e+01,-1.62948e+01};
		double A[3] = {1.12732e+00,5.45759e-01,3.25500e-01};
		double pwr[3] = {4.82216e+00,4.37727e+00,4.08808e+00};
		
		int ridx=R*10-2; //0:R=0.2 1:R=0.3 2:R=0.4
		
		TF1 *fjettrue=new TF1("fjettrue",hardjet,3,40,9);
		fjettrue->SetParameters(jet_norm,B[ridx],T[ridx],n[ridx],m0[ridx],mu[ridx],A[ridx],pwr[ridx],RAAtoy);
		fjettrue->SetNpx(10000);
		fjettrue->SetLineColor(kMagenta);
		
		TH1 *frame = new TH1I("frame", "", 1000, -100, +100);
		
		TCanvas *cspectra = new TCanvas("cspectra","cspectra",10,10,can_x,can_y);
		cspectra->cd();
		cspectra->SetGrid();
		cspectra->SetLogy();
		frame->GetXaxis()->SetTitle("p_{T, reco}^{charged} (GeV/#it{c})");
		frame->GetYaxis()->SetTitle("1/N_{events} 1/2#pi d^{2}N/dp_{T}^{ch}d#eta (GeV/#it{c})^{-1}");
		frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
		frame->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
		//frame->GetYaxis()->SetRangeUser(1E-3, 1E5);
		str="";
		frame->SetTitle(str);
		frame->DrawCopy("");
		for(int i=0; i<nhistos; i++){
			hmeasured[i]->DrawCopy("E same");
			hmeasured_toy[i][r]->DrawCopy("E same");
		}
		if(showHardJet)hmeasured_htoy->DrawCopy("E same");
		hreq->DrawCopy("E same");
		hgen->DrawCopy("E same");
		fjettrue->DrawCopy("same");
		
		
		
		TLegend *model_info = new TLegend(0.20, 0.48, 0.35, 0.900);
		model_info->SetTextSize(0.035);
		model_info->SetFillStyle(0);
		model_info->SetBorderSize(0);
		model_info->SetMargin(0.05);
		model_info->SetHeader("Run 11 Au+Au #sqrt{s_{NN}} = 200 GeV");
		model_info->AddEntry("", "0-10% Central Collisions", "");
		model_info->AddEntry("", Form("anti-k_{T}, R = %.1lf",R), "");
		model_info->AddEntry("", Form("p_{T,lead}^{min} = %.0lf GeV/#it{c}",pTlead), "");
		model_info->AddEntry("", Form("Data: %.1lfM events",mevents), "");
		model_info->AddEntry("", Form("Parametrized model: %.1lfM events",mevents_toy), "");

		model_info->DrawClone("same");
		
		TLatex *latex = new TLatex();
		latex->SetNDC();
		latex->SetTextSize(0.035);
		//latex->DrawLatex(0.4, 0.3,"STAR Preliminary");
		
		TLegend *legspectra = new TLegend(0.6467, 0.60, 0.89, 0.90);
		legspectra->SetTextSize(0.04);
		legspectra->SetFillStyle(0);
		legspectra->SetBorderSize(0);
		for(int i=0; i<nhistos; i++){
			legspectra->AddEntry(hmeasured[i], Form("STAR Data"), "lp");
			legspectra->AddEntry(hmeasured_toy[i][r], Form("TOYMODEL - reconstructed"), "lp");
		}
		if(showHardJet)legspectra->AddEntry(hmeasured_htoy, Form("TOYMODEL-reco. hard jets"), "lp");
		legspectra->AddEntry(fjettrue, Form("TM hard jet generating function"), "lp");
		legspectra->AddEntry(hreq, "TM parton level spec.", "lp");
		legspectra->AddEntry(hgen, "TM detector level spec.", "lp");
		legspectra->DrawClone("same");
		
		str = Form("%s/spectra_pTreco_STAR_vs_toymodel%s_pTlead%.0lf_R%.1lf.%s", outdir.Data(),cent_suff.Data(),pTlead,R,ext.Data());
		cspectra->SaveAs(str.Data());
		
		TH1D* hratio[nhistos];
		for(int i=0; i<nhistos; i++){
			TString str=Form("hratio%i",i);
			hratio[i]=(TH1D*)hmeasured[i]->Clone(str);
			hratio[i]->Rebin(2);
			hratio[i]->Scale(1./2.);
			hratio[i]->Divide(hmeasured_toy[i][r]);
		}
		
		TCanvas *cratio = new TCanvas("cratio","cratio",10,10,can_x,can_y);
		cratio->cd();
		cratio->SetGrid();
		for(int i=0; i<nhistos; i++){
			if(i==0){
				hratio[0]->Draw("e");
				hratio[0]->GetYaxis()->SetRangeUser(0, 2);
				hratio[0]->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
				hratio[0]->SetLineColor(kBlue);
				hratio[0]->SetMarkerColor(kBlue);
			}
			else hratio[i]->Draw("esame");
		}
		TLine *one = new TLine(spectraXmin, 1, spectraXmax, 1);
		one->SetLineWidth(2);
		one->SetLineStyle(2);
		one->SetLineColor(kBlack);
		one->DrawClone("same");
		
		str = Form("%s/ratio_pTreco_STAR_vs_toymodel%s_pTlead%.0lf_R%.1lf.%s", outdir.Data(),cent_suff.Data(),pTlead,R,ext.Data());
		cratio->SaveAs(str.Data());
}//STAR vs TOYMODEL
*/


//========================================================================
//STAR vs PYTHIA
//========================================================================

plot_measured_vs_pythia(float R=0.3, float pTlead=0, systems system=cent, TString version="GPC2", bool showSTAR=0,TString ext="pdf")
{
	/* gStyle->SetOptStat(0);
	 *  gStyle->SetPalette(1);
	 *  gStyle->SetOptDate(0);
	 *  gStyle->SetPadLeftMargin(0.12);
	 *  gStyle->SetPadRightMargin(0.09);
	 *  gStyle->SetPadTopMargin(0.1);
	 *  gStyle->SetPadBottomMargin(0.15);
	 *  gStyle->SetTitleSize(0.055,"Y");
	 *  gStyle->SetTitleOffset(0.95,"Y");
	 *  gStyle->SetTitleSize(0.06,"X");
	 *  gStyle->SetTitleOffset(0.95,"X");
	 *  gStyle->SetLabelSize(0.03,"X");
	 *  gStyle->SetLabelSize(0.03,"Y");*/
	
	double PtEmb=20; //pT of embedded probe used to delta pT calculation
	double pThard=0.0; //use pythia spectrum above this value
	
	TString outdir = Form("../../plotting_out/obr/%s/effi",version.Data()); //figure output directory
	
	const float marker_size=1.0;
	const float line_width=2;
	
	const Color_t colorList[]={kGray+1,kRed,kBlue,kBlack,kGreen+2,kMagenta,kOrange,kYellow,kBlue-4,kOrange-2,kGreen,kRed+2,7,kGreen+4,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,28};
	const int marker_list[]={21,20,22,30};
	
	const int lineStList[]={1,1,1,1,1,1,1};
	
	//canvas size
	const int can_x=1200; //1600
	const int can_y=680; //900
	
	TString trigger="MB";
	TString cent_suff="_central";
	float TAA=22.5;
	if(system==peri)
	{
		TAA=0.5;
		cent_suff="_peripheral";
	}
	else if(system==pp) 
	{
		cent_suff+="_pp";
		TAA=22.5;
	}
	trigger+=cent_suff;
	
	
	const float ppXsection=42;
	
	const  float spectraXmin=-25;
	const float spectraXmax=60;
	const float spectraYmin=1E-9;
	const float spectraYmax=1E1;
	
	
	TH1D* hmeasured;
	TH1D* hmeasured_tmp;
	int mevents=0;
	
	
	TString str = Form("../../plotting_out/root/%s/%s/histos_inclusivejet_normal.root",trigger.Data(),version.Data());
	TFile *f = new TFile(str.Data(), "OPEN");
	TH1I* hevents;
	hevents= (TH1I*) f->Get("hevents");
	int nevents=hevents->GetBinContent(2);
	mevents=(int)nevents/1000000;
	//double hole=(1.0/12.)+((R*2.)/(2.*TMath::Pi()));
	double hole=0;
	double scale_jets = 1./(2*(1.-R)*2.*TMath::Pi()*nevents*(1-hole));
	//TH2D* histo2d=(TH2D*) f->Get("hpT_pTlead");
	TH2D* histo2d=(TH2D*) f->Get(Form("hpT_pTlead_R0%.0lf",R*10));
	hmeasured_tmp=histo2d->ProjectionX("histo1d",histo2d->GetYaxis()->FindBin(pTlead),histo2d->GetYaxis()->GetNbins());
	//hmeasured->SetMarkerStyle(30);
	hmeasured_tmp->SetLineColor(G_colorFill[0]);
	hmeasured_tmp->SetLineWidth(line_width);
	//hmeasured->SetMarkerColor(colorList[0]);
	//hmeasured->SetMarkerSize(marker_size);
	hmeasured_tmp->SetFillStyle(3004);
	hmeasured_tmp->SetFillColor(G_colorFill[0]);
	hmeasured_tmp->Sumw2();
	TH1D* hmeasured=(TH1D*) hmeasured_tmp->Rebin(47,"hmeasured_rbn",new_bins1);
	hmeasured->Scale(scale_jets,"width");
	
	//TString name=Form("../../plotting_out/root/pp_pythia8/histos_pythiajet_R%.1lf%s_pTsmear_effi_dpTsmear.root",R,cent_suff.Data());
	TString name=Form("../../plotting_out/root/pp_pythia8/histos_pythiajet_R%.1lf%s_pTsmear_effi.root",R,cent_suff.Data());
	TFile *infile = new TFile(name.Data(),"OPEN");
	
	//PYTHIA spectrum - particle level
	TH1D* hptl_part=(TH1D*) infile->Get(Form("hpT_pTl%.0lf",pTlead));
	hptl_part->SetLineColor(colorList[1]);
	hptl_part->SetMarkerColor(colorList[1]);
	hptl_part->SetLineStyle(lineStList[1]);
	hptl_part->SetMarkerStyle(29);
	hptl_part->SetLineWidth(line_width);
	hptl_part->SetMarkerSize(0.9);
	
	//PYTHIA spectrum - detector level (momentum smearing, tracking efficiency, delta pT smearing)
	TH1D* htmp=(TH1D*) infile->Get(Form("hpT_pTl%.0lf_dete",pTlead));
	TH1D* hptl_dete111=(TH1D*) hmeasured->Clone("hptl_dete111");
	hptl_dete111->SetLineColor(colorList[2]);
	hptl_dete111->SetMarkerColor(colorList[2]);
	hptl_dete111->SetLineStyle(lineStList[2]);
	hptl_dete111->SetMarkerStyle(29);
	hptl_dete111->SetLineWidth(line_width);
	hptl_dete111->SetMarkerSize(0.9);
	hptl_dete111->SetFillStyle(0);
	hptl_dete111->Reset("MICE");
	
	//delta pT histogram for dpT smearing
	name = Form("../../plotting_out/intersteps/embedding/%s/%s/histos_embeddedjet_sp.root",trigger.Data(),version.Data());
	fdpt=new TFile(name.Data(), "OPEN");
	name=Form("delta_pt_BG_sp_%0.lf_R0%.0lf",pTlead,R*10);
	TH2D* h2dpt=(TH2D*)fdpt->Get(name);
	int bin=h2dpt->GetXaxis()->FindBin(PtEmb);
	TH1D* hdpt=(TH1D*)h2dpt->ProjectionY("hdpt",bin,bin);
	//hdpt->Draw();
	htmp->DrawClone();
	//smear pythia with delta-pT
	for(int bn=1; bn<=htmp->GetNbinsX();bn++)
	{
		double yval=htmp->GetBinContent(bn);
		//cout<<"y="<<yval<<endl;
		double stat=10000;
		double weight=yval/stat;
		//if(yval<100) continue;
		double ptb=htmp->GetBinCenter(bn);
		if(ptb<pThard)continue;
		//cout<<"pT:"<<ptb<<endl;
		for (int i=0; i<stat; i++)
		{
			
			double dpt=hdpt->GetRandom();
			if(dpt>13 || dpt<-10)continue; //these are dpT values with too large statistical errors
			hptl_dete111->Fill(ptb+dpt,weight);
			//if(ptb+dpt>40) cout<<"dpt:"<<dpt<<" pT:"<<ptb<<" weight:"<<yval/stat<<endl;
		}
	}
	delete htmp;
	
	
	//PYTHIA spectrum - detector level (momentum smearing)
	name=Form("../../plotting_out/root/pp_pythia8/histos_pythiajet_R%.1lf_pTsmear.root",R);
	TFile *infile2 = new TFile(name.Data(),"OPEN");
	TH1D* hptl_dete100=(TH1D*) infile2->Get(Form("hpT_pTl%.0lf_dete",pTlead));
	hptl_dete100->SetLineColor(colorList[3]);
	hptl_dete100->SetMarkerColor(colorList[3]);
	hptl_dete100->SetLineStyle(lineStList[3]);
	hptl_dete100->SetMarkerStyle(29);
	hptl_dete100->SetLineWidth(line_width);
	hptl_dete100->SetMarkerSize(0.9);
	
	//PYTHIA spectrum - detector level (tracking efficiency)
	name=Form("../../plotting_out/root/pp_pythia8/histos_pythiajet_R%.1lf%s_effi.root",R,cent_suff.Data());
	TFile *infile3 = new TFile(name.Data(),"OPEN");
	TH1D* hptl_dete010=(TH1D*) infile3->Get(Form("hpT_pTl%.0lf_dete",pTlead));
	hptl_dete010->SetLineColor(colorList[4]);
	hptl_dete010->SetMarkerColor(colorList[4]);
	hptl_dete010->SetLineStyle(lineStList[4]);
	hptl_dete010->SetMarkerStyle(29);
	hptl_dete010->SetLineWidth(line_width);
	hptl_dete010->SetMarkerSize(0.9);
	
	//PYTHIA spectrum - detector level (momentum smearing, tracking efficiency)
	name=Form("../../plotting_out/root/pp_pythia8/histos_pythiajet_R%.1lf%s_pTsmear_effi.root",R,cent_suff.Data());
	TFile *infile4 = new TFile(name.Data(),"OPEN");
	TH1D* hptl_dete110=(TH1D*) infile4->Get(Form("hpT_pTl%.0lf_dete",pTlead));
	hptl_dete110->SetLineColor(colorList[5]);
	hptl_dete110->SetMarkerColor(colorList[5]);
	hptl_dete110->SetLineStyle(lineStList[5]);
	hptl_dete110->SetMarkerStyle(29);
	hptl_dete110->SetLineWidth(line_width);
	hptl_dete110->SetMarkerSize(0.9);
	
	
	TH1I* heventsp=(TH1I*) infile->Get("hevts");
	int neventsp=heventsp->GetBinContent(2);
	
	TH1I* heventsp2=(TH1I*) infile2->Get("hevts");
	int neventsp2=heventsp2->GetBinContent(2);
	TH1I* heventsp3=(TH1I*) infile3->Get("hevts");
	int neventsp3=heventsp3->GetBinContent(2);
	TH1I* heventsp4=(TH1I*) infile4->Get("hevts");
	int neventsp4=heventsp4->GetBinContent(2);
	
	float jetscale=TAA/(neventsp);
	float jetscale2=TAA/(neventsp2);
	float jetscale3=TAA/(neventsp3);
	float jetscale4=TAA/(neventsp4);
	
	//hptl_dete111->Rebin(20);
	hptl_dete111->Sumw2();
	hptl_dete111->Scale(jetscale,"width");
	hptl_dete100->Rebin(20);
	hptl_dete100->Scale(jetscale2,"width");
	hptl_dete010->Rebin(20);
	hptl_dete010->Scale(jetscale3,"width");
	hptl_dete110->Rebin(20);
	hptl_dete110->Scale(jetscale4,"width");
	hptl_part->Rebin(20);
	hptl_part->Scale(jetscale,"width");
	
	//make a special histo to hide the artifical line of the hstogram between zero and 1st value
	TH1D* hwhite=hptl_part->Clone("hwhite");
	hwhite->Reset("MICE");
	int firstbin=hptl_part->GetXaxis()->FindBin(pTlead+0.01);
	double firstYval=hptl_part->GetBinContent(firstbin);
	hwhite->SetBinContent(firstbin,firstYval);
	hwhite->SetLineColor(kWhite);
	hwhite->SetLineWidth(line_width+1);
	
	TH1D* hwhite2=hptl_dete100->Clone("hwhite2");
	hwhite2->Reset("MICE");
	int firstbin2=hptl_dete100->GetXaxis()->FindBin(pTlead+0.01);
	double firstYval2=hptl_dete100->GetBinContent(firstbin2);
	hwhite2->SetBinContent(firstbin2,firstYval2);
	hwhite2->SetLineColor(kWhite);
	hwhite2->SetLineWidth(line_width+1);
	
	TH1D* hwhite3=hptl_dete110->Clone("hwhite3");
	hwhite3->Reset("MICE");
	int firstbin3=hptl_dete110->GetXaxis()->FindBin(0);
	for(int i=0;i<6;i++)//hide first few bins
	{
		double yval=hptl_dete110->GetBinContent(firstbin3+i);
		hwhite3->SetBinContent(firstbin3+i,yval);
	}
	hwhite3->SetLineColor(kWhite);
	hwhite3->SetLineWidth(line_width+1);
	
	TH1 *frame = new TH1I("frame", "", 1000, -100, +100);
	TH1 *frame2 = new TH1I("frame2", "", 1000, -100, +100);
	
	if(!showSTAR)
	{
		TCanvas *cspectra = new TCanvas("cspectra","cspectra",10,10,can_x,can_y);
		cspectra->cd();
		cspectra->SetLogy();
		
		frame->GetXaxis()->SetTitle("p_{T, jet}^{charged} (GeV/#it{c})");
		frame->GetYaxis()->SetTitle("1/N_{events} 1/2#pi d^{2}N/dp_{T}^{ch}d#eta (GeV/#it{c})^{-1}");
		frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
		frame->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
		
		str="";
		frame->SetTitle(str);
		frame->DrawCopy("");
		hptl_part->DrawCopy("histo L same");
		//hptl_dete100->DrawCopy("E same");
		hptl_dete010->DrawCopy("histo L same");
		hwhite->DrawCopy("histo L same");
		hwhite2->DrawCopy("histo L same");
		hptl_dete110->DrawCopy("histo L same");
		hwhite3->DrawCopy("histo L same");
		hptl_dete111->DrawCopy("histo L same");
		
		TLegend *model_info = new TLegend(0.20, 0.5, 0.35, 0.900);
		model_info->SetTextSize(0.04);
		model_info->SetFillStyle(0);
		model_info->SetBorderSize(0);
		model_info->SetMargin(0.05);
		model_info->SetHeader("PYTHIA pp #otimes T_{AA}, #sqrt{s}=200 GeV");
		model_info->AddEntry("", "central Au+Au background and instr. response", "");
		model_info->AddEntry("", Form("anti-k_{T}, R = %.1lf",R), "");
		model_info->AddEntry("", Form("p_{T,lead}^{min} = %.0lf",pTlead),"");
		
		//model_info->AddEntry("", Form("%iM events",mevents), "");
		model_info->DrawClone("same");
		
		TLegend *legspectra = new TLegend(0.6, 0.60, 0.89, 0.90);
		legspectra->SetTextSize(0.04);
		legspectra->SetFillStyle(0);
		legspectra->SetBorderSize(0);
		legspectra->AddEntry(hptl_part, "particle level ", "l");
		//legspectra->AddEntry(hptl_dete100, "PYTHIA p_{T}-smear", "l");
		legspectra->AddEntry(hptl_dete010, "tracking eff.", "l");
		legspectra->AddEntry(hptl_dete110, "tr. eff. + p_{T}-smear", "l");
		legspectra->AddEntry(hptl_dete111, "tr. eff. + p_{T}-smear + #delta p_{T}", "l");
		legspectra->DrawClone("same");
		
		str=Form("%s/PYTHIA_EffSize_R0%.0lf_pTlead%.0lf%s.%s",outdir.Data(),10*R,pTlead,cent_suff.Data(),ext.Data());
		cspectra->SaveAs(str.Data());
	}
	
	else{ //show also  STAR data and ratio PYTHIA to STAR
		TCanvas *cspectra = new TCanvas("cspectra","cspectra",10,10,can_x,can_y*1.3);
		cspectra->cd();
		//cspectra->SetGrid();
		//cspectra->SetLogy();
		double eps=0.02;
		TPad* p1 = new TPad("p1","p1",0,0.35-eps,1,1,0); p1->Draw();
		p1->SetBottomMargin(eps);
		//p1->SetGrid();
		p1->SetLogy();
		
		TPad* p2 = new TPad("p2","p2",0,0,1,0.35*(1.-eps),0); p2->Draw(); 
		p2->SetTopMargin(0);
		p2->SetBottomMargin(0.25);
		//p2->SetGrid();
		p2->SetFillColor(0);
		p2->SetFillStyle(0);
		
		p1->cd();
		
		frame->GetXaxis()->SetTitle("p_{T, jet}^{charged} (GeV/#it{c})");
		frame->GetYaxis()->SetTitle("1/N_{events} 1/2#pi d^{2}N/dp_{T}^{ch}d#eta (GeV/#it{c})^{-1}");
		frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
		frame->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
		//frame->GetYaxis()->SetRangeUser(1E-3, 1E5);
		str="";
		frame->SetTitle(str);
		frame->DrawCopy("");
		hptl_part->DrawCopy("E same");
		hmeasured->DrawCopy("histo same");
		//hptl_dete100->DrawCopy("E same");
		hptl_dete010->DrawCopy("E same");
		hptl_dete110->DrawCopy("E same");
		hptl_dete111->DrawCopy("E same");
		
		TLegend *model_info = new TLegend(0.20, 0.48, 0.35, 0.900);
		model_info->SetTextSize(0.035);
		model_info->SetFillStyle(0);
		model_info->SetBorderSize(0);
		model_info->SetMargin(0.05);
		model_info->SetHeader("Au+Au #sqrt{s_{NN}} = 200 GeV");
		model_info->AddEntry("", "0-10% Central Collisions", "");
		model_info->AddEntry("", Form("anti-k_{T}, R = %.1lf",R), "");
		model_info->AddEntry("", Form("p_{T,lead}^{min} = %.0lf",pTlead),"");
		
		//model_info->AddEntry("", Form("%iM events",mevents), "");
		model_info->DrawClone("same");
		
		TLegend *legspectra = new TLegend(0.5, 0.60, 0.89, 0.90);
		legspectra->SetTextSize(0.04);
		legspectra->SetFillStyle(0);
		legspectra->SetBorderSize(0);
		legspectra->AddEntry(hmeasured, "STAR data (uncorrected)", "f");
		legspectra->AddEntry(hptl_part, "PYTHIA (a): particle level ", "lp");
		//legspectra->AddEntry(hptl_dete100, "PYTHIA p_{T}-smear", "lp");
		legspectra->AddEntry(hptl_dete010, "PYTHIA (b): tracking eff.", "lp");
		legspectra->AddEntry(hptl_dete110, "PYTHIA (c): tr. eff. + p_{T}-smear", "lp");
		legspectra->AddEntry(hptl_dete111, "PYTHIA (d): tr. eff. + p_{T}-smear + #delta p_{T}", "lp");
		legspectra->DrawClone("same");
		
		
		p2->cd();
		
		//crat->SetLogy();
		frame2->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
		frame2->GetYaxis()->SetRangeUser(0.0, 1.5);
		frame2->GetXaxis()->SetTitle("p_{T, jet}^{charged} (GeV/#it{c})");
		frame2->GetYaxis()->SetTitle("STAR/PYTHIA (d)");
		frame2->GetYaxis()->SetLabelSize(0.06);
		frame2->GetXaxis()->SetLabelSize(0.06);
		frame2->GetYaxis()->SetTitleSize(0.09);
		frame2->GetYaxis()->SetTitleOffset(0.4);
		frame2->GetXaxis()->SetTitleSize(0.10);
		frame2->DrawCopy("");
		
		TH1D* hrat=(TH1D*) hmeasured->Clone("hrat");
		hrat->Reset("MICE");
		for(int bin=1; bin< hrat->GetNbinsX();bin++)
		{
			double pt=hmeasured->GetBinCenter(bin);
			int bin2=hptl_dete111->FindBin(pt);
			double val=hptl_dete111->GetBinContent(bin2);
			if(val>0)  hrat->SetBinContent(bin,hmeasured->GetBinContent(bin)/val);
		}
		hrat->Draw("esame");
		
		TLine *one = new TLine(spectraXmin, 1, spectraXmax, 1);
		one->SetLineWidth(2);
		one->SetLineStyle(2);
		one->SetLineColor(kBlack);
		one->DrawClone("same");
		
		str=Form("%s/STAR_vs_PYTHIA_R0%.0lf_pTlead%.0lf%s.%s",outdir.Data(),10*R,pTlead,cent_suff.Data(),ext.Data());
		cspectra->SaveAs(str.Data());
	}
}

//--------------------------------------------------
//Area and rho
//--------------------------------------------------
void plot_area_and_rho(systems system=cent, bool doToymodel=0, int toytype=0, TString version="GPC2", TString label="THIS THESIS", bool is4paper=0,bool pythiaEmb=0, TString ext="pdf")
{
	
	//bool kTR04=1;
	//bool alexpico=1;
	
	int nr=3; //(is4paper) ? 2 : 3;
	const int nR=nr;
	
	float R[]={0.2,0.3,0.4}; //jet resolution parameter
	float Acut[]={0.07,0.2,0.4}; //jet area cut
	float area_scale[]={0.37,0.25,0.22}; //for normalization range 
	//float area_scale[]={1,1,1}; //for normalization range 
	
	bool showInfo[]={1,1,1};
	bool showLeg[]={1,1,1};
	bool showLabel[]={0,0,0};
	bool showWatermark[]={0,0,0};
	
	Color_t cut_col=kGray+1; //area cut line color
	//for the paper we don't need R=0.3 -> replace it by R=0.4
	/*
	if(is4paper)
	{
		R[1]=R[2];
		Acut[1]=Acut[2];
		area_scale[1]=area_scale[2];
		showInfo[1]=showInfo[2];
		showLeg[1]=showLeg[2];
		showLabel[1]=showLabel[2];
		showWatermark[1]=showWatermark[2];
	}
*/
		
	
	TFile *fin[nR]; 
	TFile *fin2[nR]; 
	TH1D* hjetarea[nR];
	TH2D* hjetpTarea[nR];
	TH1D* hrho[nR];
	TH1D* harea_emb5[nR];
	TH1D* harea_emb10[nR];
	TH1D* harea_emb20[nR];
	TString base_name[2]={"star","toy"};
	TString cent_name[2]={"central","peripheral"};
	//TString subdir[2]={"kTRr","kTR04"};
	TString subdirToy[3]={"jet_plus_bg","jetonly","boltzman"};
	//TString pico[2]={"mypico","alexpico"};
	TString path=Form("../../plotting_out/intersteps/area_and_rho/%s_%s",base_name[doToymodel].Data(),cent_name[system].Data());
	TString out_dir=Form("../../plotting_out/obr/%s/area/",version.Data());
	if(!doToymodel)
	{
		//path+=Form("/%s_%s",subdir[kTR04].Data(),pico[alexpico].Data());
		//out_dir+=Form("/%s_%s",subdir[kTR04].Data(),pico[alexpico].Data());
		path+=Form("/%s",version.Data());
	}
	else if(doToymodel)
	{
		path+=Form("/%s",subdirToy[toytype].Data());
		out_dir+=Form("/toy/%s",subdirToy[toytype].Data());
	}

	for(int r=0; r<nR; r++)
	{
		//TString fname=Form("%s/histos_inclusivejet_R%.1lf.root",path.Data(),R[r]);
		TString fname=Form("%s/histos_inclusivejet_RRho0%.0lf.root",path.Data(),10*R[r]);
		if(doToymodel) fname=Form("%s/histos_jets_R%.1lf_pTcut0.2.root",path.Data(),R[r]);
		cout<<"opening file:"<<fname.Data()<<endl;
		fin[r]=new TFile(fname.Data(), "OPEN");

		TString sarea="hjetarea";
		TString sjetpTarea="hjetpTcorrArea";
		TString srho="hrho";
		if(!doToymodel)
		{
			TString sR=Form("_R0%.0lf",10*R[r]);
			sarea+=sR;
			sjetpTarea+=sR;
			srho+=sR;
		}
		
		hjetarea[r]=(TH1D*) fin[0]->Get(sarea);
		hjetarea[r]->Sumw2();
		hjetarea[r]->Scale(1./hjetarea[r]->Integral());
		TH2D* htmp=(TH2D*) fin[0]->Get(sjetpTarea);

		hjetpTarea[r]=new TH2D(sjetpTarea,"jet p_{T} vs. area; jet area (Sr); jet p_{T}",htmp->GetNbinsY(),htmp->GetYaxis()->GetBinLowEdge(1),htmp->GetYaxis()->GetBinUpEdge(htmp->GetNbinsY()),htmp->GetNbinsX(),htmp->GetXaxis()->GetBinLowEdge(1),htmp->GetXaxis()->GetBinUpEdge(htmp->GetNbinsX()));
		for(int i=1; i<=htmp->GetNbinsX();i++){
		for(int j=1; j<=htmp->GetNbinsY();j++){
			double val=htmp->GetBinContent(i,j);
			hjetpTarea[r]->SetBinContent(j,i,val);
		}}
		delete htmp;
		
		hrho[r]=(TH1D*) fin[r]->Get(srho);
		float rhoint=hrho[r]->Integral();
		hrho[r]->Sumw2();
		hrho[r]->SetLineColor(G_colorList[r]);
		hrho[r]->SetLineWidth(G_linewidth);
		hrho[r]->SetLineStyle(G_lineStyle[r]);
		hrho[r]->SetFillColor(G_colorFill[r]);
		hrho[r]->SetFillStyle(G_styleFill[r]);

		hrho[r]->Scale(1./rhoint);
		
		if(!doToymodel){
		fname=Form("%s/histos_embeddedjet.root",path.Data());
		if(pythiaEmb) fname=Form("%s/histos_embeddedjet_pythiaU.root",path.Data());
		//fname=Form("%s/histos_embeddedjet_R%.1lf.root",path.Data(),R[r]);
		fin2[r]=new TFile(fname.Data(), "OPEN");
		TString sembarea=Form("hjetpTembArea%s",sR.Data());
		//TString sembarea="hjetpTembArea";
		TH2D* hareaemb2d=(TH2D*) fin2[r]->Get(sembarea);
		int bin5=hareaemb2d->GetXaxis()->FindBin(5);
		int bin10=hareaemb2d->GetXaxis()->FindBin(10);
		int bin20=hareaemb2d->GetXaxis()->FindBin(20);
		harea_emb5[r]=hareaemb2d->ProjectionY(Form("harea_emb5_r%i",r),bin5,bin5);
		harea_emb10[r]=hareaemb2d->ProjectionY(Form("harea_emb10_r%i",r),bin10,bin10);
		harea_emb20[r]=hareaemb2d->ProjectionY(Form("harea_emb20_r%i",r),bin20,bin20);
		harea_emb5[r]->Sumw2();
		harea_emb5[r]->Scale(area_scale[r]/harea_emb5[r]->Integral());
		harea_emb10[r]->Sumw2();
		harea_emb10[r]->Scale(area_scale[r]/harea_emb10[r]->Integral());
		harea_emb20[r]->Sumw2();
		harea_emb20[r]->Scale(area_scale[r]/harea_emb20[r]->Integral());
		delete hareaemb2d;
			
		}
	}
	
	
	const int can_x=500; //1600
	const int can_y=500; //900
	const  float pTXmin=-20;
	const float pTXmax=50;
	const  float areaXmin=0;
	const float areamax[nR]={0.4,0.6,1.0};
	const float areaXmax[nR]={0.25,0.5,0.9};
	const float areaYmin=1E-4;
	float areaYmax[nR];
	float rhoXmin=15;
	float rhoXmax=45;
	float rhoYmin=0.0;
	float rhoYmax=0.2;
	if(system==peri)
	{
		rhoXmin=-15;
		rhoXmax=15;
		rhoYmax=1.1;
	}

	
	Color_t colorList[]={kRed,kBlue,kBlack,kMagenta};
	
	TString title_2d_y="p_{T,jet}^{reco} (GeV/#it{c})";
	TString title_2d_x="Area (Sr)";
	
	TString title_area_x="Area (Sr)";
	TString title_area_y="arb. units";
		
	TString title_rho_x="#rho (GeV/(c.Sr))";
	TString title_rho_y="arb. units";
	
	TLatex *latex = new TLatex();
	latex->SetNDC();
	latex->SetTextSize(0.04);
      
	gStyle->SetPadRightMargin(0.10);

	//areas vs jet pT
	    
	TLatex *latexW = new TLatex();
	latexW->SetNDC();
	latexW->SetTextSize(G_latex_sz);
	latexW->SetTextColor(G_latex_cl);
	
	//area vs pTjet
	TCanvas *c1[nR];
	//c1->Divide(nR,1);
	for(int r=0;r<nR;r++)
	{
		c1[r]= new TCanvas(Form("c1_%i",r),"area vs. pT",10,10,1.2*can_x,can_y);
		//c1->cd(r+1);
		c1[r]->cd();
		
		double lx1,lx2,ly1,ly2;
		lx1=0.20;  lx2=0.45; ly1=0.7; ly2=0.9;
		if (r==0) {lx1=0.55;  lx2=0.8;}
		TLegend *model_info = new TLegend(lx1,ly1,lx2,ly2);
		model_info->SetTextSize(0.045);
		model_info->SetFillStyle(0);
		model_info->SetBorderSize(0);
		model_info->SetMargin(0.05);
		if(doToymodel)
			model_info->SetHeader("Parametrized model");
		else
			model_info->SetHeader("Au+Au #sqrt{s_{NN}} = 200 GeV");
		model_info->AddEntry("",G_system_info[system].Data(),"");
		if(showLabel[r]) model_info->AddEntry("", "anti-k_{T}", "");
		else model_info->AddEntry("", Form("anti-k_{T}, R = %.1lf",R[r]), "");
			//model_info->AddEntry("", Form("%iM events",mevents), "");
		
		gPad->SetLogz();
		hjetpTarea[r]->GetYaxis()->SetRangeUser(pTXmin, pTXmax);
		hjetpTarea[r]->GetXaxis()->SetRangeUser(areaXmin, areamax[r]);
		//hjetpTarea[r]->GetYaxis()->SetTitleOffset(1.2);
		hjetpTarea[r]->GetXaxis()->SetTitleSize(0.06);
		hjetpTarea[r]->GetYaxis()->SetTitleSize(0.06);
		hjetpTarea[r]->GetXaxis()->SetTitle(title_2d_x);
		hjetpTarea[r]->GetYaxis()->SetTitle(title_2d_y);
		hjetpTarea[r]->SetTitle("");
		hjetpTarea[r]->Draw("COLZ");
		if(showLabel[r]) latex->DrawLatex(0.7, 0.8,Form("R = %.1lf",R[r]));
		
		TLine *l_areacut = new TLine(Acut[r], pTXmin, Acut[r], pTXmax);
		l_areacut->SetLineWidth(2);
		l_areacut->SetLineStyle(2);
		l_areacut->SetLineColor(cut_col);
		l_areacut->DrawClone("same");
		
		if(showInfo[r]) model_info->DrawClone();
		if(showWatermark[r]) latexW->DrawLatex(0.55, 0.8,label);

	
		TString ostr=Form("%s/%s_%s_area_vs_pT_R0%.0lf.%s",out_dir.Data(),base_name[doToymodel].Data(),cent_name[system].Data(),R[r]*10,ext.Data());
		c1[r]->SaveAs(ostr.Data());
	}
		int mstyle=29;
		
	//gStyle->SetPadRightMargin(0.02);

	//jet areas
	if(!doToymodel)
	{
	//TCanvas *c2 = new TCanvas("c2","area",10,10,nR*can_x,can_y);
	//c2->Divide(nR,1);
	TCanvas *c2[nR];
	//xy coordinates of the latex label
	float latx[nR]={0.2,0.3,0.3};
	float laty[nR]={0.60,0.62,0.6};
	
	for(int r=0;r<nR;r++)
	{
		c2[r]= new TCanvas(Form("c2_%i",r),"area",10,10,1.2*can_x,can_y);
		c2[r]->cd();
		
		TLegend *model_info = new TLegend(0.20, 0.7, 0.45, 0.900);
		model_info->SetTextSize(0.042);
		model_info->SetFillStyle(0);
		model_info->SetBorderSize(0);
		model_info->SetMargin(0.05);
		if(doToymodel)
			model_info->SetHeader("Parametrized model");
		else
			model_info->SetHeader("Au+Au #sqrt{s_{NN}} = 200 GeV");
		model_info->AddEntry("",G_system_info[system].Data(),"");
		if(showLabel[r]) model_info->AddEntry("", "anti-k_{T}", "");
		else model_info->AddEntry("", Form("anti-k_{T}, R = %.1lf",R[r]), "");
		
		
		//gPad->SetLogy();
		int binmax = hjetarea[r]->GetMaximumBin();
		areaYmax[r]=1.75*hjetarea[r]->GetBinContent(binmax);
/*
		harea_emb5[r]->GetYaxis()->SetRangeUser(areaYmin, areaYmax[r]);
		harea_emb5[r]->GetXaxis()->SetRangeUser(areaXmin, areaXmax[r]);
		harea_emb5[r]->GetXaxis()->SetTitle(title_area_x);
		harea_emb5[r]->GetYaxis()->SetTitle(title_area_y);
		harea_emb5[r]->GetYaxis()->SetTitleOffset(1.25);
		harea_emb5[r]->GetXaxis()->SetTitleSize(0.05);
		harea_emb5[r]->GetYaxis()->SetTitleSize(0.045);
		harea_emb5[r]->SetTitle("");*/
		harea_emb5[r]->SetLineWidth(2);
		harea_emb5[r]->SetLineStyle(2);
		harea_emb5[r]->SetLineColor(kRed+2);
		harea_emb5[r]->SetMarkerSize(1);
		harea_emb5[r]->SetMarkerColor(kGray);
		harea_emb5[r]->SetMarkerStyle(mstyle);
		//harea_emb5[r]->SetFillColor(kGray);
		//harea_emb5[r]->SetFillStyle(3001);
		
		harea_emb10[r]->SetLineColor(kRed);
		harea_emb10[r]->SetLineWidth(1);
		harea_emb10[r]->SetLineStyle(1);
		harea_emb10[r]->SetMarkerSize(1);
		harea_emb10[r]->SetMarkerColor(kRed);
		harea_emb10[r]->SetMarkerStyle(mstyle);
		harea_emb10[r]->SetFillColor(kGray);
		harea_emb10[r]->SetFillStyle(3003);
		
		harea_emb20[r]->SetLineColor(kMagenta);
		harea_emb20[r]->SetLineWidth(2);
		harea_emb20[r]->SetLineStyle(4);
		harea_emb20[r]->SetMarkerSize(1);
		harea_emb20[r]->SetMarkerColor(kMagenta);
		harea_emb20[r]->SetMarkerStyle(mstyle);
		
		hjetarea[r]->GetYaxis()->SetRangeUser(areaYmin, areaYmax[r]);
		hjetarea[r]->GetXaxis()->SetRangeUser(areaXmin, areaXmax[r]);
		hjetarea[r]->GetXaxis()->SetTitle(title_area_x);
		hjetarea[r]->GetYaxis()->SetTitle(title_area_y);
		hjetarea[r]->GetYaxis()->SetTitleOffset(1.3);
		hjetarea[r]->GetXaxis()->SetTitleSize(0.05);
		hjetarea[r]->GetYaxis()->SetTitleSize(0.055);
		hjetarea[r]->SetTitle("");
		hjetarea[r]->SetLineWidth(1);
		hjetarea[r]->SetLineColor(kBlue);
		hjetarea[r]->SetMarkerSize(1);
		hjetarea[r]->SetMarkerColor(kBlue);
		hjetarea[r]->SetMarkerStyle(mstyle);
		
		hjetarea[r]->DrawCopy("e");
		//harea_emb5[r]->DrawCopy("histo same");
		harea_emb10[r]->DrawCopy("histo same");
		harea_emb20[r]->DrawCopy("histo same");
		hjetarea[r]->DrawCopy("esame");

		
		
		if(showLabel[r]) latex->DrawLatex(latx[r], laty[r],Form("R = %.1lf",R[r]));
		
		float ymax=harea_emb5[r]->GetBinContent(hjetarea[r]->FindBin(R[r]*R[r]*3.14));
		TLine *l_areacut = new TLine(Acut[r], areaYmin, Acut[r], ymax);
		l_areacut->SetLineWidth(2);
		l_areacut->SetLineStyle(2);
		l_areacut->SetLineColor(cut_col);
		l_areacut->DrawClone("same");
		
		if(showInfo[r]) model_info->DrawClone();
		if(showWatermark[r]) latexW->DrawLatex(0.4, 0.2,label);
		
		if(showLeg[r]){
		  TLegend *legspectra = new TLegend(0.5, 0.6, 0.95, 0.90);
			legspectra->SetTextSize(0.040);
			legspectra->SetFillStyle(0);
			legspectra->SetBorderSize(0);
			legspectra->AddEntry("",   "STAR data:", "");
			legspectra->AddEntry(hjetarea[r],   "  inclusive jets", "lp");
			//legspectra->AddEntry(harea_emb5[r], "embedded jets, p_{T}^{emb}=5GeV/#it{c}", "l");
			TString embdesc=(pythiaEmb) ?  "Embed. (PYTHIA):" : "Embedding (SP):";
			legspectra->AddEntry("", embdesc , "");
			legspectra->AddEntry(harea_emb10[r], "  p_{T}^{emb}=10 GeV/#it{c}", "l");  //Jana
			legspectra->AddEntry(harea_emb20[r], "  p_{T}^{emb}=20 GeV/#it{c}", "l");  //Jana
			

			legspectra->DrawClone("same");}
  
			TString frag_suf=(pythiaEmb) ? "_pythia" : "";
		ostr=Form("%s/%s_%s_area_pTemb_R0%.0lf%s.%s",out_dir.Data(),base_name[doToymodel].Data(),cent_name[system].Data(),R[r]*10,frag_suf.Data(),ext.Data());
		c2[r]->SaveAs(ostr.Data());
	}//R loop
	
	
	}//!doToymodel
	
	
	//rho
	TCanvas *c3 = new TCanvas("c3","rhos",10,10,1.5*can_x,can_y);
	c3->cd();
	
	TH1 *frame = new TH1I("frame", "", 1000, -100, +100);
	frame->GetXaxis()->SetRangeUser(rhoXmin, rhoXmax);
	frame->GetYaxis()->SetRangeUser(rhoYmin, rhoYmax);
	frame->GetXaxis()->SetTitle(title_rho_x);
	frame->GetYaxis()->SetTitle(title_rho_y);
	frame->GetYaxis()->SetTitleOffset(0.95);
	frame->SetTitle("");
	frame->Draw();
	for(int r=0;r<nR;r++)
	{
		hrho[r]->Draw(G_histoDraw[r]);
	}
	
	TLegend *model_info = new TLegend(0.20, 0.7, 0.35, 0.900);
	model_info->SetTextSize(0.04);
	model_info->SetFillStyle(0);
	model_info->SetBorderSize(0);
	model_info->SetMargin(0.05);
	if(doToymodel)
		model_info->SetHeader("Parametrized model");
	else
		model_info->SetHeader("Au+Au #sqrt{s_{NN}} = 200 GeV");
	model_info->AddEntry("",G_system_info[system].Data(),"");
	model_info->AddEntry("", "k_{T} algorithm", "");
	model_info->DrawClone();
	
	TLegend *legspectra = new TLegend(0.6467, 0.60, 0.89, 0.90);
	legspectra->SetTextSize(0.04);
	legspectra->SetFillStyle(0);
	legspectra->SetBorderSize(0);
	for(int r=0;r<nR;r++)
	{
		legspectra->AddEntry(hrho[r], Form("R = %.1lf",R[r]), G_legDraw[r]);
	}
	legspectra->DrawClone("same");
	latexW->DrawLatex(0.45, 0.2,label);
	
	ostr=Form("%s/%s_%s_rho.%s",out_dir.Data(),base_name[doToymodel].Data(),cent_name[system].Data(),ext.Data());
	c3->SaveAs(ostr.Data());
	
	
}



//--------------------------------------------------
//ratios of spectra with different RHO calculation
//--------------------------------------------------

void plot_ratios_rho(float R=0.2, TString suf="kTR04_alexpico", TString ext="gif")
{
	
	TString out_dir="../../plotting_out/obr/intersteps/area";
	TFile *fin[3];
	TH1D* hspec[3];
	TString str[3];
	str[0]=Form("../../plotting_out/intersteps/area_and_rho/star_central/kTR04_alexpico/histos_inclusivejet_R%.1lf.root",R);
	str[1]=Form("../../plotting_out/intersteps/area_and_rho/star_central/kTRr_mypico/histos_inclusivejet_R%.1lf.root",R);
	str[2]=Form("../../plotting_out/intersteps/area_and_rho/star_central/%s/histos_inclusivejet_R%.1lf.root",suf.Data(),R);
	
	for(int i=0;i<3;i++)
	{
		fin[i]=new TFile(str[i].Data(), "OPEN");
		hspec[i]=(TH1D*) fin[i]->Get("hpT_pTl5");
		hspec[i]->Scale(1.0/hspec[i]->Integral());
	}
	
	TH1D* hrat_good=(TH1D*) hspec[0]->Clone("hratgood");
	TH1D* hrat_bad=(TH1D*) hspec[1]->Clone("hratbad");
	hrat_good->Divide(hspec[2]);
	hrat_bad->Divide(hspec[2]);
	
	TH1D* hrat_goodbad=(TH1D*) hspec[0]->Clone("hratgoodbad");
	hrat_goodbad->Divide(hspec[1]);
	
	const int can_x=1200; //1600
	const int can_y=680; //900
	
	TCanvas *c1 = new TCanvas("c1","ratio",10,10,can_x,can_y);
	c1->cd();
	hrat_good->GetXaxis()->SetRangeUser(-20, 50);
	hrat_good->Draw();
	hrat_bad->SetLineColor(kRed);
	hrat_bad->Draw("same");
	
	TLegend *legspectra = new TLegend(0.6467, 0.60, 0.89, 0.90);
	legspectra->SetTextSize(0.04);
	legspectra->SetFillStyle(0);
	legspectra->SetBorderSize(0);
	legspectra->AddEntry(hrat_bad, "bad/selected", "lp");
	legspectra->AddEntry(hrat_good, "good/selected", "lp");
	legspectra->DrawClone("same");
	
	TCanvas *c2 = new TCanvas("c2","ratio",10,10,can_x,can_y);
	c2->cd();
	TH1 *frame = new TH1I("frame", "", 1000, -100, +100);
	frame->GetXaxis()->SetRangeUser(-20, 50);
	frame->GetYaxis()->SetRangeUser(0.1, 3);
	frame->GetXaxis()->SetTitle("p_{T}-A#rho (GeV/#it{c})");
	frame->SetTitle("");
	frame->Draw();
	hrat_goodbad->SetLineColor(kBlue);
	hrat_goodbad->SetLineStyle(1);
	hrat_goodbad->SetLineWidth(2);
	hrat_goodbad->Rebin(2);
	hrat_goodbad->Scale(1./2.);
	hrat_goodbad->Draw("same");
	
	TLine *one = new TLine(-20, 1, 50, 1);
	one->SetLineWidth(1);
	one->SetLineStyle(2);
	one->SetLineColor(kBlack);
	one->DrawClone("same");
	
	TLegend *legspectra = new TLegend(0.4, 0.60, 0.89, 0.90);
	legspectra->SetTextSize(0.04);
	legspectra->SetFillStyle(0);
	legspectra->SetBorderSize(0);
	legspectra->AddEntry(hrat_goodbad, Form("(p_{T}-A#rho_{R = 0.4})/(p_{T}-A#rho_{R = %.1lf})",R), "lp");
	legspectra->DrawClone("same");
	
	TLegend *model_info = new TLegend(0.20, 0.48, 0.35, 0.900);
	model_info->SetTextSize(0.035);
	model_info->SetFillStyle(0);
	model_info->SetBorderSize(0);
	model_info->SetMargin(0.05);
	model_info->SetHeader("Run 11 Au+Au #sqrt{s_{NN}} = 200 GeV");
	model_info->AddEntry("", "0-10% most central", "");
	model_info->AddEntry("", Form("anti-k_{T}, R = %.1lf",R), "");
	model_info->DrawClone("same");
	
	TString ostr=Form("%s/rho_ratio_R0%.0lf_pTlead0.%s",out_dir.Data(),R*10,ext.Data());
	c2->SaveAs(ostr.Data());
}

//==========================================================
//jet reconstruction efficiency
//==========================================================

plot_epsilon(float R=0.4, float pTlead=5.0, systems system=cent, TString version="GPC1", TString label="",TString ext="gif" /*output figure format*/)
{
	
	bool useNewEpsilon=1; //set to 0 for older files
	bool showTrackingEff=1; //show also tracking efficiency
	
	const int can_x=1000; //1600
	const int can_y=800; //900
	const  float spectraXmin=0;
	const float spectraXmax=50;
	const float spectraYmin=0.6;
	const float spectraYmax=1.3;
	
	const int nPanels=4; //max number of columns
	const int nRows=3; //maximum number of rows
	int ridx=calc_R_index(R);
	
	//in which panels do we want to show info, legend and labels?
	bool showInfo[nRows][nPanels]={/*row1*/{1,1,1,0},
	/*row2*/{1,1,1,0},
	/*row3*/{1,1,1,0}};
	bool showLeg[nRows][nPanels]={/*row1*/{1,1,1,0},
	/*row2*/{1,1,1,0},
	/*row3*/{1,1,1,0}};
	bool showLabel[nRows][nPanels]={/*row1*/{1,1,1,0},
	/*row2*/{1,1,1,0},
	/*row3*/{1,1,1,0}};
	bool showWatermark[nRows][nPanels]={/*row1*/{0,0,0,0},
	/*row2*/{0,0,0,0},
	/*row3*/{0,0,0,0}};
	
	short panelType=1; //what does show different panels? 0: different Rs | 1: different centralities
	int pidx=(panelType==0)? ridx : system; //panel index
	int rowidx=(panelType==1)? ridx : system; //row index
	
	//set description of panels and descriptor in information legend based on the panel type
	TString descR=Form("R = %.1lf",R); 
	TString descC=G_system_info[system];
	/*TString descPanel=descR;//(panelType==0)? descR : descC;
	 *	TString descInfo=descC; //(panelType==1)? descR : descC;
	 *	float panx[]={0.5,0.5,0.5,0.5}; //x-position of the panel desc.
	 *	float pany=0.8;//y-position of the panel desc.*/
	
	TString out_dir=Form("../../plotting_out/obr/%s/effi",version.Data());
	
	int binmerge=4; //how many bins merge together
	TString jettype="u+g (2:1)"; //jet fragmentation model
	
	const int ntypes=5;
	TString type[ntypes]={"normal","u","g","m5","p5"};
	TString cent_name[]={"central","peripheral"};
	Color_t color[ntypes]={kBlue, kRed-7, kMagenta, kTeal,kCyan}; //line color
	Color_t fillColor[ntypes]={kGray+1,kGray,kGray,kGray,kGray};
	int lstyle[ntypes]={1,2,2,4,4}; //line style
	int linewidth=3;
	int fillStyle[ntypes]={0,0,0,0,0};
	TString legDraw[]={"l","l","l","l","l"};
	
	TString leg_desc[ntypes]={"nominal","u-quark","gluon","#epsilon_{track}-5%","#epsilon_{track}+5%"};
	
	TFile *fin[ntypes];
	TGraph* hepsilon[ntypes];
	
	double pTmax=spectraXmax;
	double pTmin=pTlead+1; //show only bins ~1 GeV above pTlead cut, lower bins are dropping to zero and we don't want to show them
	double binWidth=0.5;
	int nbins=(int) (pTmax-pTmin)/binWidth; //number of bins, 0.5 GeV bin size
	const int npoints=nbins;
	double grY[npoints];
	double grX[npoints];
	
	for(int t=0;t<ntypes;t++)
	{
		cout<<"opening file"<<endl;
		TString str;
		if(!useNewEpsilon) str = Form("../../plotting_out/intersteps/epsilon/%s/%s/%s/epsilon_R%.1lf_pTlead%.0lf.root",version.Data(),cent_name[system].Data(),type[t].Data(),R,pTlead);
		else str = Form("../../plotting_out/intersteps/epsilon/%s/%s/%s/pythia_emb_R%.1lf.root",version.Data(),cent_name[system].Data(),type[t].Data(),R);
		fin[t] = new TFile(str.Data(), "OPEN");
		TH1D* htmp;
		if(!useNewEpsilon)htmp=(TH1D*) fin[t]->Get("hepsilon_unfolded");
		else htmp=(TH1D*) fin[t]->Get(Form("hepsilon_pTl%.0lf",pTlead));
		//hepsilon[t]=new TH1D(Form("hepsilon_%i",t),"efficiency",nbins,pTmin,pTmax);
		
		/*
		 *		//correct error scale (Sumw2 was not called before scaling original histograms by the number of jobs=20)
		 *		hepsilon[t]->Scale(20);
		 *		hepsilon[t]->Sumw2();
		 *		hepsilon[t]->Scale(1.0/20);*/
		htmp->Rebin(binmerge);
		htmp->Scale(1./binmerge);
		
		//for the final histogram, copy only bins ~1 GeV above pTlead cut, lower bins are dropping to zero and we don't want to show them
		for(int bn=1;bn<nbins+1;bn++)
		{
			double pTc=pTmin+(bn-0.5/*start in the center of the bin*/)*binWidth;
			int bntmp=htmp->GetXaxis()->FindBin(pTc);
			double yval=htmp->GetBinContent(bntmp);
			//hepsilon[t]->SetBinContent(bn,yval);
			grY[bn-1]=yval;
			grX[bn-1]=pTc;
		}
		delete htmp;
		hepsilon[t]=new TGraph(npoints,grX,grY);
		hepsilon[t]->SetLineColor(color[t]);
		hepsilon[t]->SetLineStyle(lstyle[t]);
		hepsilon[t]->SetLineWidth(linewidth);
		//hepsilon[t]->SetFillStyle(fillStyle[t]);
		//hepsilon[t]->SetFillColor(fillColor[t]);
	}
	
	TCanvas *ceps = new TCanvas("ceps","jet reco. efficiency",10,10,can_x,can_y);
	ceps->cd();
	//ceps->SetLogy();
	//ceps->SetGrid();
	
	TLatex *latexW = new TLatex();
	latexW->SetNDC();
	latexW->SetTextSize(G_latex_sz);
	latexW->SetTextColor(G_latex_cl);
	float labx=0.42;
	float laby=0.25;
	
	TLatex *latex = new TLatex();
	latex->SetNDC();
	latex->SetTextSize(0.04);
	
	TH1 *frame = new TH1I("frame", "", 1000, -100, +100);
	//frame->GetXaxis()->SetTitle("p_{T, true}^{ch jet} (GeV/#it{c})");
	frame->GetXaxis()->SetTitle("p_{T, jet}^{part} (GeV/#it{c})");
	frame->GetYaxis()->SetTitle("#varepsilon_{jet}");
	frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
	frame->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
	//frame->GetYaxis()->SetRangeUser(1E-3, 1E5);
	TString str="";
	frame->SetTitle(str);
	frame->DrawCopy("");
	
	TLine *one = new TLine(spectraXmin, 1, spectraXmax, 1);
	one->SetLineWidth(2);
	one->SetLineStyle(2);
	one->SetLineColor(kGray+1);
	one->DrawClone("same");
	
	bool central=0;
	if(system==cent)central=1;
	TF1* feff_Alex=efficiency11(central,"AuAu");
	feff_Alex->SetLineColor(kOrange+7);
	if(showTrackingEff) feff_Alex->Draw("same");
	
	for(int t=0;t<ntypes;t++)
	{
		hepsilon[t]->Draw("L");
	}
	
	//legends and labels
	TLegend *legspectra = new TLegend(0.6467, 0.60, 0.89, 0.90);
	legspectra->SetTextSize(0.04);
	legspectra->SetFillStyle(0);
	legspectra->SetBorderSize(0);
	for(int t=0;t<ntypes;t++)
	{
		legspectra->AddEntry(hepsilon[t], leg_desc[t].Data(), legDraw[t]);
	}
	legspectra->AddEntry(""," ","");
	if(showTrackingEff)legspectra->AddEntry(feff_Alex, "tracking eff.", "l");
	
	TLegend *model_info = new TLegend(0.20, 0.6, 0.35, 0.900);
	model_info->SetTextSize(0.038);
	model_info->SetFillStyle(0);
	model_info->SetBorderSize(0);
	model_info->SetMargin(0.05);
	if(showInfo[rowidx][pidx])model_info->SetHeader(Form("PYTHIA %s jets",jettype.Data()));
	else model_info->SetHeader(" ");
	model_info->AddEntry("",descC,"");
	model_info->AddEntry("", Form("anti-k_{T}, %s",descR.Data()), "");
	//	if(showInfo[rowidx][pidx])model_info->AddEntry("", Form("p_{T,lead}^{min} = %.1lf GeV/#it{c}",pTlead), ""); //Jana
	if(showInfo[rowidx][pidx])model_info->AddEntry("", Form("p_{T,lead}^{min} = %.0lf GeV/#it{c}",pTlead), ""); //Jana
	
	if(showLeg[rowidx][pidx]) legspectra->DrawClone("same");
	if(showWatermark[rowidx][pidx]) latexW->DrawLatex(labx,laby,label);
	if(showLabel[rowidx][pidx]) model_info->DrawClone("same");
	//if(showLabel[rowidx][pidx])latex->DrawLatex(panx[ridx],pany,descPanel);
	
	str=Form("%s/epsilon_R0%.0lf_pTlead%.0lf_%s.%s",out_dir.Data(),R*10,pTlead,cent_name[system].Data(),ext.Data());
	ceps->SaveAs(str.Data());
}

//==========================================================
//delta-pT
//==========================================================
plot_dpT(float Rpar=0.4, float pTlead=0, bool doToymodel=0, systems system=cent,TString version="GPC2", TString label="THIS THESIS", TString ext="pdf")
{
	
	const int can_x=1000; //1600
	const int can_y=800; //900
	const  float spectraXmin=-25;
	const float spectraXmax=30;
	const float spectraYmin=1E-5;
	const float spectraYmax=2.0;
	const int pyt_emb=2; //which probe pT to draw for PYTHIA jets
	int binmerge=4; //how many bins merge together
	
	const int nR=4;
	float R[nR]={0.2,0.3,0.4,0.5};
	int ridx=calc_R_index(Rpar);
	
	const int nPanels=4; //maximum number of panels in one row
	const int nRows=3; //maximum number of rows
	
	//in which panels do we want to show info, legend and labels?
	bool showInfo[nRows][nPanels]={/*row1*/{1,1,1,0},
	/*row2*/{1,1,1,0},
	/*row3*/{1,1,1,0}};
	bool showLeg[nRows][nPanels]={/*row1*/{1,1,1,0},
	/*row2*/{1,1,1,0},
	/*row3*/{1,1,1,0}};
	bool showLabel[nRows][nPanels]={/*row1*/{1,1,1,0},
	/*row2*/{1,1,1,0},
	/*row3*/{1,1,1,0}};
	bool showWatermark[nRows][nPanels]={/*row1*/{1,1,1,0},
	/*row2*/{1,1,1,0},
	/*row3*/{1,1,1,0}};
	short panelType=0; //what does show different panels? 0: different Rs | 1: different centralities
	int pidx=(panelType==0)? ridx : system; //panel index
	int rowidx=(panelType==1)? ridx : system; //row index
	
	//set description of panels and descriptor in information legend based on the panel type
	TString descR=Form("R = %.1lf",Rpar); 
	TString descC=G_system_info[system];
	/*TString descPanel=descR;//(panelType==0)? descR : descC;
	 *	TString descInfo=descC;//(panelType==1)? descR : descC;
	 *	float panx[]={0.75,0.25,0.25,0.25}; //x-position of the panel desc.
	 *	float pany=0.8;//y-position of the panel desc.*/
	
	//compare deltapT for different pTemb
	const int nemb =3;
	double fEmbPt[] = {5.0, 10.0, 20.0};
	
	//compare deltapT for different probes (SP, PYTHIA, v2 corrected SP)
	const int ntypes=3;
	double pTEmb_show[ntypes]={20.0,20,20.0}; //which pTembed use for these delta-pTs
	TString type[]={"sp", "sp", "pythiaU"};
	TString leg_desc[]={"SP, #it{v}_{2}=0", "SP #otimes v_{2}=0.04", "PYlq, #it{v}_{2}=0"};
	short v2type=1; //which type is the one with v2 correction?
	
	//lines
	const float linewidth=3.0;
	Color_t colorList[]={kRed,kGray+3,kBlue,kBlack};
	int lineStyle[]={1,3,2,4};
	
	//fills
	Color_t colorFill[]={kGray,kGray,kGray,kGray};
	int styleFill[]={0,0,0,0};
	//markers
	const float marksize=1.0;
	int markerList[]={20,21,22,23,29,33,34};
	
	//draw options
	TString histoDraw[nemb]={"histo same","histo same","histo same"};
	TString legDraw[nemb]={"l","l","l"};
	
	TFile *fin[ntypes];
	TFile *fv2;
	
	TH2D *h2dpt[nR][ntypes];
	TH2D *h2dpt2[nR];
	TH2D *hv2[nR]; //v2 corrections
	TH1D* hdpt1[nemb];
	TH1D* hdpt2[ntypes];
	
	TString out_dir=Form("../../plotting_out/obr/%s/dpT",version.Data());
	TString data_names[]={"MB","toymodel"};
	TString str;
	
	TString name=data_names[doToymodel];
	if(system==peri)name+="_peripheral";
	else if(system==cent) name+="_central";
	
	//open files
	str="../../plotting_out/intersteps/embedding/v2corr.root"; //file with v2 correction factors for delta-pT
	fv2=new TFile(str.Data(), "OPEN");
	for(int t=0;t<ntypes;t++){
		str = Form("../../plotting_out/intersteps/embedding/%s/%s/histos_embeddedjet_%s.root",name.Data(),version.Data(),type[t].Data());
		cout<<"opening file "<<str.Data()<<endl;
		fin[t] = new TFile(str.Data(), "OPEN");
		for(int r=0;r<nR;r++)
		{
			//TString str = Form("../../plotting_out/intersteps/embedding/%s/histos_embeddedjet_R%.1lf.root",name.Data(),R[r]);
			
			str=Form("delta_pt_BG_sp_%0.lf_R0%.0lf",pTlead,R[r]*10);
			//cout<<"opening histogram "<<str.Data()<<endl;
			h2dpt[r][t]=(TH2D*)fin[t]->Get(str);
			
			if(t>0) continue;
			//v2 corrections
			short rv2=r;
			if(r>2) rv2=2; //we have v2 correction factors only up to R=0.4
			str=Form("dpt_v2corrVSuncorr_R0%.0lf",R[rv2]*10);
			if(system==cent) str+="_central";
			else str+="_peripheral";
			hv2[r]=(TH2D*)fv2->Get(str.Data());
			
			
		}}
		
		//delta pT histograms for one particular R
		for(int i=0; i<nemb; i++)
		{
			if(fEmbPt[i]<pTlead)continue;
			int bin=h2dpt[ridx][0]->GetXaxis()->FindBin(fEmbPt[i]);
			hdpt1[i]=(TH1D*)h2dpt[ridx][0]->ProjectionY(Form("hdpt1_%i",i),bin,bin);
			hdpt1[i]->SetLineColor(colorList[i]);
			hdpt1[i]->SetLineWidth(linewidth);
			hdpt1[i]->SetLineStyle(lineStyle[i]);
			hdpt1[i]->SetFillColor(colorFill[i]);
			hdpt1[i]->SetFillStyle(styleFill[i]);
			//hdpt1[i]->SetMarkerColor(colorList[i]);
			//hdpt1[i]->SetMarkerSize(marksize);
			//hdpt1[i]->SetMarkerStyle(markerList[i]);
			//hdpt1[i]->Rebin(binmerge);
			//hdpt1[i]->Scale(1./(binmerge*hdpt1[i]->Integral()));
			hdpt1[i]->Scale(1./(hdpt1[i]->Integral()));
		}
		
		
		
		for(int i=0; i<ntypes; i++)
		{
			int bin=h2dpt[ridx][i]->GetXaxis()->FindBin(pTEmb_show[i]);
			hdpt2[i]=(TH1D*)h2dpt[ridx][i]->ProjectionY(Form("hdpt2_%i",i),bin,bin);
			hdpt2[i]->SetLineColor(colorList[i]);
			hdpt2[i]->SetLineWidth(linewidth);
			hdpt2[i]->SetLineStyle(lineStyle[i]);
			hdpt2[i]->SetFillColor(colorFill[i]);
			hdpt2[i]->SetFillStyle(styleFill[i]);
			//hdpt2[i]->SetMarkerColor(colorList[i]);
			//hdpt2[i]->SetMarkerSize(marksize);
			//hdpt2[i]->SetMarkerStyle(markerList[i]);
			//hdpt2[i]->Rebin(binmerge);
			//hdpt2[i]->Scale(1./(binmerge*hdpt2[i]->Integral()));
			hdpt2[i]->Scale(1./(hdpt2[i]->Integral()));
			
			if(i!=v2type) continue;
			for(Int_t bin = 1; bin <= hdpt2[v2type]->GetNbinsX(); bin++)
			{
				Float_t pTreco=hdpt2[v2type]->GetBinCenter(bin);
				Float_t prob=hdpt2[v2type]->GetBinContent(bin);
				if(!prob>0)continue;
				Float_t prob_err=hdpt2[v2type]->GetBinError(bin);
				Int_t binx=hv2[ridx]->GetXaxis()->FindBin(fEmbPt[i]);
				Int_t biny=hv2[ridx]->GetYaxis()->FindBin(pTreco);
				Float_t scaler=hv2[ridx]->GetBinContent(binx,biny);
				if(!scaler>0)scaler=1;
				//cout<<"v2 scaling:"<<scaler<<endl;
				hdpt2[v2type]->SetBinContent(bin,prob*scaler);
				hdpt2[v2type]->SetBinError(bin,prob_err*scaler);
			}
		}
		
		
		//----------------------------------
		//Draw delta pT only for R=Rpar
		//----------------------------------
		
		//watermark
		TLatex *latexW = new TLatex();
		latexW->SetNDC();
		latexW->SetTextSize(G_latex_sz);
		latexW->SetTextColor(G_latex_cl);
		float labx=0.43;
		if(system==peri) labx=0.25;
		float laby=0.2;
		
		//panel desc.
		TLatex *latexP = new TLatex();
		latexP->SetNDC();
		latexP->SetTextSize(0.04);
		
		
		//*********************
		//different pTemb
		//*********************
		TCanvas *cdelta = new TCanvas("cdelta","delta pT",10,10,can_x,can_y);
		cdelta->cd();
		//cdelta->SetGrid();
		cdelta->SetLogy();
		
		TH1 *frame = new TH1I("frame", "", 1000, -100, +100);
		//frame->GetXaxis()->SetTitle("#delta p_{T} = p_{T, jet}^{measured}-#rho A_{jet} - p_{T, jet}^{emb}  (GeV/#it{c})");
		frame->GetXaxis()->SetTitle("#delta p_{T} (GeV/#it{c})");
		frame->GetYaxis()->SetTitle("probability/(GeV/#it{c})");
		frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
		frame->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
		//frame->GetYaxis()->SetRangeUser(1E-3, 1E5);
		TString str="";
		frame->SetTitle(str);
		frame->DrawCopy("");
		
		for(int i=0; i<nemb; i++)
		{
			if(fEmbPt[i]<pTlead)continue;
			hdpt1[i]->DrawCopy(histoDraw[i]);	
		}
		
		TLegend *legspectra = new TLegend(0.58, 0.60, 0.89, 0.90);
		legspectra->SetTextSize(0.04);
		legspectra->SetFillStyle(0);
		legspectra->SetBorderSize(0);
		for(int i=0; i<nemb; i++)
		{
			if(fEmbPt[i]<pTlead)continue;
			else legspectra->AddEntry(hdpt1[i], Form("p_{T}^{emb} = %.0lf GeV/#it{c}",fEmbPt[i]), legDraw[i]);
		}
		
		TLegend *model_info = new TLegend(0.20, 0.6, 0.35, 0.900);
		model_info->SetTextSize(0.038);
		model_info->SetFillStyle(0);
		model_info->SetBorderSize(0);
		model_info->SetMargin(0.05);
		if(showInfo[rowidx][pidx])
		{
			if(doToymodel)model_info->SetHeader("Parametrized model");
			else model_info->SetHeader("STAR Au+Au #sqrt{s_{NN}} = 200 GeV");
		}
		else model_info->SetHeader(" ");
		model_info->AddEntry("", descC, "");
		model_info->AddEntry("", Form("anti-k_{T}, %s",descR.Data()), "");
		if(pTlead>0) model_info->AddEntry("", Form("p_{T,lead}^{min} = %.0lf GeV/#it{c}",pTlead), "");
		
		if(showLabel[rowidx][pidx]) model_info->DrawClone("same");
		if(showWatermark[rowidx][pidx]) latexW->DrawLatex(labx, laby,label);
		if(showLeg[rowidx][pidx])legspectra->DrawClone("same");
		//if(showLabel[rowidx][pidx])latexP->DrawLatex(panx[pidx], pany,descPanel);
		
		str=Form("%s/dpT_%s_R0%.0lf_pTlead%.0lf.%s",out_dir.Data(),name.Data(),Rpar*10,pTlead,ext.Data());
		cdelta->SaveAs(str.Data());
		
		
		//*********************
		//different probes
		//*********************
		TCanvas *cdelta2 = new TCanvas("cdelta2","delta pT",10,10,can_x,can_y);
		cdelta2->cd();
		//cdelta->SetGrid();
		cdelta2->SetLogy();
		frame->DrawCopy("");
		for(int i=0; i<ntypes; i++)
		{
			hdpt2[i]->DrawCopy(histoDraw[i]);
		}
		
		TLegend *legspectra2 = new TLegend(0.60, 0.60, 0.89, 0.90);
		legspectra2->SetTextSize(0.04);
		legspectra2->SetFillStyle(0);
		legspectra2->SetBorderSize(0);
		legspectra2->SetHeader(Form("p_{T}^{emb} = %.0lf GeV/#it{c}",pTEmb_show[0]));
		for(int i=0; i<ntypes; i++)
		{
			legspectra2->AddEntry(hdpt2[i], Form("%s",leg_desc[i].Data()), legDraw[i]);
		}
		if(showLabel[rowidx][pidx]) model_info->DrawClone("same");
		if(showWatermark[rowidx][pidx]) latexW->DrawLatex(labx, laby,label);
		if(showLeg[rowidx][pidx])legspectra2->DrawClone("same");
		//if(showLabel[rowidx][pidx])latexP->DrawLatex(panx[ridx], pany,descPanel);
		
		str=Form("%s/dpT_%s_R0%.0lf_pTlead%.0lf_probes.%s",out_dir.Data(),name.Data(),Rpar*10,pTlead,ext.Data());
		cdelta2->SaveAs(str.Data());
		
		//ratios
		TH1D* hprobe_ratio=(TH1D*) hdpt2[0]->Clone("hprobe_ratio");
		hprobe_ratio->Divide(hdpt2[2]);
		TCanvas *crat = new TCanvas("crat","probe ratio",10,10,can_x,can_y);
		crat->cd();
		hprobe_ratio->SetTitle(Form("R = %.1lf",Rpar));
		hprobe_ratio->GetXaxis()->SetTitle("#delta p_{T} (GeV/#it{c})");
		hprobe_ratio->GetYaxis()->SetTitle("SP/Pythia");
		hprobe_ratio->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
		hprobe_ratio->GetYaxis()->SetRangeUser(0.1, 2);
		hprobe_ratio->Draw();
		
		TLine *one = new TLine(spectraXmin, 1, spectraXmax, 1);
		one->SetLineWidth(2);
		one->SetLineStyle(2);
		one->SetLineColor(kGray);
		one->DrawClone("same");
		str=Form("%s/dpT_%s_R0%.0lf_pTlead%.0lf_probe_ratio.%s",out_dir.Data(),name.Data(),Rpar*10,pTlead,ext.Data());
		crat->SaveAs(str.Data());
		
}

//***********************************************
//Response matrix
//***********************************************
plot_RM(float R=0.4, float pTlead=5.0, bool doToymodel=0, systems system=cent, TString version="GPC2",TString label="", TString ext="gif")
{
	
	//histo sizes
	const int can_x=500; //1600
	const int can_y=500; //900
	const  float rmXmin=-30;
	const float rmXmax=50;
	const float rmYmin=0;
	const float rmYmax=50;
	const float projectYmin=0;
	const float projectYmax=1.5;
	const float projectYmin2[]={0,0.001,0.1};
	const float projectYmax2[]={1.5,1E4,1E4};
	
	gStyle->SetPadRightMargin(0.15);
	
	int prior=5; //prior distribution (5=pT^-5, 2=biased pythia)
	TString prior_name[]={"",/*1*/"flat",/*2*/"pythia",/*3*/"powlaw3",/*4*/"powlaw45",/*5 */"powlaw5",/*6*/"powlaw55"};
	
	TString out_dir=Form("../../plotting_out/obr/%s/RM",version.Data());
	TString data_names[]={"MB","toymodel"};
	if(system==peri)data_names[0]+="_peripheral";
	else if(system==cent)data_names[0]+="_central";
	
	const int nrm=3; //numbre of RM types (3: dpT, dete, dpT x dete)
	const int nrm2=3; //numbre of RM-BGxD types (3: normalized, weighted by prior, weighted by prior + coarse binning)
	TString rm_ltype[nrm]={"#delta p_{T} RM", "detector RM", "#delta p_{T} RM #times detector RM"};
	TString rm_stype[nrm]={"BG", "D", "BGD"};
	TString rm_ltype2[nrm2]={"probability", "# of events", "# of ev., coarse bins"};
	TString rm_stype2[nrm2]={"Norm", "Weighted", "Coarse"};
	
	//in which panels do we want to show info, legend and labels?
	bool showInfo[2][nrm]={/*RM :*/{1,0,0},
	/*Yprojection*/{1,1,1}};
	bool showWatermark[2][nrm]={/*RM :*/{0,0,0},
	/*Yprojection*/{0,0,0}};
	
	TFile *fin[nrm];
	TFile *fin2[nrm2];
	TString dname=Form("%s/%s",data_names[doToymodel].Data(),version.Data());
	TString str = Form("../../plotting_out/intersteps/RM/%s/response_matrix_BG_sp_R%.1lf_pTlead%.0lf.root",dname.Data(),R,pTlead);
	fin[0] = new TFile(str.Data(), "OPEN");
	str=Form("../../plotting_out/intersteps/RM/%s/pythia_emb_R%.1lf.root",dname.Data(),R,pTlead);
	fin[1] = new TFile(str.Data(), "OPEN");
	str=Form("../../plotting_out/intersteps/RM/%s/response_matrix_BGD_R%.1lf_pTlead%.0lf.root",dname.Data(),R,pTlead);
	fin[2] = new TFile(str.Data(), "OPEN");
	
	fin2[0] = fin[2]; //hrm2[0]=hrm[2]
	str=Form("../../plotting_out/intersteps/RM/%s/response_matrix_BGD_R%.1lf_pTlead%.1lf_allPriors.root",dname.Data(),R,pTlead);
	fin2[1] = new TFile(str.Data(), "OPEN");
	str=Form("../../plotting_out/root/%s/Unfolded_R%.1lf_Bayes_VARbins_bining1_BGD_%s_normal/%s/unfolded_SignalSpectrum_R%.1lf_pTthresh%.1lf.root",dname.Data(),R,version.Data(),prior_name[prior].Data(),R,pTlead);
	fin2[2] = new TFile(str.Data(), "OPEN");
	TDirectoryFile* dir = (TDirectoryFile*)fin2[2]->Get("input");
	
	TH2D* hrm[nrm];
	TH1D* hrm_projectY[nrm];
	for(int i=0;i<nrm;i++)
	{
		if(i==1)hrm[i]=(TH2D*)fin[i]->Get(Form("hresponse_pTl%.0lf",pTlead));
		else hrm[i]=(TH2D*)fin[i]->Get("hResponse_1E9");
		TString hname=Form("hproject_%i",i);
		hrm_projectY[i]=(TH1D*) hrm[i]->ProjectionY(hname,1,hrm[i]->GetXaxis()->GetNbins());
	}
	
	TH2D* hrm2[nrm];
	TH1D* hrm_projectY2[nrm];
	for(int i=0;i<nrm2;i++)
	{
		if(i==1) hrm2[i]=(TH2D*)fin2[i]->Get(Form("hResponse_%s",prior_name[prior].Data()));
		else if(i==2) hrm2[i]=(TH2D*)dir->Get("hresponse");
		else hrm2[i]=(TH2D*)fin2[i]->Get("hResponse_1E9");
		TString hname=Form("hproject2_%i",i);
		hrm_projectY2[i]=(TH1D*) hrm2[i]->ProjectionY(hname,1,hrm2[i]->GetXaxis()->GetNbins());
	}
	
	TLatex *latex = new TLatex();
	latex->SetNDC();
	latex->SetTextSize(0.035);
	float desx=0.4;
	float desy=0.93;
	
	//watermark
	TLatex *latexW = new TLatex();
	latexW->SetNDC();
	latexW->SetTextSize(G_latex_sz);
	latexW->SetTextColor(G_latex_cl);
	float labx=0.65;
	float laby=0.2;
	
	//crm->Divide(nrm,1);
	
	TLegend *model_info = new TLegend(0.20, 0.68, 0.35, 0.88);
	model_info->SetTextSize(0.04);
	model_info->SetFillStyle(0);
	model_info->SetBorderSize(0);
	model_info->SetMargin(0.05);
	model_info->SetHeader("Au+Au #sqrt{s_{NN}} = 200 GeV");
	if(doToymodel)model_info->SetHeader("Parametrized model");
	
	TH1 *frame = new TH1I("frame", "", 1000, -100, +100);
	frame->GetXaxis()->SetTitle("p_{T, jet}^{det}  (GeV/#it{c})");
	frame->GetYaxis()->SetTitle("p_{T, jet}^{part}  (GeV/#it{c})");
	frame->GetZaxis()->SetTitle("probability");
	frame->GetXaxis()->SetRangeUser(rmXmin, rmXmax);
	frame->GetYaxis()->SetRangeUser(rmYmin, rmYmax);
	//frame->GetYaxis()->SetRangeUser(1E-3, 1E5);
	TString str="";
	frame->SetTitle(str);
	
	//-------------------------------------
	//plot 3 RM - dpT, dete, dpt x dete
	//-------------------------------------
	TCanvas *crm[nrm];
	
	for(int i=0; i<nrm; i++)
	{
		crm[i]= new TCanvas(Form("crm_%i",i),"response matrix",10,10,can_x,can_y);
		crm[i]->cd();
		gPad->SetLogz();
		frame->SetTitle(rm_ltype[i]);
		frame->DrawCopy("");
        hrm[i]->GetZaxis()->SetTitle("");
		hrm[i]->DrawCopy("COLZsame");
		//latex->DrawLatex(desx, desy,rm_ltype[i]);
		if(i==0) 
		{
			model_info->AddEntry("", G_system_info[system],"");
			model_info->AddEntry("", Form("anti-k_{T}, R = %.1lf",R), "");
			if(pTlead>0) model_info->AddEntry("", Form("p_{T,lead}^{min} = %.0lf GeV/#it{c}",pTlead), "");
		}
		if(showInfo[0][i])model_info->DrawClone("same");
		if(showWatermark[0][i])latexW->DrawLatex(labx, laby,label);
		
		str=Form("%s/rm%s_%s_R0%.0lf_pTlead%.0lf.%s",out_dir.Data(),rm_stype[i].Data(),data_names[doToymodel].Data(),R*10,pTlead,ext.Data());
		crm[i]->SaveAs(str.Data());
	}//RM type loop
	
	//-------------------------------------
	//plot 3 RM BGxD - normalize, prior scaled, prior scaled + coarse binning
	//-------------------------------------
	TCanvas *crm2[nrm2];
	
	for(int i=0; i<nrm2; i++)
	{
		crm2[i]= new TCanvas(Form("crm2_%i",i),"response matrix",10,10,can_x,can_y);
		crm2[i]->cd();
		if(i!=1)gPad->SetLogz();
		frame->SetTitle(rm_ltype2[i]);
		frame->DrawCopy("");
		hrm2[i]->DrawCopy("COLZsame");
		//latex->DrawLatex(desx, desy,rm_ltype2[i].Data());
		
		if(showInfo[0][i])model_info->DrawClone("same");
		if(showWatermark[0][i])latexW->DrawLatex(labx, laby,label);
		
		str=Form("%s/rm%s_%s_R0%.0lf_pTlead%.0lf.%s",out_dir.Data(),rm_stype2[i].Data(),data_names[doToymodel].Data(),R*10,pTlead,ext.Data());
		crm2[i]->SaveAs(str.Data());
	}//RM type loop
	
	
	//----------------------------------------------------
	//plot Y projjections of RMs (dpT, dete, dpt x dete)
	//----------------------------------------------------	
	TString projectYtitleY1="prob. dens.";
	TString projectYtitleY2="arb. units";
	
	TCanvas *crmp[nrm];
	//crmp->Divide(nrm,1);
	
	frame->GetXaxis()->SetTitle("p_{T, jet}^{generated}  (GeV/#it{c})");
	frame->GetYaxis()->SetTitle(projectYtitleY1);
	frame->GetXaxis()->SetRangeUser(rmYmin, rmYmax);
	frame->GetYaxis()->SetRangeUser(projectYmin, projectYmax);
	
	for(int i=0; i<nrm; i++)
	{
		crmp[i] = new TCanvas(Form("crmp_%i",i),"response matrix Y projection",10,10,can_x,can_y);
		crmp[i]->cd();
		frame->SetTitle(rm_ltype[i]);
		frame->DrawCopy();
		hrm_projectY[i]->SetLineWidth(1);
		hrm_projectY[i]->SetLineColor(kBlue+1);
		hrm_projectY[i]->DrawCopy("same");
		
		//latex->DrawLatex(desx, desy,rm_ltype[i]);
		
		if(showInfo[1][i])model_info->DrawClone("same");
		if(showWatermark[1][i])latexW->DrawLatex(labx, laby,label);
		
		str=Form("%s/rm%s_Yp_%s_R0%.0lf_pTlead%.0lf.%s",out_dir.Data(),rm_stype[i].Data(), data_names[doToymodel].Data(),R*10,pTlead,ext.Data());
		crmp[i]->SaveAs(str.Data());
	}//RM type loop
	
	//----------------------------------------------------
	//plot Y projjections of RMs (normalize, prior scaled, prior scaled + coarse binning)
	//----------------------------------------------------
	TCanvas *crmp2[nrm];
	//crmp->Divide(nrm,1);
	
	for(int i=0; i<nrm2; i++)
	{
		crmp2[i] = new TCanvas(Form("crmp2_%i",i),"response matrix Y projection",10,10,can_x,can_y);
		crmp2[i]->cd();
		if(i>0) crmp2[i]->SetLogy();
		frame->GetYaxis()->SetTitle(projectYtitleY2);
		frame->GetYaxis()->SetRangeUser(projectYmin2[i], projectYmax2[i]);
		frame->SetTitle(rm_ltype2[i]);
		frame->DrawCopy();
		hrm_projectY2[i]->SetLineWidth(1);
		hrm_projectY2[i]->SetLineColor(kBlue+1);
		hrm_projectY2[i]->DrawCopy("same");
		//latex->DrawLatex(desx, desy,rm_ltype2[i].Data());
		
		if(showInfo[1][i])model_info->DrawClone("same");
		if(showWatermark[1][i])latexW->DrawLatex(labx, laby,label);
		
		str=Form("%s/rm%s_Yp_%s_R0%.0lf_pTlead%.0lf.%s",out_dir.Data(),rm_stype2[i].Data(), data_names[doToymodel].Data(),R*10,pTlead,ext.Data());
		crmp2[i]->SaveAs(str.Data());
	}//RM type loop
	
}


//================================================================
//plot number of trucks as a function of BBC rate (pileup study)
//================================================================
void show_bbc_ntr(int cut=0, TString ext="gif")
{
	//TString suff[3]={"nocut","ranking","bemctof"};
	TString suff[3]={"nocut","badrank","bemc"};
	//TString desc[3]={"no cuts","|z(vpd)-z(TPC)|<4cm, DCA<1cm, ranking>0","VPD cut, DCA cut, TOF or BEMC match"};
	TString desc[3]={"no cuts","|z(vpd)-z(TPC)|<4cm, DCA<1cm, ranking<0","VPD cut, DCA cut, BEMC match"};
	TFile* fin=new TFile(Form("../../plotting_out/intersteps/pileup/qa_%s.root",suff[cut].Data()),"OPEN");
	TString out_path=("../../plotting_out/obr/intersteps/pileup");
	
	TH2D* h2=(TH2D*)fin->Get("hbbcrate_ntr");
	TH1D* h1=h2->ProjectionX("h1",1,h2->GetNbinsY());
	h1->Reset("MICE");
	h1->SetTitle("");
	h1->GetYaxis()->SetTitle("<N> tracks");
	
	//make also 3 example slices of #of tracks	
	float bbc_example_slice[3]={500000,1500000,2500000}; //bbc values for which we show the example slices
	TH1D* hntr_example[3];
	Color_t color[3]={kRed,kBlue,kBlack};
	int binex[3];
	for(int j=0;j<3;j++)
	{
		binex[j]=h2->GetXaxis()->FindBin(bbc_example_slice[j]);
	}
	for(int bn=1; bn<h2->GetNbinsX();bn++)
	{
		TH1D* hslice=(TH1D*) h2->ProjectionY(Form("hslice_%i",bn),bn,bn);
		float mean=hslice->GetMean();
		h1->SetBinContent(bn,mean);
		
		//save example slices
		for(int j=0;j<3;j++)
		{
			if(binex[j]==bn) 
			{
				hntr_example[j]=(TH1D*) hslice->Clone(Form("hexample_%i",j));
				hntr_example[j]->SetLineColor(color[j]);
				hntr_example[j]->Scale(1.0/hntr_example[j]->Integral());
				hntr_example[j]->SetLineWidth(2);
			}
		}
		delete hslice;
		
	}	
	
	
	h1->Rebin(5);
	h1->Scale(1.0/3.0);
	TCanvas* c1=new TCanvas("c1","<N> tracks vs. BBC rate",10,10,900,600);
	c1->cd();
	int binmax = h1->GetMaximumBin();
	float Ymax=1.3*h1->GetBinContent(binmax);
	h1->GetYaxis()->SetRangeUser(0, Ymax);
	h1->Draw("");
	
	TLatex *latex = new TLatex();
	latex->SetNDC();
	latex->SetTextSize(0.045);
	latex->DrawLatex(0.2, 0.8,desc[cut]);
	
	if(cut==1)
	{
		TLine *one = new TLine(0, 8.5, 500000, 8.5);
		one->SetLineWidth(2);
		one->SetLineStyle(2);
		one->SetLineColor(kGray);
		one->DrawClone("same");
		
		TLine *two = new TLine(0, 10.25, 2500000, 10.25);
		two->SetLineWidth(2);
		two->SetLineStyle(2);
		two->SetLineColor(kGray);
		two->DrawClone("same");
		
	}	
	
	c1->SaveAs(Form("%s/bbcrate_vs_ntr_%s.%s",out_path.Data(),suff[cut].Data(),ext.Data()));
	
	TCanvas* c2=new TCanvas("c2","<N> tracks examples",10,10,900,600);
	c2->cd();
	c2->SetLogy();
	hntr_example[0]->SetTitle("");
	hntr_example[0]->Draw("");
	for(int j=1;j<3;j++)
		hntr_example[j]->Draw("same");
	
	TLegend *legraa = new TLegend(0.2,0.2,0.7,0.35);
	legraa->SetTextSize(0.035);
	legraa->SetFillStyle(0);
	legraa->SetBorderSize(0);	
	for(int j=0;j<3;j++)
	{
		legraa->AddEntry(hntr_example[j],Form("BBC rate: %.0lf",bbc_example_slice[j]),"l");
	}
	legraa->Draw("same");
	
	float xtx[3]={0.7,0.5,0.4};
	latex->DrawLatex(xtx[cut], 0.8,desc[cut]);
	
	c2->SaveAs(Form("%s/ntr_%s.%s",out_path.Data(),suff[cut].Data(),ext.Data()));
}

//================================================================
//plot toymodel background distribution corrected for the 
//================================================================
void show_toy_BG(float R=0.3, float pTlead=5.0,int skipNbins=0,TString suffix="")
{
	//load files and histograms
	//measured spectrum
	TString indir="../../plotting_out/intersteps/toy_background";
	TFile* fspec=new TFile(Form("%s/histos_jets_R%.1lf_pTcut0.2%s.root",indir.Data(),R,suffix.Data()),"OPEN");
	TH1I* hevents=(TH1I*) fspec->Get("hevents");
	int nevents=hevents->GetBinContent(2);
	float scaler=1.0/(2*(1.0-R)*2*TMath::Pi()*nevents);
	TH2D* hspec2D=(TH2D*) fspec->Get("fhDSpTleading");
	int firstbin=hspec2D->GetXaxis()->FindBin(pTlead);
	int lastbin=hspec2D->GetXaxis()->GetNbins();
	TH1D* hspectrum=(TH1D*) hspec2D->ProjectionY("hspectrum",firstbin,lastbin);
	hspectrum->Scale(scaler,"width");
	
	//efficiency histogram
	TFile* feff=new TFile(Form("%s/epsilon_R%.1lf_pTlead%.0lf.root",indir.Data(),R,pTlead),"OPEN");
	TH1D* hepsilon=(TH1D*)feff->Get("hepsilon_unfolded");
	
	//correct for efficiency
	int nbins=hspectrum->GetNbinsX();
	for(int bin=1; bin<nbins;bin++)
	{
		float val=hspectrum->GetBinContent(bin);
		float pT=hspectrum->GetBinCenter(bin);
		
		int effbin=hepsilon->FindBin(pT);
		float eff=hepsilon->GetBinContent(effbin);
		
		float corr=(eff>0) ? (1.0/eff) : 0;
		hspectrum->SetBinContent(bin,val*corr);
	}
	
	TCanvas* c1=new TCanvas("c1","background",10,10,800,600);
	c1->cd();
	c1->SetLogy();
	hspectrum->SetTitle("jet reco. efficiency corrected BG");
	hspectrum->Draw();
	
	//calculate yield
	double integral=0;
	int firstbin=hspectrum->FindBin(pTlead);
	for(int bn=firstbin+skipNbins; bn<=nbins; bn++)
	{
		double y=hspectrum->GetBinContent(bn);
		double width=hspectrum->GetBinWidth(bn);
		integral+=width*y;
	}
	
	cout<<"Yield (x1E6): "<<integral*1E6<<endl;
	
	
}


//========================================================================
//additional utility functions
//========================================================================

double hardjet(double *x, double *par)
{
	double pT = *x;
	
	double C = 1/(1 + TMath::Exp(-(pT - 2.8)/0.3)); //cutoff function
	double J,y1,y2, B, T, n, m0, mu,pT0,R,a,b,A,pwr;
	B    = par[1];
	T    = par[2];
	n    = par[3];
	m0   = par[4];
	mu   = par[5];
	A=par[6];
	pwr=par[7];
	double mT = TMath::Sqrt((pT-mu)*(pT-mu)+m0*m0);
	
	y2 = B/TMath::Power(1.0+(mT-m0)/(n*T),n);
	y1=A*TMath::Power(pT,-pwr);
	
	if(pT<10)
		J=y1;
	else if(pT<20)
	{
		double c1=(pT-10.0)/10.0;
		J=((1-c1)*y1+c1*y2);
	}
	else J=y2;
	if(J<0)J=0;
	
	//pT-dependent RAA
	double RAA;
	/*
	 * double RAA_lowpT=0.5;
	 * if(par[8]>RAA_lowpT)RAA_lowpT=par[8];
	 * if(pT<5)RAA=RAA_lowpT;
	 * //else if(pT<10)RAA=1.0;
	 * else if(pT<15)RAA=(RAA_lowpT-par[8])*(15-pT)/10+par[8];
	 * else RAA=par[8];*/
	RAA=par[8];
	
	return par[0]*J*RAA/*C*/; //hard jet spectrum (with suppressed low pT part)
}

int calc_R_index(float R)
{
	int ridx=0;
	if(R>0.25) ridx=1;
	if(R>0.35) ridx=2;
	if(R>0.45) ridx=3;
	return ridx;
}
