//this macro calculates JET ENERGY SCALE (error in jet pT (x-axis) arrising from error in jet yield (y-axis))
//and pT shift between central collsions and peripheral collisins jet spectra (an alternative to RAA)

#include "../utils/common.h" //several constants
#include "../utils/utils.C" //some useful functions

enum ratios { RAA, RCP };

void calc_pTshift(int pTlead=5, int bining=1, TString evo="omicron", bool print_JES=1, bool print_pTshift=1)
{
	bool useTAAerror=1; //include TAA errors in the RCP systematic uncertainty
	//colors and markes fo pTshift graphs
	Color_t cjan=kBlue;
	Color_t cjan2=kGray;
	Color_t calex=kRed;
	int mjan=29;
	int malex=21;

	//between which values do we want to evaluate?
	float xmin=15;
	float xmax=30;
	float xmax_peripheral=25;
	
	//y range of plotted graphs
	float ymin=0;
	float ymax=8;
	
	const float NBIN[/*systems*/]={955.4,20.4,42.0/30.0}; //average number of binary collisions

	TCanvas* c1=new TCanvas("c1","spectra",10,10,1200,600);
	c1->Divide(C_nR,1);
	
	TCanvas* c2=new TCanvas("c2","shift",10,10,1600,600);
	c2->Divide(C_nR,1);
	
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.04);
	
	//declar arrays
	TFile* fin1[C_nR];
	TFile* fin2[C_nR];
	Graph* gstat1[C_nR];
	Graph* gunf1[C_nR];
	Graph* gcorr1[C_nR];
	TF1* fspec1[C_nR];	
	Graph* gstat2[C_nR];
	Graph* gunf2[C_nR];
	Graph* gcorr2[C_nR];
	TF1* fspec2[C_nR];	
	Graph* gstat_dpT[C_nR];
	Graph* gunf_dpT[C_nR];
	Graph* gcorr_dpT[C_nR];
	Graph* gsys1[C_nR];
	Graph* gsys2[C_nR];
	Graph* gsys_dpT[C_nR];
	TGraph *grinter[C_nR];
	TH1D* hspec2[C_nR];
	
	//pT shift graphs
	Graph* gshift[C_nR];
	Graph* gshift_sys[C_nR];
	Graph* gshift_avg[C_nR];
	Graph* gshift_avg_sys[C_nR];
	//Alex's pT shift for recoil jets
	TGraphErrors* gshift_alex[C_nR];
	TGraphErrors* gshift_alex_sys[C_nR];
	double sh_alex_x[1]={15};
	double sh_alex_xerr[1]={5};
	double sh_alex_r02[1]={4.4};
	double sh_alex_r03[1]={5.0};
	double sh_alex_r04[1]={5.1};
	double sh_alex_r02_stat[1]={0.2};
	double sh_alex_r03_stat[1]={0.5};
	double sh_alex_r04_stat[1]={0.5};
	double sh_alex_r02_sys[1]={1.2};
	double sh_alex_r03_sys[1]={1.2};
	double sh_alex_r04_sys[1]={1.2};
	gshift_alex[0]=new TGraphErrors(1,sh_alex_x,sh_alex_r02,sh_alex_xerr,sh_alex_r02_stat);
	gshift_alex_sys[0]=new TGraphErrors(1,sh_alex_x,sh_alex_r02,sh_alex_xerr,sh_alex_r02_sys);
	gshift_alex[1]=new TGraphErrors(1,sh_alex_x,sh_alex_r03,sh_alex_xerr,sh_alex_r03_stat);
	gshift_alex_sys[1]=new TGraphErrors(1,sh_alex_x,sh_alex_r03,sh_alex_xerr,sh_alex_r03_sys);
	gshift_alex[2]=new TGraphErrors(1,sh_alex_x,sh_alex_r04,sh_alex_xerr,sh_alex_r04_stat);
	gshift_alex_sys[2]=new TGraphErrors(1,sh_alex_x,sh_alex_r04,sh_alex_xerr,sh_alex_r04_sys);
	
	for(int r=0; r<C_nR; r++)
	{
	float R=C_Rs[r];
	
	fin1[r]=open_sys_file(cent, R, pTlead, bining, evo, 0);
	fin2[r]=open_sys_file(peri, R, pTlead, bining, evo, 0);
	
	//input for JES calculation
	gstat1[r]=(Graph*) fin1[r]->Get("spectrum_BSL");
	/*Graph*/ gunf1[r]=(Graph*) fin1[r]->Get("spectrum_unfold");
	/*Graph*/ gcorr1[r]=(Graph*) fin1[r]->Get("spectrum_corr");
	fspec1[r]=(TF1*) fin1[r]->Get("spectrum_fit");	

	gstat2[r]=(Graph*) fin2[r]->Get("spectrum_BSL");
	/*Graph*/ gunf2[r]=(Graph*) fin2[r]->Get("spectrum_unfold"); //we take errors from RCP for the pT-shift error calculation since they take into account all the correlations between central end peripheral col.
	/*Graph*/ gcorr2[r]=(Graph*) fin2[r]->Get("spectrum_corr"); //we take errors from RCP for the pT-shift error calculation since they take into account all the correlations between central end peripheral col.
	fspec2[r]=(TF1*) fin2[r]->Get("spectrum_fit");	

	//input for pT-shift calculation
	gstat_dpT[r]=(Graph*) fin1[r]->Get("RCP_BSL");
	/*Graph*/ gunf_dpT[r]=(Graph*) fin1[r]->Get("RCP_unfold");
	/*Graph*/ gcorr_dpT[r]=(Graph*) fin1[r]->Get("RCP_corr");
	
	//create graphs with combined errors
	gsys1[r]=CombineGraphErrors((Graph*) gunf1[r],(Graph*) gcorr1[r]);
	gsys2[r]=CombineGraphErrors((Graph*) gunf2[r],(Graph*) gcorr2[r]);
	gsys_dpT[r]=CombineGraphErrors((Graph*) gunf_dpT[r],(Graph*) gcorr_dpT[r]);
/*
	delete gunf1;
	delete gunf2;
	delete gcorr1;
	delete gcorr2;
*/
		
	float scale_RCP=NBIN[0]/NBIN[1];


	bool analyzed=false;
	int npoints=gstat1[r]->GetN();
	
	//calculate also average values over all bins
	double delta_avg=0;
	double delta_avg_sys_up=0;
	double delta_avg_sys_down=0;
	double delta_avg_stat_up=0;
	double delta_avg_stat_down=0;
	int N=0;
	double JESc_up=0;
	double JESc_down=0;
	double JESp_up=0;
	double JESp_down=0;
	int NJESc=0;
	int NJESp=0;
	
	
	//arrays for graph
	const int npts=npoints;
	double sh_x[npts];
	double sh_xerr_l[npts];
	double sh_xerr_h[npts];
	double sh_y[npts];
	double sh_stat_up[npts];
	double sh_sys_up[npts];
	double sh_stat_down[npts];
	double sh_sys_down[npts];
	double sh_avg_x[1];
	double sh_avg_xerr_l[1];
	double sh_avg_xerr_h[1];
	double sh_avg_y[1];
	double sh_avg_stat_up[1];
	double sh_avg_sys_up[1];
	double sh_avg_stat_down[1];
	double sh_avg_sys_down[1];
	
	for(int pt=0; pt<npoints; pt++)
	{
		//==========================
		double xval1;
		double yval;
		gstat1[r]->GetPoint(pt,xval1,yval);
		//statistical error
		double errStat1_up=gstat1[r]->GetErrorYhigh(pt);
		double errStat1_down=gstat1[r]->GetErrorYlow(pt);
		//systematic error
		double errSys1_up=gsys1[r]->GetErrorYhigh(pt);
		double errSys1_down=gsys1[r]->GetErrorYlow(pt);
		if(useTAAerror)
		{
			errSys1_up=TMath::Sqrt(errSys1_up*errSys1_up+(yval*C_RCP_NormErr)*(yval*C_RCP_NormErr));
			errSys1_down=TMath::Sqrt(errSys1_down*errSys1_down+(yval*C_RCP_NormErr)*(yval*C_RCP_NormErr));
		}
		//for peripheral collisins
		double xval2;
		double yval2;
		gstat2[r]->GetPoint(pt,xval2,yval2);
		double errSys2_up=gsys2[r]->GetErrorYhigh(pt);
		double errSys2_down=gsys2[r]->GetErrorYlow(pt);
		
		if(print_JES || print_pTshift) cout<<"analyzing point "<<pt<<", x="<<xval1<<endl;
		
		//move also the yvalue within the errors
		double yval_stat_up=yval+errStat1_up;
		double yval_stat_down=yval-errStat1_down;
		double yval_sys_up=yval+errSys1_up;
		double yval_sys_down=yval-errSys1_down;
		//peripheral
		double yval2_sys_up=yval2+errSys2_up;
		double yval2_sys_down=yval2-errSys2_down;
		
		//for pT shift calculation: rescale y values by corresponding <NBIN> ratio so the central collisins can be compared to peripheral
		float yval_scl=yval/scale_RCP;
		float yval_stat_up_scl=yval_stat_up/scale_RCP;
		float yval_stat_down_scl=yval_stat_down/scale_RCP;
		float yval_sys_up_scl=yval_sys_up/scale_RCP;
		float yval_sys_down_scl=yval_sys_down/scale_RCP;

		//*********************************
		//JET ENRGY SCALE UNCERTAINTY
		//*********************************

		if(xval1>xmin && xval1<xmax )
		{
		if(print_JES)cout<<"JES - central"<<endl;
		//cout<<"x, y, y+y_sys_up, y-y_sys_down:"<<xval1<<" | "<<yval<<" | "<<yval_sys_up<<" | "<<yval_sys_down<<endl;
		//!!!Falling spectrum => y_err_UP ~ x_err_DOWN
		float x_err_down=interpolate_x(gstat1[r], yval_sys_up, none, NULL,1,0);
		float x_err_up=interpolate_x(gstat1[r], yval_sys_down, none, NULL,1,0);
		//this is a crosscheck
		float x_err_upB=interpolate_x(gstat1[r], yval, up, gsys1[r],1,0);
		float x_err_downB=interpolate_x(gstat1[r], yval, down, gsys1[r],1,0);
		if(print_JES)cout<<"pT:"<<xval1<<", sys. err. up [GeV]:+"<<x_err_up-xval1<<", down [GeV]:-"<<xval1-x_err_down<<endl;
		//cout<<"x-check: sys. err. up (B) [GeV]:+"<<x_err_upB-xval1<<", down (B) [GeV]:-"<<xval1-x_err_downB<<endl;
		if(print_JES)cout<<"JES [%]: +"<<(x_err_up-xval1)*100/xval1<<", -"<<(xval1-x_err_down)*100/xval1<<endl;
		if(print_JES)cout<<"JES (B) [%]: +"<<(x_err_upB-xval1)*100/xval1<<", -"<<(xval1-x_err_downB)*100/xval1<<endl;
		if(print_JES)cout<<"==================================================="<<endl;
		//calculate average
		JESc_up+=((x_err_up+x_err_upB)/2-xval1)/xval1;
		JESc_down+=(xval1-(x_err_down+x_err_downB)/2)/xval1;
		NJESc++;
		}
		
		if(xval2>xmin && xval2<xmax_peripheral)
		{
		if(print_JES)cout<<"JES - peripheral"<<endl;
		//cout<<"x, y, y+y_sys_up, y-y_sys_down:"<<xval1<<" | "<<yval<<" | "<<yval_sys_up<<" | "<<yval_sys_down<<endl;
		//!!!Falling spectrum => y_err_UP ~ x_err_DOWN
		float x_err_down2=interpolate_x(gstat2[r], yval2_sys_up, none, NULL,1,0);
		float x_err_up2=interpolate_x(gstat2[r], yval2_sys_down, none, NULL,1,0);
		//this is a crosscheck
		float x_err_up2B=interpolate_x(gstat2[r], yval2, up, gsys2[r],1,0);
		float x_err_down2B=interpolate_x(gstat2[r], yval2, down, gsys2[r],1,0);
		if(print_JES)cout<<"pT:"<<xval2<<", sys. err. up [GeV]:+"<<x_err_up2-xval2<<", down [GeV]:-"<<xval2-x_err_down2<<endl;
		//cout<<"x-check: sys. err. up (B) [GeV]:+"<<x_err_up2B-xval2<<", down (B) [GeV]:-"<<xval2-x_err_down2B<<endl;
		if(print_JES)cout<<"JES [%]: +"<<(x_err_up2-xval2)*100/xval2<<", -"<<(xval2-x_err_down2)*100/xval2<<endl;
		if(print_JES)cout<<"JES (B) [%]: +"<<(x_err_up2B-xval2)*100/xval2<<", -"<<(xval2-x_err_down2B)*100/xval2<<endl;
		if(print_JES)cout<<"==================================================="<<endl;
		
		JESp_up+=((x_err_up2+x_err_up2B)/2-xval2)/xval2;
		JESp_down+=(xval2-(x_err_down2+x_err_down2B)/2)/xval2;
		NJESp++;
		}
		

		if(xval1<xmin || xval1>xmax_peripheral)continue;

		//*********************************
		//pT-shift central-peripheral
		//*********************************
		
		//A) from the fit function
		//==========================

		float xval2_fit=fspec2[r]->GetX(yval_scl);
		float dx_fit=xval1-xval2_fit;
		
		float yval2_fit=fspec2[r]->Eval(xval1);
		float rcp=yval_scl/yval2_fit;
		if(print_pTshift) cout<<"RCP:"<<rcp<<endl;
			
		if(print_pTshift) cout<<"pT (central):"<<xval1<<" delta pT (from fit):"<<dx_fit<<endl;


		//B) linear interpolation
		//==========================
		
		double xval2_lin=interpolate_x(gstat2[r], yval_scl, none);
		double xval2_stat_up=interpolate_x(gstat2[r], yval_scl/*yval_stat_up_scl*/, up, gstat_dpT[r],1,1);
		double xval2_stat_down=interpolate_x(gstat2[r], yval_scl/*yval_stat_down_scl*/, down, gstat_dpT[r],1,1);
		double xval2_sys_up=interpolate_x(gstat2[r], yval_scl/*yval_sys_up_scl*/, up, gsys_dpT[r],1,1);
		double xval2_sys_down=interpolate_x(gstat2[r], yval_scl/*yval_sys_down_scl*/, down, gsys_dpT[r],1,1);
		double dx_lin=xval1-xval2_lin;
		
		//calculate total error as maximum of left and right error
		double err_stat1=TMath::Abs(xval2_lin-xval2_stat_up);
		double err_stat2=TMath::Abs(xval2_lin-xval2_stat_down);
		double err_stat=(err_stat1>err_stat2) ? err_stat1 : err_stat2;
		double err_sys1=TMath::Abs(xval2_lin-xval2_sys_up);
		double err_sys2=TMath::Abs(xval2_lin-xval2_sys_down);
		double err_sys=(err_sys1>err_sys2) ? err_sys1 : err_sys2;
		
		
		if(print_pTshift)
		{
			cout<<"==================================================="<<endl;
			cout<<"interpolated point (x,y):"<<endl;
			cout<<"("<<xval2_lin<<","<<yval_scl<<")"<<endl;
		
			cout<<"pT (central):"<<xval1<<" delta pT (from interpolation):"<<dx_lin<<endl;
			cout<<"pT (central):"<<xval1<<" stat. err.:"<<err_stat<<endl;
			cout<<"pT (central):"<<xval1<<" sys. err.:"<<err_sys<<endl;
			cout<<"==================================================="<<endl;
		}
		
	
		delta_avg+=dx_lin;
		delta_avg_stat_down+=(err_stat);
		delta_avg_stat_up+=(err_stat);
		delta_avg_sys_down+=(err_sys);
		delta_avg_sys_up+=(err_sys);
		
		
		sh_x[N]=xval1;
		sh_xerr_l[N]=gstat1[r]->GetErrorXlow(pt);
		sh_xerr_h[N]=gstat1[r]->GetErrorXhigh(pt);;
		sh_y[N]=-dx_lin;
		sh_stat_up[N]=err_stat;
		sh_stat_down[N]=err_stat;
		sh_sys_up[N]=err_sys;
		sh_sys_down[N]=err_sys;
		
		N++;



	}//loop over points we want to analyze	
		
	//average shift over all bins
	delta_avg=delta_avg/N;
	delta_avg_stat_down=delta_avg_stat_down/N;
	delta_avg_stat_up=delta_avg_stat_up/N;
	delta_avg_sys_down=delta_avg_sys_down/N;
	delta_avg_sys_up=delta_avg_sys_up/N;
	JESc_up=100*JESc_up/NJESc;
	JESc_down=100*JESc_down/NJESc;
	JESp_up=100*JESp_up/NJESp;
	JESp_down=100*JESp_down/NJESp;
	
	sh_avg_x[0]=15;
	sh_avg_xerr_l[0]=5;
	sh_avg_xerr_h[0]=5;
	sh_avg_y[0]=TMath::Abs(delta_avg);
	sh_avg_stat_up[0]=TMath::Abs(delta_avg_stat_up);
	sh_avg_stat_down[0]=TMath::Abs(delta_avg_stat_down);
	sh_avg_sys_up[0]=TMath::Abs(delta_avg_sys_up);
	sh_avg_sys_down[0]=TMath::Abs(delta_avg_sys_down);
	
	//fill graphs
	gshift[r]=new TGraphAsymmErrors(N,sh_x,sh_y,sh_xerr_l,sh_xerr_h,sh_stat_down,sh_stat_up);
	gshift_sys[r]=new TGraphAsymmErrors(N,sh_x,sh_y,sh_xerr_l,sh_xerr_h,sh_sys_down,sh_sys_up);
	tgraph_set_atributes(gshift[r], 5, 25,  ymin, ymax, cjan2, cjan2, cjan2, 2, 1.5, mjan, 0);
	tgraph_set_atributes(gshift_sys[r], 5, 25,  ymin, ymax, cjan2, cjan2, cjan2, 1, 1.5, mjan, 0);
	
	gshift_avg[r]=new TGraphAsymmErrors(1,sh_avg_x,sh_avg_y,sh_avg_xerr_l,sh_avg_xerr_h,sh_avg_stat_down,sh_avg_stat_up);
	gshift_avg_sys[r]=new TGraphAsymmErrors(1,sh_avg_x,sh_avg_y,sh_avg_xerr_l,sh_avg_xerr_h,sh_avg_sys_down,sh_avg_sys_up);
	tgraph_set_atributes(gshift_avg[r], 5, 25,  ymin, ymax, cjan, cjan, cjan, 2, 2.0, mjan, 0);
	tgraph_set_atributes(gshift_avg_sys[r], 5, 25,  ymin, ymax, cjan, cjan, cjan, 1, 2.0, mjan, 0);
	
	tgraph_set_atributes(gshift_alex[r], 5, 25,  ymin, ymax, calex, calex, calex, 2, 1.5, malex, 0);
	tgraph_set_atributes(gshift_alex_sys[r], 5, 25,  ymin, ymax, calex, calex, calex, 1, 1.5, malex, 0);
	
	cout<<""<<endl;
	cout<<"************************************"<<endl;
	cout<<"R="<<R<<" delta="<<delta_avg<<"+-"<<delta_avg_stat_up<<"/"<<delta_avg_stat_down<<"+-"<<delta_avg_sys_up<<"/"<<delta_avg_sys_down<<endl;
	cout<<"JES - central [%]: +"<<JESc_up<<", -"<<JESc_down<<endl;
	cout<<"JES - peripheral [%]: +"<<JESp_up<<", -"<<JESp_down<<endl;
	cout<<"************************************"<<endl;
	cout<<""<<endl;

	c1->cd(r+1);
	gstat1[r]->Draw("AP");
	fspec1[r]->Draw("same");

	//rescale by corresponding <NBIN>
	hspec2[r]=(TH1D*) fspec2[r]->GetHistogram();
	hspec2[r]->Scale(scale_RCP);
	hspec2[r]->SetLineColor(kBlue);
	hspec2[r]->DrawCopy("same");

	//rescale by corresponding <NBIN>
	for (int i=0;i<gstat2[r]->GetN();i++) gstat2[r]->GetY()[i] *= scale_RCP;
	gstat2[r]->Draw("P");
	
	
	c2->cd(r+1);
	gPad->SetMargin(0.13,0.02,0.17,0.02);

	
  TH1D *htmp1=new TH1D("htmp1","",20,xmin-10,xmax+12);
  htmp1->SetXTitle("p_{T,jet}^{ch} (central) [GeV/c]");htmp1->SetYTitle("#Delta p_{T} (central -> peripheral)");
  htmp1->SetTitleOffset(1.2,"x");htmp1->SetTitleOffset(1.0,"y");
  htmp1->SetTitleSize(0.06,"x");htmp1->SetTitleSize(0.05,"y");
  htmp1->SetLabelSize(0.05,"x");htmp1->SetLabelSize(0.04,"y");
  //htmp1->SetNdivisions(505,"x");htmp1->SetNdivisions(505,"y");
  htmp1->SetMinimum(ymin);htmp1->SetMaximum(ymax);
  htmp1->DrawCopy();

	gshift[r]->Draw("P");
	gshift_sys[r]->Draw("2");
	
	gshift_alex[r]->Draw("P");
	gshift_alex_sys[r]->Draw("2");
	
	gshift_avg[r]->Draw("P");
	gshift_avg_sys[r]->Draw("2");
	
	latex->DrawLatex(0.45, 0.85,Form("R=%.1lf",R));
	
	if(r==2)
	{
	double posx1=0.7;
	double posx2=0.95;
	double posy1=0.7;
	double posy2=0.95;
	TLegend *legraa = new TLegend(posx1,posy1,posx2,posy2);
	legraa->SetTextSize(0.032);
	legraa->SetFillStyle(0);
	legraa->SetBorderSize(0);
	legraa->AddEntry(gshift[0], "inclusive jets","lp");
	legraa->AddEntry(gshift_avg[0], "inclusive jets","lp");
	legraa->AddEntry(gshift_avg[0], "  (average)","");
	legraa->AddEntry(gshift_alex[0], "recoil jets","lp");
	legraa->AddEntry(gshift_alex[0], "(w/o p_{T}^{lead} cut)","");
	legraa->Draw("same");
	}
	
	else if(r==0)
	{
		TLegend *model_info = new TLegend(0.15, 0.7, 0.5, 0.95);
		model_info->SetTextSize(0.035);
		model_info->SetFillStyle(0);
		model_info->SetBorderSize(0);
		model_info->SetMargin(0.05);
		model_info->AddEntry("", "STAR Au+Au","");
		model_info->AddEntry("", "#sqrt{s_{NN}}=200GeV","");
		model_info->AddEntry("", "charged jets","");
		model_info->AddEntry("", Form("p_{T}^{lead}>%i GeV/c",pTlead),"");
		model_info->AddEntry("", "anti-k_{T}","");
		model_info->Draw("same");
	}
	

}//loop over R

c2->SaveAs(Form("../../plotting_out/obr/final/central/%s/bining%i/pTshift_pTl%i.gif",evo.Data(),bining,pTlead));
c2->SaveAs(Form("../../plotting_out/obr/final/central/%s/bining%i/pTshift_pTl%i.pdf",evo.Data(),bining,pTlead));
	

	//grinter[r]->Draw("P");

	/*
	float A=fspec2[r]->GetParameter(0);
	A=A*scale_RCP;
	fspec2[r]->SetParameter(0,A);
	fspec2[r]->Update();
	fspec2[r]->SetLineColor(kBlue);
	fspec2[r]->Draw("same");*/
	return;
}

//===========================================================================

TFile* open_sys_file(C_systems system=cent, float R=0.3, float pTlead=5.0,short bining=1, TString evolution="omicron", bool doToy=0)
{
	TString sdir="../../plotting_out/systematics";
	if(doToy)sdir+="/toymodel";
	if(system==peri) sdir+="/peripheral";
	else if(system==pp) sdir+="/pp";
	else sdir+="/central";
	sdir+=Form("/%s",evolution.Data());
	sdir+=Form("/bining%i",bining);
	//sdir+=Form("/%s",test_name[test].Data());
	TString str=Form("%s/systematics_R%.1lf_pT%0.lf.root",sdir.Data(),R,pTlead);
	cout<<"opening file:"<<str.Data()<<endl;
	TFile *fsys;
	fsys= new TFile(str.Data(), "OPEN");
	return fsys;
}

int find_point(TGraphAsymmErrors* graph, float x_value)
{
	int point=0;
	for(int pt=0; pt<graph->GetN(); pt++)
	{
		double xpoint,ypoint,xerr_left,xerr_right;
		graph->GetPoint(pt,xpoint,ypoint);
		xerr_left=graph->GetErrorXlow(pt);
		xerr_right=graph->GetErrorXhigh(pt);
		if(x_value>=xpoint-xerr_left && x_value<=xpoint+xerr_right)
		{
			point=pt;
			break;
		}		
	}	
	return point;
}


