#include "../utils/utils.C" //some useful functions

typedef TGraphAsymmErrors Graph;
enum ratios {RAA, RCP, RR, plpl};
TString G_system_info[]={"Central (0-10%)", "Peripheral (60-80%)", "pp"};

/*
 * TFile* open_sys_file(C_systems system=cent, float R=0.3, float pTlead=5.0,short bining=1, TString evo="eta", bool doToy=0, tests test=camb);
 * TString output_path(C_systems system=cent, short bining=1, TString evo="eta", bool doToy=0);
 * void SetGraphTitles(TGraphAsymmErrors* graph, TString title, TString xtitle, TString ytitle);
 * void SetGraphRangeY(TGraphAsymmErrors* graph, float min, float max);
 * void ShrinkGraph(TGraphAsymmErrors* graph, float min, float max);
 */
void plot_FINAL(float pTlead=5.0, float pTlead_denom=5.0, C_systems system=cent, TString evo="GPC2", TString evo_sys_corr="GPC2", short binning=1, bool showOnePanel=0, bool doToy=0, bool ppBiased=0, TString label="",TString ext="pdf", Int_t plotStyle=0)
{
	
	
	//basic settings
    bool scale_by_pT=1; //scale jet spectrum by 1/pT
	bool drawpp=0; //draw pp baseline in spectra plot
	bool showCorrErr=1; //show correlated errors as well
	bool showTMClosure=1; //show TOYMODEL closure test
	bool showppErr=1; //show pp baseline uncertainty
	bool showSpectrum=1; //plot spectra
	bool spectraInOne=1; //plot spectra in one plot
	//bool showOnePanel=1; //plot the ratio plots into one panel
	
	//label style
	float label_size=0.03;
	Color_t label_color(kGray+1);
	Color_t biasColor=kRed-7; //color of "UNBIASED" marker

	TString trigger="MB";

	float shift=0.9; //0.9; //[0 | ~1 | ~2 ...]; how many first bins do we want to skip in the plots; if shift~1 -> do not plot the first bin in the plots (due to a bad closure test)
	if(doToy)
	{
		showCorrErr=0;
		showppErr=0;
		shift=0; //however we want to see the bad bin on the closure test
	}
	
	const int nr=3; //number of R parameters for spectra
	float Rshow[nr]={0.2,0.3,0.4}; //which R values plot for the spectra
	bool Rshow_onePanel[nr]={1,0,1}; //which R do we want to show if showOnePanel=true -> we can plot several panels into one
	int firstPanel=1; //what is the number of the first panel that will be plotted
	int lastPanel=3; //what is the number of the last panel that will be plotted
	bool showMoreThanOneR;
		int sum=0;
		for(int i=0; i<nr; i++)
		{
			sum+=Rshow_onePanel[i];
			if(showOnePanel) //calculate the position of the first and last panel to plot (we can plot several panels into one)
			{
				if(Rshow_onePanel[i]==0 && firstPanel==(i+1)) firstPanel++;
				if(Rshow_onePanel[i]==1) lastPanel=i+1;
			}
		}
		if(sum>0)showMoreThanOneR=true;
		else showMoreThanOneR=false;
		
		
	
	const int nratios=4; //RAA, RCP, R/R, pTl/pTl
	TString ratio_name[]={"RAA","RCP","RR",Form("pTlpTl%.0lf",pTlead_denom)};
	//if(doToy)ratio_name[3]="pTlpTl";
	if(doToy && showTMClosure)ratio_name[1]="Closure"; //show closure test instead of RCP for toymodel
	int nplots[nratios]={3,3,3,3};
	
	//normalization uncertainty for RAA and RCP - from TAA uncertinty
	float normErr_TAA[2/*RAA,RCP*/][2/*central,peripheral*/]={/*RAA:*/{C_RAA_NormErr[0],C_RAA_NormErr[1]},/*RCP:*/{C_RCP_NormErr,0.0}};
	//normalization uncertainty from Pythia
	float normErr_pyt[nr]={0.22,0.20,0.18};
	bool logRatio[nratios]={1,1,0,1}; //logscale for ratios
	if(doToy && showTMClosure)logRatio[1]=0;
	if(system==pp)logRatio[0]=0;
	//which ratios to show
	bool showRatio[nratios]={1,1,1,1};
	if((/*pTlead>6.0 ||*/ system==peri || system==pp) && !doToy) showRatio[RCP]=0; //show RCP only for central collisions
	if(ppBiased) //show only RAA because other ratios are same as with ppBiased=0
	{
		showSpectrum=0;
		showRatio[RCP]=0;
		showRatio[RR]=0;
		showRatio[plpl]=0;
	}
	
	//Graphics settings
	//using variables from rootlogon.C
	//--------------------------------------------------------
	bool showGrid=0;
	int graphics=1; //1 - better for pdf 0 - better for gif
	gStyle->SetGridColor(gridColorG[graphics]);
	gStyle->SetHatchesLineWidth(hatchesLineWidthG[graphics]);
	gStyle->SetGridWidth(gridLineWidthG[graphics]);
	gStyle->SetGridColor(gridColorG[graphics]);
	gStyle->SetGridStyle(gridStyleG[graphics]);
	
	const Int_t lineLineWidth=lineLineWidthG[graphics];
	const Int_t graphLineWidth=graphLineWidthG[graphics];
	const Int_t fillStyle=fillStyleG[graphics];
	
	
	//colors
	Color_t cspectrum[nr]={kBlue,kMagenta,kRed};
	Color_t cspectrum_corr[nr]={kBlue,kMagenta,kRed};
	Color_t cspectrum_unf[nr]={kGray+1,kGray+1,kGray+1};
	
	Color_t colratio[nr]={kBlue,kBlue,kBlue};
	Color_t colratio_corr[nr]={kMagenta,kMagenta,kMagenta};
	Color_t colratio_unf[nr]={kGray,kGray,kGray};
		if(showOnePanel) //we want to have different colors for different Rs in case of one panel plot
		{
			for (int i=0; i<nr;i++)
			{
				colratio[i]=cspectrum[i];
				colratio_corr[i]=cspectrum[i];
				colratio_unf[i]=cspectrum[i];
			}
		}
	Color_t biasAreaColor=kGray; //color of RAA points in biased region
	
	//markers
	int msz_spec=2;
	int mt_spec[nr]={29,22,33};
	
	int msz_rat=2;
	int mt_rat[nr]={29,29,29};
	
	//lines 
	//width
	int lw_spec=2;
	int lw_spec_corr=1;
	int lw_spec_unf=1;
	//style
	int ls_spec=4;
	int ls_spec_corr=1;
	int ls_spec_unf=1;
	
	//width
	int lw_rat=1;
	int lw_rat_corr=1;
	int lw_rat_unf=1;
	//style
	int ls_rat=1;
	int ls_rat_corr=1;
	int ls_rat_unf=1;
	
	//fills
	int fill_spec_corr=0;
	int fill_spec_unf=1001;
	
	int fill_rat_corr=0;
	int fill_rat_unf=1001;
		if(showOnePanel)fill_rat_unf=0;
	
	//draw styles 
	//spectra
	TString drawspec="P"; //draw style
	TString drawspec_shape="5"; //draw style
	TString drawspec_corr="2"; //draw style
	//ratios
	TString drawrat="P"; //draw style
	TString drawrat_shape="5"; //draw style
	TString drawrat_corr="2"; //draw style
	TString leg_corr="f"; //draw style in legend
	TString leg_shape="f"; //draw style in legend
		if(showOnePanel) 
		{
			drawrat_shape="2";
			drawrat_corr="[]";
			leg_corr="l";
		}

	
	//canvas size
	Int_t can_x=600; //1600
	Int_t can_y=600; //900
	
	//x axis range
	float x_spectra_max=40;
	float x_spectra_min=0;
	if(system==peri || system==pp)
		x_spectra_max=30;
	
	//y axis ranges
	float y_spectra_max=1E-2;
	float y_spectra_min=1E-9;
	if(system==peri || system==pp)
	{
		float y_spectra_max=1E-3;
		float y_spectra_min=1E-9;
	}
	if(scale_by_pT)
    {
        y_spectra_max=y_spectra_max*0.1;
        y_spectra_min=y_spectra_min*0.1;
    }
	
	double spec_scaler[nr]={0.1,1,10}; //for scaling spectra if plotted in one plot
	if(spectraInOne)
	{
		y_spectra_max=y_spectra_max*10;
		y_spectra_min=y_spectra_min*0.1;
	}
	
	float y_ratio_max[C_nsystems][nratios]={
		//"RAA","RCP","RR","pTlpTl"
		{2.0,2.0,1.6,2.0},	//central
		{2.0,2.0,1.6,2.0},	//peripheral
		{3.0,10.0,1.4,2.0}};	//pp
	float y_ratio_min[C_nsystems][nratios]={
		{0.05,0.05,0.4,0.05},	//central
		{0.05,0.05,0.4,0.05},	//peripheral
		{0.0,0.1,0.4,0.0}};	//pp
		if(doToy && showTMClosure) 
		{
			y_ratio_max[0][1]=2.8;
			y_ratio_min[0][1]=0.0;
		}
			
		//y position of "x10" and "x0.1" labels on spectra plots
		float pTscl=(scale_by_pT)? 0.95 : 1.0;
		float yTimes01[C_nsystems]={0.28*pTscl,0.23*pTscl,0.2*pTscl};
		float yTimes10[C_nsystems]={0.48*pTscl,0.44*pTscl,0.4*pTscl};
		float xTimes[C_nsystems]={0.75,0.8,0.8};
		
		float x_ratio_min[/*panels*/]={0,0.05,0.05}; //range x axis of ratio plots; we don't want to plot the "0" axis lable on panels 2 and 3
		if(doToy && showTMClosure && showOnePanel) 
			x_ratio_min[0]=pTlead-1;
		float x_ratio_max[C_nsystems][nratios]={
			{37,37,37,37}, //central
			{37,37,37,37}, //peripheral
			{32,32,32,32}}; //pp
		float xlineMin=x_ratio_min[0]; //range for drawing the line in the ratio plots
		float xlineMax; //range for drawing the line in the ratio plots
		float xbiasLine[nr]={14,16,16}; //x position of a line showing pTlead bias range in RAA
		//y range of a line showing pTlead bias range in RAA and plpl
		float ybiasLine_max[C_nsystems][nratios]={
			{0.3,0,0,1.0}, //central
			{1.0,0,0,1.0}, //peripheral
			{0,0,0,0}}; //pp
		float ybiasLine_min[C_nsystems][nratios]={
			{0.05,0.0,0.0,0.3},	//central
			{0.3,0.0,0.0,0.3},	//peripheral
			{0.0,0.0,0.0,0.0}};	//pp
		float xbiasDesc[nr]={0.55,0.52,0.52}; //x position of the bias textbox
		//x position of the bias textbox
		float ybiasDesc[C_nsystems][nratios]={
			{0.3,0.0,0.0,0.6}, 	//central
			{0.58,0.0,0.0,0.6},	//peripheral
			{0.0,0.0,0.0,0.0}};	//pp
		float boxX[nratios]={30,2,30,30}; //where to draw normalization error box
		float yline[nratios]={0.2,0.35,0.4,0}; //y value of the line in the ratio plots
		bool drawLine=0; //draw line (in addition to y=1), which marks e.g. some previous results
		float xshowMax[C_nsystems][nratios]={ //maxmimal value up to which draw the points in ratios 
			{30.0, 25.0, 30.0, 30.0}, //central
			{25.0, 20.0, 25.0, 25.0}, //peripheral
			{30.0, 20.0, 30.0, 30.0}}; //pp
		if(doToy && showTMClosure) xshowMax[0][1]=30;
		//description of graphs and axis
		TString xtitle="p_{T,jet}^{charged} (GeV/#it{c})";
		TString ytitle="1/N_{events} 1/2#pi d^{2}N/dp_{T}^{ch}d#eta (GeV/#it{c})^{-1}";
        if(scale_by_pT) ytitle="1/(p_{T}N_{events}) 1/2#pi d^{2}N/dp_{T}^{ch}d#eta (GeV/#it{c})^{-1}";
		TString ytitle_rat[nratios]={"R_{AA}^{Pythia}","R_{CP} (0-10%)/(60-80%)","R1/R2","pTlead1/pTlead2"};
		TString leg_desc_rat[2][nratios]={"STAR Au+Au/PYTHIA6","STAR (0-10%)/(60-80%)","STAR data",
			"STAR data","Param. model/PYTHIA","unfolded/generated","Parametrized model","Parametrized model"};
		if(system==pp) ytitle_rat[0]="STAR pp/PYTHIA";
		if(doToy && showTMClosure)ytitle_rat[1]="unfolded/generated";
		//float pTlead_base[C_nsystems]={5.0,4.0,4.0}; //pTlead used as denominator in pTlead ratios. This value is used only for y axis title.
		TString outdir=output_path(system, binning, evo, doToy,ppBiased);
		
		int pyt_rebin=1; //rebin pythia spectra used for pp pTlpTl ratio
					
					
		//TGraphs
		Graph* gspectrum_avg[nr];
		Graph* gspectrum_sys_corr[nr];
		Graph* gspectrum_sys_unf[nr];
		
		Graph* gratio_avg[nr][nratios];
		Graph* gratio_sys_corr[nr][nratios];
		Graph* gratio_sys_unf[nr][nratios];
		Graph* gratio_sys_norm[nr];
		
		TFile* fin[nr];
		TFile* fin_sys_corr[nr];
		TString str;
		
		//------------------------------------------
		//load files
		//------------------------------------------
		for(int r=0; r<nr; r++)
		{
			float R=C_Rs[r];
			TString suff=(doToy) ? "" : "_BSL";
			fin[r]=open_sys_file(system, R, pTlead,binning,evo,doToy, ppBiased); //file with input histograms/graphs
			fin_sys_corr[r]=open_sys_file(system, R, pTlead,binning,evo_sys_corr,doToy, ppBiased); //file with correlated errors
			gspectrum_avg[r]=(Graph*)fin[r]->Get(Form("spectrum%s",suff.Data()));
			if(showCorrErr) gspectrum_sys_corr[r]=(Graph*)fin_sys_corr[r]->Get("spectrum_corr");
			gspectrum_sys_unf[r]=(Graph*)fin[r]->Get("spectrum_unfold");
			
			//if(showppErr)gratio_sys_norm[r]=(Graph*)fin[r]->Get(Form("RAA_norm%s",suff.Data()));
			for(int ratio=0; ratio<nratios; ratio++)
			{
				if(!showRatio[ratio]) continue;
				gratio_avg[r][ratio]=(Graph*)fin[r]->Get(Form("%s%s",ratio_name[ratio].Data(),suff.Data()));
				cout<<"R="<<R<<", ratio:"<<ratio_name[ratio].Data()<<", N="<<gratio_avg[r][ratio]->GetN()<<endl;
				if(showCorrErr) gratio_sys_corr[r][ratio]=(Graph*)fin_sys_corr[r]->Get(Form("%s_corr%s",ratio_name[ratio].Data(),suff.Data()));
				gratio_sys_unf[r][ratio]=(Graph*)fin[r]->Get(Form("%s_unfold%s",ratio_name[ratio].Data(),suff.Data()));
				cout<<"Nsys="<<gratio_sys_unf[r][ratio]->GetN()<<endl;
			}
		}
		
		//=================================
		//DRAW
		//=================================
		
		TLatex *latex = new TLatex();
		latex->SetNDC();
		latex->SetTextSize(0.045);
		
		TLatex *latexL = new TLatex();
		latexL->SetNDC();
		latexL->SetTextSize(label_size);
		latexL->SetTextColor(label_color);
		
		TString sNplots;
		
		//------------------
		//spectra
		//------------------
		gStyle->SetTitleSize(0.055,"Y");
		gStyle->SetTitleOffset(1.30,"Y");
		gStyle->SetTitleSize(0.06,"X");
		gStyle->SetTitleOffset(1.0,"X");
		gStyle->SetLabelSize(0.042,"X");
		gStyle->SetLabelSize(0.042,"Y");
        if(scale_by_pT)
        {
            gStyle->SetTitleOffset(1.4,"Y");
            gStyle->SetTitleSize(0.05,"Y");
        }
		if(showSpectrum){
			TCanvas *cspec;
			if(spectraInOne) //all spectra (for each R) in one plot
				cspec= new TCanvas("cspec","cspec",10,10,can_x,can_y);
			else
			{
				if(nr==4)
				{
					cspec= new TCanvas("cspec","cspec",10,10,2*can_x,2*can_y);
					cspec->Divide(2,2);
				}
				else 
				{
					cspec= new TCanvas("cspec","cspec",10,10,nr*can_x,can_y);
					cspec->Divide(nr,1);
				}
			}//spectra in separate plots for each R
			for(int r=0; r<nr; r++)
			{
				if(spectraInOne) cspec->cd();
				else cspec->cd(r+1);
				if(showGrid) gPad->SetGrid();
				gPad->SetLogy();
				
				//rescale the spectra so we can distinguish them
				if(spectraInOne)
				{
					scale_graph(gspectrum_avg[r],spec_scaler[r]);
					scale_graph(gspectrum_sys_unf[r],spec_scaler[r]);
					if(showCorrErr)scale_graph(gspectrum_sys_corr[r],spec_scaler[r]);
				}
				
				//remove high pT points which we don't want to show
				ShrinkGraph(gspectrum_avg[r],pTlead+shift,xshowMax[system][0]); //remove the first datapoint (due to a bad closure test) and a few last points with bad statistics
				ShrinkGraph(gspectrum_sys_unf[r],pTlead+shift,xshowMax[system][0]);
				if(showCorrErr) ShrinkGraph(gspectrum_sys_corr[r],pTlead+shift,xshowMax[system][0]);

                //scale the whole spectrum by 1/pT in order to make it invariant
                if(scale_by_pT)
                {
                    double* xarr=gspectrum_avg[r]->GetX();
                    scale_graph(gspectrum_avg[r],xarr,1);
                    if(showCorrErr) scale_graph(gspectrum_sys_corr[r],xarr,1);
                    scale_graph(gspectrum_sys_unf[r],xarr,1);
                }
                
				//--------------------
				//set graph properties
				//--------------------
				int idx=(spectraInOne)? r : 0;
				tgraph_set_atributes(gspectrum_avg[r], x_spectra_min, x_spectra_max, y_spectra_min, y_spectra_max, cspectrum[idx], cspectrum[idx], cspectrum[idx], lw_spec, msz_spec, mt_spec[idx], 0, ls_spec);
				if(showCorrErr) tgraph_set_atributes(gspectrum_sys_corr[r], x_spectra_min, x_spectra_max, y_spectra_min, y_spectra_max, cspectrum_corr[idx], cspectrum_corr[r], cspectrum_corr[r], lw_spec_corr, msz_spec, mt_spec[idx], fill_spec_corr, ls_spec_corr);
				tgraph_set_atributes(gspectrum_sys_unf[r], x_spectra_min, x_spectra_max, y_spectra_min, y_spectra_max, cspectrum_unf[idx], cspectrum_unf[r],  cspectrum_unf[r], lw_spec_unf, msz_spec, mt_spec[idx], fill_spec_unf, ls_spec_unf);
				
				
				SetGraphTitles(gspectrum_sys_unf[r],"",xtitle,ytitle);
				//SetGraphRangeY(gspectrum_sys_unf[r],y_spectra_min, y_spectra_max);
				
				TString draw_opt="a"+drawspec_shape;
				if(spectraInOne && r>0) draw_opt=drawspec_shape;
				gspectrum_sys_unf[r]->Draw(draw_opt);
				if(showCorrErr) gspectrum_sys_corr[r]->Draw(drawspec_corr);
				gspectrum_avg[r]->Draw(drawspec);
				
				if(!spectraInOne) latex->DrawLatex(0.7, 0.75,Form("R=%.1lf",Rshow[r]));
				else
				{
					 latex->DrawLatex(xTimes[system], yTimes10[system],"#times 10");
					 latex->DrawLatex(xTimes[system], yTimes01[system],"#times 0.1");
				}
				latexL->DrawLatex(0.18, 0.85,label);
				
				//legends
				if((r==1 && !spectraInOne) || (r==nr-1 && spectraInOne))
				{
					double posx0=0.2,posx1=0.16, posx2=0.85, posy0=0.2,posy1=0.35,posy2=0.5;
					if(spectraInOne) {posx0=0.6, posx1=0.55;posx2=0.95;posy0=0.6;posy1=0.75;posy2=0.9;}
					if(doToy) posy0+=0.07;
					TLegend *legspec = new TLegend(posx1,posy1,posx2,posy2);
					legspec->SetTextSize(0.04);
					legspec->SetFillStyle(0);
					legspec->SetBorderSize(0);
					
					
					if(spectraInOne)
					{
						for(int rr=0; rr<nr;rr++)
							legspec->AddEntry(gspectrum_avg[rr], Form("R=%.1lf",Rshow[rr]/*,spec_scaler[rr]*/),"lp");
					}
					else
						legspec->AddEntry(gspectrum_avg[1], "STAR Au+Au","lp");
					legspec->DrawClone("same");
					
					TLegend *legspec2 = new TLegend(posx0,posy0,posx2,posy1);
					legspec2->SetTextSize(0.035);
					legspec2->SetFillStyle(0);
					legspec2->SetBorderSize(0);
					if(showCorrErr) 
					{
						legspec2->AddEntry(gspectrum_sys_corr[0], "correlated unc.", "f");
						legspec2->AddEntry(gspectrum_sys_unf[0], "shape unc.", "f");
					}
					else legspec2->AddEntry(gspectrum_sys_unf[0], "unfolding unc.", "f"); 
					legspec2->DrawClone("same");
					
				}
				
				if(r==0 || spectraInOne){
					TLegend *model_info = new TLegend(0.18, 0.2, 0.4, 0.45);
					model_info->SetTextSize(0.035);
					model_info->SetFillStyle(0);
					model_info->SetBorderSize(0);
					model_info->SetMargin(0.05);
					if(!doToy) //data
					{
						if(system!=pp)
						{
							model_info->SetHeader("STAR Au+Au");
							//model_info->AddEntry("", Form("Run11, %s",trigger.Data()), "");
						}
						else 
						{
							model_info->SetHeader("STAR p+p");
							//model_info->AddEntry("", "Run12, MB+HT", "");
						}
						model_info->AddEntry("", "#sqrt{s_{NN}}=200 GeV", "");
					}
					else //toymodel
					{
						model_info->SetHeader("PARAMETRIZED MODEL");
					}
					model_info->AddEntry("", "Charged jets, anti-k_{T}", "");
					model_info->AddEntry("", G_system_info[system], "");
					model_info->AddEntry("", Form("p_{T,lead}^{min} = %.0lf GeV/#it{c}",pTlead), "");
					model_info->DrawClone("same");
				}
				
			}//R loop
			
			sNplots="all_R";
			if(spectraInOne) sNplots="all_R_in1";
			str = Form("%s/%s_spec_pTlead%.0lf_bin%i.%s", outdir.Data(),sNplots.Data(),pTlead,binning,ext.Data());
			cspec->SaveAs(str.Data());
		}//show spectrum
		//------------------
		//Ratio
		//------------------
		//the following gStyle settings are probably not working at this place 
		gStyle->SetTitleSize(0.07,"Y"); 
		gStyle->SetTitleOffset(0.8,"Y");
		gStyle->SetTitleSize(0.07,"X");
		gStyle->SetTitleOffset(1.05,"X");
		gStyle->SetLabelSize(0.05,"X");
		gStyle->SetLabelSize(0.05,"Y");
		
		TFile* fpyt=new TFile("../../plotting_out/root/pp_PYTHIA_Jana/jethistos.root","OPEN"); //file with baseline PYTHIA for pythia pTlead ratios
		
		//fit function for pp pTlpTl ratio
		TF1 *fpyfit = new TF1("fpyfit","[0]*(1.0-TMath::Exp([1]*x+[2]))",6,35);
      fpyfit->SetRange(pTlead,xshowMax[system][plpl]);
      fpyfit->SetParameters(1.0,-0.3,1.0);
		//histograms for pp pTlpTl ratio
		TH1D* hnum[nr];
		TH1D* hdenom[nr];
		TH1D* hdivide[nr];
		//pp pTlead ratio for comparison
		for(int rr=0; rr<nr; rr++)
		{
			hdivide[rr]=(TH1D*) fpyt->Get(Form("hjet_pT_R0%0.lf_pTl%.0lf_sum_fudge",Rshow[rr]*10,pTlead));
			hdivide[rr]->Rebin(pyt_rebin);
			hdenom[rr]=(TH1D*) fpyt->Get(Form("hjet_pT_R0%0.lf_pTl%.0lf_sum_fudge",Rshow[rr]*10,pTlead_denom));
			hdenom[rr]->Rebin(pyt_rebin);
			hdivide[rr]->Divide(hdenom[rr]);
			hdivide[rr]->GetXaxis()->SetRangeUser(pTlead, xshowMax[system][plpl]);
		}
		TCanvas *cratio[nratios];
		for(int ratio=0; ratio<nratios; ratio++)
		{
			if(!showRatio[ratio]) continue;
			str=Form("cratio_%i",ratio);
			
			xlineMax=x_ratio_max[system][ratio];
			//if(system==peri) boxX[ratio]=boxX[ratio]-5;
			
			if(showOnePanel) //single-panel plot
			{
                cratio[ratio]= new TCanvas(str,ratio_name[ratio],10,10,1.0*can_x,1.0*can_y);
				cratio[ratio]->cd();
				
			}
			else if(nplots[ratio]==4) 
			{
				cratio[ratio]= new TCanvas(str,ratio_name[ratio],10,10,2*can_x,2*can_y);
				cratio[ratio]->Divide(2,2,0,0);
			}
			else
			{	   
				cratio[ratio]= new TCanvas(str,ratio_name[ratio],10,10,nplots[ratio]*can_x,can_y);
				cratio[ratio]->Divide(nplots[ratio],1,0,0);
			}
			TLegend *legratio_R;
			
			for(int panel=1; panel<nr+1; panel++)//loop over panels
			{
				int ridx=panel-1;
				
				if(ratio==RR) 
				{
					ridx=panel;
					if(panel==lastPanel)ridx=0;
				}
				if(panel>nplots[ratio]) continue; //last panel has been already reached
				if(showOnePanel && (Rshow_onePanel[ridx]==0)) continue; //we don't want to display this R
				float rshow=Rshow[ridx];
				if((!showOnePanel) || panel==1) //skip 2nd and 3rd panel in showOnePanel plot
				{
					cratio[ratio]->cd(panel);
					if(showGrid) gPad->SetGrid();
					if(logRatio[ratio]) gPad->SetLogy();
					//Set margins (zero margins between panels)
					gPad->SetLeftMargin(0.13);
					if(ratio==plpl) gPad->SetLeftMargin(0.16);
					gPad->SetRightMargin(0.02);
					if(ratio!=RR && !showOnePanel)
					{
						if(panel==firstPanel) 
							gPad->SetRightMargin(0.00);
						else if ((panel==nplots[ratio])) 
							gPad->SetLeftMargin(0.00);
						else
						{
							gPad->SetLeftMargin(0.00);
							gPad->SetRightMargin(0.00);
						}
						
					}        
					else if(showOnePanel)
                    {
                        gratio_sys_unf[ridx][ratio]->GetXaxis()->SetTitleSize(0.065);
                        gratio_sys_unf[ridx][ratio]->GetXaxis()->SetLabelSize(0.045);
						gratio_sys_unf[ridx][ratio]->GetXaxis()->SetTitleOffset(1.0);
                        gratio_sys_unf[ridx][ratio]->GetYaxis()->SetTitleSize(0.065);
                        gratio_sys_unf[ridx][ratio]->GetYaxis()->SetLabelSize(0.045);
						gratio_sys_unf[ridx][ratio]->GetYaxis()->SetTitleOffset(0.85);
                    }
					if(ratio==RR) 
					{
						ytitle_rat[2]=Form("R=0.2/R=%.1lf",rshow);
						if(panel==nplots[RR]/*last panel*/) ytitle_rat[2]="R=0.3/R=0.4";
						if(evo=="eta_RR0403") ytitle_rat[2]=Form("R=0.3/R=%.1lf",rshow);
					}
					else if(ratio==plpl) 
					{
						ytitle_rat[ratio]=Form("p_{T,lead}^{min}=%.0lf GeV/#it{c} / p_{T,lead}^{min}=%.0lf GeV/#it{c}",pTlead,pTlead_denom);
						gratio_sys_unf[ridx][ratio]->GetYaxis()->SetTitleSize(0.055);
						gratio_sys_unf[ridx][ratio]->GetYaxis()->SetTitleOffset(1.2);
					}
				}//skip 2nd and 3rd panel in showOnePanel plot
				ShrinkGraph(gratio_avg[ridx][ratio],pTlead+shift,xshowMax[system][ratio]); //remove the first datapoint (due to a bad closure test) and a few last points with bad statistics
				ShrinkGraph(gratio_sys_unf[ridx][ratio],pTlead+shift,xshowMax[system][ratio]);
				if(showCorrErr) ShrinkGraph(gratio_sys_corr[ridx][ratio],pTlead+shift,xshowMax[system][ratio]);
				//--------------------
				//set graph properties
				//--------------------
				int idx=ridx;
				//show xtitle only on last panel
				TString xtl=xtitle;
				if(panel<nr && ratio!=RR && !showOnePanel) xtl="";
				tgraph_set_atributes(gratio_avg[ridx][ratio], x_ratio_min[panel-1], x_ratio_max[system][ratio], y_ratio_min[system][ratio], y_ratio_max[system][ratio], colratio[idx], colratio[idx], colratio[idx], lw_rat, msz_rat, mt_rat[idx], 0, ls_rat);
				if(showCorrErr) tgraph_set_atributes(gratio_sys_corr[ridx][ratio], x_ratio_min[panel-1], x_ratio_max[system][ratio], y_ratio_min[system][ratio], y_ratio_max[system][ratio], colratio_corr[idx], colratio_corr[idx], colratio_corr[idx], lw_rat_corr, msz_rat, mt_rat[idx], fill_rat_corr, ls_rat_corr);
				tgraph_set_atributes(gratio_sys_unf[ridx][ratio], x_ratio_min[panel-1], x_ratio_max[system][ratio], y_ratio_min[system][ratio], y_ratio_max[system][ratio], colratio_unf[idx], colratio_unf[idx], colratio_unf[idx], lw_rat_unf, msz_rat, mt_rat[idx], fill_rat_unf, ls_rat_unf);
				SetGraphTitles(gratio_sys_unf[ridx][ratio],"",xtl,ytitle_rat[ratio]);
				SetGraphRangeY(gratio_sys_unf[ridx][ratio],y_ratio_min[system][ratio],y_ratio_max[system][ratio]);
				gratio_sys_unf[ridx][ratio]->GetXaxis()->SetLimits(x_ratio_min[panel-1],x_ratio_max[system][ratio]);
				
				//********************
				//draw
				//********************
				TString draw_axis="a";
				if(showOnePanel && panel>1) draw_axis="";
				gratio_sys_unf[ridx][ratio]->DrawClone(draw_axis+drawrat_shape);
				if(showCorrErr) gratio_sys_corr[ridx][ratio]->DrawClone(drawrat_corr);
				gratio_avg[ridx][ratio]->DrawClone(drawrat);
				
				//replot biased RAA region only in gray color
				if(ratio==RAA)
				{
					ShrinkGraph(gratio_sys_unf[ridx][ratio],0,xbiasLine[ridx]);
					gratio_sys_unf[ridx][ratio]->SetFillColor(biasAreaColor);
					gratio_sys_unf[ridx][ratio]->DrawClone(drawrat_shape);
					gratio_sys_unf[ridx][ratio]->SetFillColor(colratio_unf[idx]); //reset color for legend
					
					if(showCorrErr){
					ShrinkGraph(gratio_sys_corr[ridx][ratio],0,xbiasLine[ridx]);
					gratio_sys_corr[ridx][ratio]->SetLineColor(biasAreaColor+1);
					gratio_sys_corr[ridx][ratio]->DrawClone(drawrat_corr);
					gratio_sys_corr[ridx][ratio]->SetLineColor(colratio_corr[idx]); //reset color for legend
					}
					ShrinkGraph(gratio_avg[ridx][ratio],0,xbiasLine[ridx]);
					gratio_avg[ridx][ratio]->SetMarkerColor(biasAreaColor+1);
					gratio_avg[ridx][ratio]->SetLineColor(biasAreaColor+1);
					gratio_avg[ridx][ratio]->DrawClone(drawrat);
				}
				TLine *one = new TLine(xlineMin, 1, xlineMax,1);
				one->SetLineWidth(lineLineWidth);
				one->SetLineStyle(2);
				one->SetLineColor(kGray+1);
				one->DrawClone("same");
				
				TLine *line = new TLine(xlineMin, yline[ratio], xlineMax, yline[ratio]);
				line->SetLineWidth(lineLineWidth);
				line->SetLineStyle(2);
				line->SetLineColor(kGray+2);
				if(drawLine) line->DrawClone("same");
				//xbiasLine=1.75*(0.6+rshow)*(pTlead+4);
                double xb1=xbiasLine[ridx];
                double xb2=xbiasLine[ridx];
                double yb1=ybiasLine_min[system][RAA];
                double yb2=ybiasLine_max[system][RAA];
                if(ratio==plpl)
                {
						 yb1=ybiasLine_min[system][plpl];
						 yb2=ybiasLine_max[system][plpl];
                    //yb1=y_ratio_min[system][plpl];
                    //yb2=y_ratio_max[system][plpl];
                }
				TLine *bias = new TLine(xb1,yb1,xb2,yb2);
				bias->SetLineWidth(lineLineWidth);
				bias->SetLineStyle(2);
				bias->SetLineColor(biasColor);
                TLatex *latexbias = new TLatex();
				latexbias->SetNDC();
				latexbias->SetTextSize(0.05);
				latexbias->SetTextColor(biasColor);
				if(ratio==RAA && ppBiased==0) 
				{
                    bias->DrawClone("same");
                    double xbdesc=xbiasDesc[ridx];
                    if(system==peri) xbdesc=xbiasDesc[ridx]+0.1;
					latexbias->DrawLatex(xbiasDesc[ridx], ybiasDesc[system][ratio],"--> ~ UNBIASED");
					//latexbias->DrawLatex((pTlead/10.*0.22)+(0.4-rshow)/3.0, 0.5,"p_{T}^{lead} bias <--");
				}
				if(ratio==plpl) 
				{
					hdivide[ridx]->Fit("fpyfit","0",""); //fit pythia plpl ratio histogram with function
					//hdivide[ridx]->Draw("histo c same");
					fpyfit->SetLineColor(kOrange+9);
					fpyfit->DrawClone("same");
					fpyfit->SetParameters(1.0,-0.3,1.0); //reset parameters before next fit
					
					bias->DrawClone("same");
					latexbias->DrawLatex(xbiasDesc[ridx], ybiasDesc[system][ratio],"--> ~ UNBIASED");
				}
				TLatex *latexr = new TLatex();
				latexr->SetNDC();
				double lat_tex_size=0.06;
				double tex_x=0.8, tex_y=0.9;
				//if(ratio==RAA && system==peri)tex_x=0.2;
				if(ratio==plpl)tex_x=0.2;
				//if(ratio==RCP) {tex_x=0.8; tex_y=0.9;}
				if(showOnePanel){ tex_x=0.15; tex_y=0.8;lat_tex_size=0.04;}
				if(panel==firstPanel)tex_x+=0.05;
				latexr->SetTextSize(lat_tex_size);
				if(!showOnePanel && ratio!=RR) latexr->DrawLatex(tex_x, tex_y,Form("R=%.1lf",rshow));
				TBox* ppbox_TAA; //sys uncertainty from TAA
				TBox* ppbox_pyt; //sys uncertainty of PYTHIA baseline
				//Draw pp baseline error bar for RAA and RCP
				if(showppErr && (ratio==RCP || ratio==RAA))
				{
					ppbox_TAA=new TBox(boxX[ratio], 1-normErr_TAA[ratio][system], boxX[ratio]+1, 1+normErr_TAA[ratio][system]);
					ppbox_TAA->SetFillColor(kOrange+6);
					ppbox_TAA->SetFillStyle(fillStyle);
					ppbox_TAA->DrawClone("");
				}
				if(showppErr && ratio==RAA)
				{
					ppbox_pyt=new TBox(boxX[ratio]-2, 1-normErr_pyt[ridx], boxX[ratio]-1, 1+normErr_pyt[ridx]);
					ppbox_pyt->SetFillColor(kCyan+1);
					ppbox_pyt->SetFillStyle(fillStyle);
					ppbox_pyt->DrawClone("");
				}
				//legends
				double posx1, posx2, posy1, posy2;
				float leg_tex_size=0.055;
				if(showOnePanel) leg_tex_size=0.05;
																									 
				if(ratio==plpl && !logRatio[plpl]) {posx1=0.18; posy1=0.75; posx2=0.5; posy2=0.97;}
				else if(ratio==plpl && logRatio[plpl]) {posx1=0.05; posy1=0.2; posx2=0.45; posy2=0.45;}
				else if(ratio==RR) {posx1=0.40; posy1=0.75; posx2=0.8; posy2=0.97;}
				else if(ratio==RCP) {posx1=0.05; posy1=0.2; posx2=0.45; posy2=0.5;}
				else if(ratio==RAA && system==cent) {posx1=0.18; posy1=0.68; posx2=0.5; posy2=0.98;}
				else if(ratio==RAA) {posx1=0.1; posy1=0.2; posx2=0.45; posy2=0.5;}
				else {posx1=0.18; posy1=0.7; posx2=0.4; posy2=0.97;}
				if(ratio==RCP && showTMClosure && doToy) {posx1=0.5; posy1=0.75; posx2=0.8; posy2=0.8;}
				
				if((panel==2 || (showOnePanel && panel==lastPanel)) && evo!="eta_RR0403"){

					//legend - uncertainties
					TLegend *legratio_err = new TLegend(posx1,posy1,posx2,posy2);
					legratio_err->SetTextSize(leg_tex_size);
					legratio_err->SetFillStyle(0);
					legratio_err->SetBorderSize(0);
					if(showCorrErr)
					{
						legratio_err->AddEntry(gratio_sys_corr[ridx][ratio], "correlated unc.", leg_corr);
						legratio_err->AddEntry(gratio_sys_unf[ridx][ratio], "shape unc.", leg_shape);
					}
					//else legratio_err->AddEntry(gratio_sys_unf[ridx][ratio], "unfolding uncert.", leg_shape);
					if(!doToy && showppErr && (ratio==RCP || ratio==RAA))legratio_err->AddEntry(ppbox_TAA, "T_{AA} uncertainty", "f");
					if(!doToy && showppErr && ratio==RAA)legratio_err->AddEntry(ppbox_pyt, "Pythia uncertainty", "f");																		 
					legratio_err->DrawClone("same");
					//if(showppErr)legratio->AddEntry(gratio_sys_norm[0], "pp baseline uncertainty", "f");
				}
				
				if(showOnePanel && showMoreThanOneR)
				{
					if(panel==firstPanel)
					{
						legratio_R= new TLegend(posx1+0.2,posy1,posx2+0.2,posy2+0.1);
						legratio_R->SetTextSize(leg_tex_size);
						legratio_R->SetFillStyle(0);
						legratio_R->SetBorderSize(0);
					}
					legratio_R->AddEntry(gratio_avg[ridx][ratio], Form("R=%.1lf",rshow));
					if(panel==lastPanel)
						legratio_R->DrawClone("same");
				}
					
				
				if(panel==lastPanel && ratio==plpl)
				{
					posx1=0.18; posy1=0.75; posx2=0.5; posy2=0.97;
					if(logRatio[plpl]) {posx1=0.5; posy1=0.2; posx2=0.9; posy2=0.45;}
					TLegend *legplpl = new TLegend(posx1,posy1,posx2,posy2);
					legplpl->SetTextSize(leg_tex_size);
					legplpl->SetFillStyle(0);
					legplpl->SetBorderSize(0);
					legplpl->AddEntry(gratio_avg[0][plpl],"STAR Au+Au","p");
					legplpl->AddEntry(fpyfit,"PYTHIA p+p","l");
					legplpl->DrawClone("same");
				}
					
				//legend - general information
				if((panel==firstPanel) && evo!="eta_RR0403"){
					if(ratio==plpl && !logRatio[plpl]) {posx1=0.18; posy1=0.7; posx2=0.5; posy2=0.97;}
					else if(ratio==plpl && logRatio[plpl]) {posx1=0.2; posy1=0.2; posx2=0.5; posy2=0.52;}
					else if(ratio==RR) {posx1=0.18; posy1=0.7; posx2=0.4; posy2=0.97;}
					else if(ratio==RCP){posx1=0.18; posy1=0.20; posx2=0.5; posy2=0.52;}
					else if(ratio==RAA && system==cent) {posx1=0.18; posy1=0.62; posx2=0.4; posy2=0.97;}
					else if(ratio==RAA){posx1=0.18; posy1=0.2; posx2=0.4; posy2=0.55;}
					else {posx1=0.18; posy1=0.7; posx2=0.4; posy2=0.97;}
				   if(ratio==RCP && showTMClosure && doToy) {posx1=0.2; posy1=0.55; posx2=0.4; posy2=0.89;}
					
					TLegend *model_info = new TLegend(posx1,posy1,posx2,posy2);
					if(ratio==RCP && showTMClosure && doToy)leg_tex_size=0.05;
					if(showOnePanel) leg_tex_size=0.035;
					model_info->SetTextSize(leg_tex_size);
					model_info->SetFillStyle(0);
					model_info->SetBorderSize(0);
					model_info->SetMargin(0.05);
					if(!doToy) //data
					{
						if(system!=pp)
						{
							model_info->SetHeader("STAR Au+Au #sqrt{s_{NN}}=200 GeV");
							//model_info->AddEntry("", Form("Run11, %s",trigger.Data()), "");
						}
						else 
						{
							model_info->SetHeader("STAR p+p #sqrt{s_{NN}}=200 GeV");
							//model_info->AddEntry("", "Run12, MB+HT", "");
						}
					}
					else //toymodel
					{
						model_info->SetHeader("Parametrized Model");
						model_info->AddEntry("", "Au+Au #sqrt{s_{NN}}=200 GeV", "");
						model_info->AddEntry("", "Central (0-10%)", "");
					}
					model_info->AddEntry("", "Charged jets, anti-k_{T} ", "");
					if(ratio!=RCP) model_info->AddEntry("", G_system_info[system], "");
					if(showOnePanel && ratio!=RR && !showMoreThanOneR) model_info->AddEntry("",Form("R=%.1lf",rshow), "");
					//if(showOnePanel && ratio!=RR) model_info->AddEntry("",Form("anti-k_{T} , R=%.1lf",rshow), "");
					//else model_info->AddEntry("","anti-k_{T} ", "");
					if(ratio!=plpl) model_info->AddEntry("", Form("p_{T,lead}^{min} = %.0lf GeV/#it{c}",pTlead), "");
					model_info->AddEntry("", "p_{T}^{const} > 0.2 GeV/#it{c}", "");
					model_info->DrawClone("same");
				}
				
				if(panel==lastPanel || (ratio==RR && nplots[RR]==2 && panel==2))
				{
					float labx=0.5;
					float laby=0.3;
					if(ratio==RR)labx=0.35;
					if(ratio==RAA)laby=0.8;																 
					latexL->DrawLatex(labx, laby,label);
				}
			}//panel loop
			sNplots="all_R";
			if(showOnePanel) sNplots="onePanel";
			TString ppBiasSuf="";
			if(ppBiased) ppBiasSuf="_ppBiased";
			str = Form("%s/%s_%s_pTlead%.0lf_bin%i%s.%s", outdir.Data(),sNplots.Data(),ratio_name[ratio].Data(),pTlead,binning,ppBiasSuf.Data(),ext.Data());
			cratio[ratio]->SaveAs(str.Data());
			
			//delete objects
			delete model_info;
			//delete legratio;
			//delete legratio_err;
			if(showSpectrum) delete cspec;
			/*delete cratio_0;
			delete cratio_1;
			delete cratio_2;
			delete cratio_3;*/
		}//ratio loop
		
		//close files
		for(int r=0; r<nr; r++)
		{
			fin[r]->Close();
			delete fin[r];
		}
		
}


//***********************************************************************************************************
//utility functions
//***********************************************************************************************************
TFile* open_sys_file(C_systems system=cent, float R=0.3, float pTlead=5.0,short bining=1, TString evo="eta", bool doToy=0, /*tests test=camb*/ bool pp_pTlead=0)
{
	TString sdir="../../plotting_out/systematics";
	if(doToy) sdir+="/toymodel";
	if(system==peri) sdir+="/peripheral";
	else if(system==pp) sdir+="/pp";
	else sdir+="/central";
	sdir+=Form("/%s",evo.Data());
	if(pp_pTlead) sdir+="_ppBiased";
	sdir+=Form("/bining%i",bining);
	//sdir+=Form("/%s",C_test_name[test].Data());
	TString str=Form("%s/systematics_R%.1lf_pT%0.lf.root",sdir.Data(),R,pTlead);
	cout<<"opening file:"<<str.Data()<<endl;
	TFile *fsys;
	fsys= new TFile(str.Data(), "OPEN");
	return fsys;
}

TString output_path(C_systems system=cent, short bining=1, TString evo="eta", bool doToy=0, bool pp_pTlead=0)
{
	TString sdir=Form("../../plotting_out/obr/%s/results",evo.Data());
	if(doToy) sdir+="/toymodel";
	if(system==peri) sdir+="/peripheral";
	else if(system==pp) sdir+="/pp";
	else sdir+="/central";
	//if(pp_pTlead) sdir+="_ppBiased";
	//sdir+=Form("/bining%i",bining);
	return sdir;
}

void SetGraphTitles(TGraphAsymmErrors* graph, TString title, TString xtitle, TString ytitle)
{
	graph->SetTitle(title);
	graph->GetXaxis()->SetTitle(xtitle);
	graph->GetYaxis()->SetTitle(ytitle);
	return;
}

void SetGraphRangeY(TGraphAsymmErrors* graph, float min, float max)
{
	graph->GetHistogram()->SetMaximum(max);
	graph->SetMinimum(min);
	
	return;
}

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
