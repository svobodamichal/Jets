#include "../../utils/common.h"
enum comparison { all, theory, hadrons }; //what do we want to compare our results to?
TString comp_suff[]={"","_theory","_hadrons"};

void RAA_RCP(C_systems system=cent, double pTlead=5.0, int binning=1,TString evo="GPC2", TString label="", bool ppBias=0, comparison compareRAA=all/*all = show hadron RAA and jet RAA models in one plot*/, TString ext="pdf")
{
	if(ppBias==1 && ((pTlead>3.0 && pTlead<5.0)||(pTlead>5.0 && pTlead<7.0))) return; //we have pythia histograms only for pTlead 3, 5 and 7 GeV
	
	//setings
	//RAA
	bool showRAA=1; //show RAA plot at all
	bool showKrishnaAll=1; //show all 3 variants of hybrid model calculation
	bool showVitev=1; //show Ivan's RAA calculation
	bool showKrishna=1; //show Krishna's RAA calculation
		bool wake=1; //parameter of Krishna's calculation
		bool wake_positive_only=0; //parameter of Krishna's calculation
	
	bool showAliceRAA=0; //show Alice full jet RAA
	bool showPhenixhRAAchh=1; //show PHENIX's charged hadron RAA
	bool showPhenixhRAApi0=1; //show PHENIX's pi0 RAA
	bool showSTARRAAchh=1; //show STAR's charged hadron RAA
	bool plotLogX_RAA=0; //log x-axis
	bool plotLogY_RAA=1; //log y-axis
	bool drawBiasLine=0; //mark unbiased region
	
	if(showKrishnaAll) //show all 3 variants of Krishna's calculation 
	{
		showKrishna=0;
	}
	if(compareRAA==hadrons) //show only hadrons RAA, not the theoretical models
	{
		showKrishnaAll=0;
		showVitev=0;
		showKrishna=0;
		plotLogX_RAA=1;
	}
	else if(compareRAA==theory) //show only theoretical models, not the hadron RAA
	{
		showAliceRAA=0;
		showPhenixhRAAchh=0;
		showPhenixhRAApi0=0;
		showSTARRAAchh=0;
	}
	
	//RCP
	bool showRCP=1; //show RCP plot at all
	bool showAlex=0; //show Alex's h+jet I_CP
	bool showAliceRCP=1; //show Alice RCP
	bool showSTARhRCP=1; //show STAR's charged hadron RCP
	bool showAtlasRCP=1; //show ATLAS ch. hadron RCP
	bool plotLogX_RCP=1; //log x-axis
	bool plotLogY_RCP=1; //log y-axis

	if(ppBias || system==peri)showRCP=0; //show only RAA since RCP is the same as with ppBias=0
	if(system==peri)
	{
		showVitev=0; 
	 showKrishna=0; 
	}
	//Colors
	Color_t cstar=kBlue; //STAR data
	Color_t cgray=kGray+1; //STAR data in biased region
	Color_t cstarH=kAzure+1; //star hadron RCP
	Color_t cstarJ=kCyan+1; //star jet norm. unc.
	Color_t calice=kRed; //ALICE
	Color_t catlas=kRed-9; //ATLAS
	Color_t calex=kBlue+1; //STAR h+jet
	Color_t cphx=kMagenta-6; //PHENIX
	Color_t ckrsh=kMagenta+2; //Krishna's HYBRID model
		Color_t ckrsh1=kMagenta+2; //Krishna's HYBRID model
		Color_t ckrsh2=kViolet-9; //Krishna's HYBRID model
		Color_t ckrsh3=kOrange-8; //Krishna's HYBRID model
	Color_t ctheory_vit=kCyan-2; //Ivan's SCET model
	Color_t ctheory_vitNLO=kAzure-7; //Ivan's NLO
	
	//line at unity
	Color_t lineOneColor=kGray;
	int lineOneWidth=1;
	int lineOneStyle=2;
	
	//pTlead bias line
	Color_t lineBiasColor=kRed-6;
	int lineBiasWidth=2;
	float xbiasLine[C_nR]={14,16,16}; //x position of a line showing pTlead bias range in RAA
	float ybiasLine=1.0;//y range of a line showing pTlead bias range in RAA
	float xbiasDesc[C_nR]={0.52,0.52,0.52}; //x position of the bias textbox
	float ybiasDesc=0.26; //x position of the bias textbox
	
	//descriptions
	float G_latex_sz=0.032;
	Color_t G_latex_cl=kGray+1;


	
	int markstar=29;
	int msizestar=4;
	int markalice=20;
	int msizealice=2;
	int markatlasH=24;
	int msizeatlasH=2;
	int markphx=27;
	int markphx2=24;
	int msizephx=2;
	int markstarH=30;
	int msizestarH=3;
	int phxFill=1001; //phenix fill style
	int krshFill=1001; //Krishna's RAA fill style
		int krshFill1=1001; //Krishna's RAA fill style
		int krshFill2=1001; //Krishna's RAA fill style
		int krshFill3=1001; //Krishna's RAA fill style
	
	//axis ranges
	//RAA
	float xmaxRAA=46;
	float xminRAA=(plotLogX_RAA)?4.0:0.0;
	float yminRAA=(plotLogY_RAA)? 0.05:0;
	float ymaxRAA=(plotLogY_RAA)?10.0:1.0; //2.0:0.8
	//RCP
	float xminRCP=0.05;
	float xmaxRCP=120;
	float xmaxRCP_short=40; //in case we don't show LHC results
	if(plotLogX_RCP) 
	{
		xminRCP=0.2;
		xmaxRCP=150;
		xmaxRCP_short=100; 
	}

	float yminRCP=(plotLogY_RCP)? 0.05:0.05;
	float ymaxRCP=(plotLogY_RCP)?2.0:1.3;
	float ptminRCP_STAR=pTlead;//STAR data range
	float ptmaxRCP_STAR=20; //STAR data range

	
	gROOT->LoadMacro("./util.C");
	gROOT->LoadMacro("./Utility.C");
   gStyle->SetOptStat(0);
	
	TString outDir=Form("../../../plotting_out/obr/%s/results/comparison/",evo.Data());
	TString sufPPbias[2]={"","_ppBiased"};
	
	TCanvas *cfig1=new TCanvas("cfig1","fig1",10,10,1400,500);
	cfig1->Divide(C_nR,1,0,0);
	
	//Alice
	double alxRAA[3]={45,55,65};
	double alxRAA_err[3]={5,5,5};
	double alyRAA[3]={0.293,0.282,0.266};
	double alyRAA_err_shape[3]={0.015,0.01,0.015};
	double alyRAA_err_corr[3]={0.035,0.032,0.035};
	double alyRAA_err[3];
	for(int i=0;i<3;i++)
	{
		alyRAA_err[i]=TMath::Sqrt(alyRAA_err_corr[i]*alyRAA_err_corr[i]+alyRAA_err_shape[i]*alyRAA_err_shape[i]);
	}
	TGraphErrors* galiceRAA=new TGraphErrors(3,alxRAA,alyRAA,alxRAA_err,0);
	TGraphErrors* galiceRAA_sys=new TGraphErrors(3,alxRAA,alyRAA,alxRAA_err,alyRAA_err);

	//Phenix charged hadrons RAA
	/*
	double phxRAAchhad_x[5]={5,5.4,5.8,6.5,7.5};
	double phxRAAchhad_x_err[5]={0.2,0.2,0.2,0.5,0.5};
	double phxRAAchhad_y_high[5]={0.0450693909433,0.04783097845539,0.042008600309937,0.048012602512257,0.081288513948774};
	double phxRAAchhad_y_center[5]={0.23,0.237,0.189,0.226,0.209};
	double phxRAAchhad_y_low[5]={0.039385911186616,0.044640816524791,0.035884850563992,0.050923570181204,0.066811844009876};
	*/
	double phxRAAchhad_x_0to10[31]={0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95,2.1,2.3,2.5,2.7,2.9,3.125,3.375,3.625,3.875,4.2,4.6,5,5.4,5.8,6.5,7.5};
	double phxRAAchhad_y_0to10[31]={0.2508,0.2755,0.3001,0.3263,0.3525,0.3824,0.4133,0.4466,0.4765,0.506,0.5287,0.5484,0.5693,0.5841,0.5858,0.5853,0.574,0.5538,0.5236,0.4909,0.4465,0.4077,0.3715,0.3402,0.3119,0.2623,0.2417,0.2452,0.1978,0.2316,0.2214};
	double phxRAAchhad_y_stat_0to10[31]={0.0002,0.0002,0.0002,0.0003,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0010,0.0012,0.0014,0.0016,0.0018,0.0015,0.0020,0.0024,0.0029,0.0035,0.0038,0.0046,0.0055,0.0067,0.0066,0.0087,0.0117,0.0166,0.0220,0.0282,0.0629};
	double phxRAAchhad_y_sysUp_0to10[31]={0.045,0.050,0.054,0.058,0.062,0.066,0.070,0.074,0.077,0.079,0.080,0.079,0.079,0.079,0.077,0.074,0.070,0.066,0.062,0.060,0.056,0.052,0.048,0.044,0.041,0.035,0.033,0.034,0.030,0.040,0.044};
	double phxRAAchhad_y_sysDown_0to10[31]={0.031,0.034,0.037,0.040,0.043,0.047,0.050,0.054,0.057,0.061,0.063,0.065,0.067,0.069,0.068,0.068,0.066,0.063,0.059,0.055,0.051,0.046,0.042,0.038,0.035,0.030,0.028,0.029,0.025,0.034,0.038};

	double phxRAAchhad_x_60to80[30]={0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95,2.1,2.3,2.5,2.7,2.9,3.125,3.375,3.625,3.875,4.2,4.6,5,5.4,5.8,6.5};
	double phxRAAchhad_y_60to80[30]={0.56425,0.5876,0.6096,0.63255,0.6546,0.68435,0.70665,0.7356,0.7649,0.78365,0.822,0.82795,0.85095,0.8725,0.88655,0.8816,0.90985,0.9329,0.9231,0.89105,0.9143,0.8585,0.8937,0.82065,0.85855,0.8203, 0.73315,0.61165,0.76395,0.6206};
	double phxRAAchhad_y_stat_60to80[30]={0.0010,0.0012,0.0014,0.0017,0.0020,0.0024,0.0028,0.0033,0.0039,0.0047,0.0055,0.0064,0.0074,0.0085,0.0098,0.0084,0.0108,0.0138,0.0171,0.0207,0.0234,0.0291,0.0369,0.0461,0.0466,0.0618,0.0837,0.1295,0.1664,0.1694};
	double phxRAAchhad_y_sysUp_60to80[30]={0.099,0.103,0.106,0.109,0.112,0.115,0.116,0.118,0.119,0.117,0.118,0.114,0.113,0.111,0.109,0.104,0.104,0.104,0.102,0.101,0.108,0.102,0.107,0.100,0.106,0.104,0.095,0.082,0.110,0.103};
	double phxRAAchhad_y_sysDown_60to80[30]={0.067,0.069,0.072,0.074,0.076,0.079,0.082,0.085,0.087,0.089,0.093,0.093,0.095,0.097,0.098,0.097,0.099,0.101,0.099,0.094,0.098,0.092,0.095,0.087,0.091,0.088,0.080,0.069,0.093,0.088};
	
	TGraphAsymmErrors *gphxRAAchhad[2];
	TGraphAsymmErrors *gphxRAAchhad_sys[2];
	gphxRAAchhad[0]=new TGraphAsymmErrors(31,phxRAAchhad_x_0to10,phxRAAchhad_y_0to10,0,0,phxRAAchhad_y_stat_0to10,phxRAAchhad_y_stat_0to10);
	gphxRAAchhad[1]=new TGraphAsymmErrors(30,phxRAAchhad_x_60to80,phxRAAchhad_y_60to80,0,0,phxRAAchhad_y_stat_60to80,phxRAAchhad_y_stat_60to80);
	gphxRAAchhad_sys[0]=new TGraphAsymmErrors(31,phxRAAchhad_x_0to10,phxRAAchhad_y_0to10,0,0,phxRAAchhad_y_sysDown_0to10,phxRAAchhad_y_sysUp_0to10);
	gphxRAAchhad_sys[1]=new TGraphAsymmErrors(30,phxRAAchhad_x_60to80,phxRAAchhad_y_60to80,0,0,phxRAAchhad_y_sysDown_60to80,phxRAAchhad_y_sysUp_60to80);
	
	
	//Phenix pi0 RAA from arXiv:1208.2254 
	double phxRAApi0_x_0to10[15]={5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 11, 13, 15, 17, 19};
	double phxRAApi0_y_0to10[15]={0.18594, 0.18552, 0.18822, 0.19223, 0.19077, 0.19671, 0.19897, 0.197, 0.22411, 0.22495, 0.22526, 0.24026, 0.32436, 0.3763, 0.26386};
	double phxRAApi0_y_stat_0to10[15]={0.0013671, 0.0017548, 0.0022771, 0.0029342, 0.0036067, 0.0045438, 0.0056201, 0.0066363, 0.0088675, 0.010915, 0.0076901, 0.015137, 0.036734, 0.072167, 0.099973};
	double phxRAApi0_y_sys_0to10[15]={0.023779, 0.023779, 0.024178, 0.024751, 0.024628, 0.02547, 0.024184, 0.024053, 0.027549, 0.028444, 0.033706, 0.05017, 0.090579, 0.13431, 0.11517};
	
	double phxRAApi0_x_60to80[14]={5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 11, 13, 15, 17};
	double phxRAApi0_y_60to80[14]={0.82138, 0.82911, 0.818905, 0.8443, 0.83944, 0.82398, 0.82851, 0.85778, 0.78218, 0.98742, 0.73433, 1.001065, 0.951795, 1.7073};
	double phxRAApi0_y_stat_60to80[14]={0.0116605, 0.0159, 0.021306, 0.028004, 0.0360805, 0.046562, 0.055896, 0.0677415, 0.08523, 0.1115945, 0.067189, 0.148925, 0.28749, 0.76182};
	double phxRAApi0_y_sys_60to80[14]={0.1050385, 0.106275, 0.1051945, 0.10871, 0.1083715, 0.106691, 0.100705, 0.104728, 0.09615, 0.12485, 0.10988, 0.209045, 0.26579, 0.60936};
	
	
	TGraphErrors *gphxRAApi0[2];
	TGraphErrors *gphxRAApi0_sys[2];
	gphxRAApi0[0]=new TGraphErrors(15,phxRAApi0_x_0to10,phxRAApi0_y_0to10,0,phxRAApi0_y_stat_0to10);
	gphxRAApi0_sys[0]=new TGraphErrors(15,phxRAApi0_x_0to10,phxRAApi0_y_0to10,0,phxRAApi0_y_sys_0to10);
	gphxRAApi0[1]=new TGraphErrors(14,phxRAApi0_x_60to80,phxRAApi0_y_60to80,0,phxRAApi0_y_stat_60to80);
	gphxRAApi0_sys[1]=new TGraphErrors(14,phxRAApi0_x_60to80,phxRAApi0_y_60to80,0,phxRAApi0_y_sys_60to80);
	
	//TGraphAsymmErrors *gphxRAA=new TGraphAsymmErrors(5,phxRAAchhad_x,phxRAAchhad_y_center,phxRAAchhad_x_err,phxRAAchhad_x_err,phxRAAchhad_y_low,phxRAAchhad_y_high);
	//TGraph *gphxRAAchhad_sys=graph_band(5,phxRAAchhad_x,phxRAAchhad_y_high,phxRAAchhad_y_low,50);
	
	//STAR charged hadrons RAA
	/*double strRAA_x[32]={0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.49, 2.7, 2.9, 3.16, 3.55, 4.06, 4.7, 5.48, 6.42, 7.43, 8.43, 9.44};
	double strRAA_y_0to5[32]={0.29, 0.32, 0.36, 0.39, 0.43, 0.46, 0.49, 0.53, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.67, 0.67, 0.69, 0.68, 0.7, 0.69, 0.7, 0.65, 0.56, 0.49, 0.44, 0.38, 0.32, 0.27, 0.26, 0.18, 0.22, 0.19};
	double strRAA_y_0to5_err[32]={0.03, 0.03, 0.03, 0.03, 0.03, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.06, 0.06, 0.06, 0.07, 0.08, 0.08, 0.08, 0.09, 0.1, 0.09, 0.06, 0.05, 0.04, 0.04, 0.04, 0.03, 0.04, 0.03, 0.06, 0.06};
	double strRAA_y_60to80[32]={0.7, 0.74, 0.77, 0.81, 0.85, 0.88, 0.91, 0.93, 0.96, 0.98, 0.98, 0.99, 1.03, 1.06, 1.08, 1.11, 1.17, 1.18, 1.23, 1.29, 1.36, 1.38, 1.28, 1.2, 1.22, 1.21, 1.19, 1.19, 1.24, 0.95, 0.82, 0.78};
	double strRAA_y_60to80_err[32]={0.06, 0.06, 0.06, 0.06, 0.07, 0.07, 0.07, 0.07, 0.07, 0.08, 0.08, 0.08, 0.09, 0.1, 0.1, 0.11, 0.13, 0.13, 0.15, 0.17, 0.2, 0.19, 0.15, 0.11, 0.12, 0.13, 0.14, 0.16, 0.22, 0.24, 0.34, 0.42};*/
	double strRAA_x[18]={0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.25, 2.7, 2.9, 3.16, 3.55, 4.06, 4.7, 5.48, 6.42, 7.43, 8.43, 9.44};
	double strRAA_y_0to5[18]={0.29, 0.39, 0.49, 0.58, 0.64, 0.67, 0.7, 0.65,0.56, 0.49, 0.44, 0.38, 0.32, 0.27, 0.26, 0.18, 0.22, 0.19};
	double strRAA_y_0to5_err[18]={0.03, 0.03, 0.04, 0.04, 0.06, 0.07, 0.08, 0.09, 0.06, 0.05, 0.04, 0.04, 0.04, 0.03, 0.04, 0.03, 0.06, 0.06};
	double strRAA_y_60to80[18]={0.7, 0.81, 0.91, 0.98, 1.03, 1.11, 1.23, 1.38, 1.28, 1.2, 1.22, 1.21, 1.19, 1.19, 1.24, 0.95, 0.82, 0.78};
	double strRAA_y_60to80_err[18]={0.06, 0.06, 0.07, 0.08, 0.09, 0.11, 0.15, 0.19, 0.15, 0.11, 0.12, 0.13, 0.14, 0.16, 0.22, 0.24, 0.34, 0.42};
	//normalization errors
	double normErrStr[2]={0.14,0.18}; //STAR hadrons
	double normErrStrJet[2]={0.21,0.35}; //STAR jets
	double normErrPhx[2]={0.12,0.22}; //Phenix hadrons
	
	TGraphErrors *gstrRAA[3];
	gstrRAA[0]=new TGraphErrors(18,strRAA_x,strRAA_y_0to5,0,strRAA_y_0to5_err);
	gstrRAA[1]=new TGraphErrors(18,strRAA_x,strRAA_y_60to80,0,strRAA_y_60to80_err);
	
	//Alex Icp
	double icp_R02_x[9]={10.0,11.9,13.3,15.8,17.339,18.6,20.0,22.647,25.3};
	double icp_R02_y_low[9]={0.174,0.189,0.234,0.206,0.182,0.160,0.204,0.1376,0.169};
	double icp_R02_y_high[9]={0.49,0.397,0.426,0.44,0.463,0.406,0.574,0.420,0.84};
	double icp_R02_y_center[9]={0.313,0.285,0.312,0.302,0.293,0.269,0.340,0.242,0.306};
	double icp_R03_x[11]={10.0,10.559,12.527,15.496,16.529,19.5,20.516,21.578,22.837,24.53,25.511};
	double icp_R03_y_low[11]={0.165,0.149,0.118,0.142,0.173,0.17,0.188,0.15,0.157,0.133,0.133};
	double icp_R03_y_high[11]={0.358,0.313,0.307,0.399,0.519,0.610,0.648,0.574,0.628,0.793,0.724};
	double icp_R03_y_center[11]={0.244,0.225,0.201,0.247,0.305,0.3,0.351,0.3,0.302,0.264,0.287};

	
	TGraph *gIcp[C_nR];
	TGraph *gIcp_sys[C_nR];
	gIcp[0]=new TGraph(9,icp_R02_x,icp_R02_y_center);
	gIcp[1]=new TGraph(11,icp_R03_x,icp_R03_y_center);
	gIcp_sys[0]=graph_band(9,icp_R02_x,icp_R02_y_high,icp_R02_y_low,50);
	gIcp_sys[1]=graph_band(11,icp_R03_x,icp_R03_y_high,icp_R03_y_low,50);
	
	//ATLAS ch. hadrons RCP 	arXiv:1504.04337
	const int nATLAShRCPbins=35;
	double atlashRCP_x[]={0.5365,0.615,0.705,0.808,0.926,1.0595,1.21,1.385,1.59,1.825,2.095,2.405,2.755,3.155,3.62,4.15,
	4.755,5.455,6.255,7.17,8.22,9.39,10.75,12.35,14.15,16.2,18.6,21.35,24.45,28.05,33.85,42.6,53.65,67.55,85.05};
	double atlashRCP_xw[]={0.0365,0.042,0.048,0.055,0.063,0.0705,0.08,0.095,0.11,0.125,0.145,0.165,
		0.185,0.215,0.25,0.28,0.325,0.375,0.425,0.49,0.56,0.61,0.75,0.85,0.95,1.1,1.3,1.45,1.65,1.95,3.85,4.9,6.15,7.75,9.75};
	double atlashRCP_y[]={0.464,0.481,0.499,0.517,0.541,0.557,0.575,0.591,0.605,0.603,0.599,0.575,0.53,0.468,
	0.399,0.33,0.275,0.24,0.221,0.221,0.23,0.242,0.26,0.283,0.305,0.332,0.358,0.402,0.428,0.436,0.523,0.576,0.667,0.677,0.606};
	double atlashRCP_stat[]={4.59E-05,4.18E-05,4.49E-05,4.91E-05,5.41E-05,6.13E-05,6.90E-05,8.27E-05,9.68E-05,
	1.15E-04,1.38E-04,1.61E-04,1.96E-04,2.20E-04,2.51E-04,2.84E-04,3.30E-04,3.84E-04,
	4.86E-04,6.41E-04,8.97E-04,1.26E-03,1.22E-03,1.78E-03,2.53E-03,3.65E-03,5.01E-03,7.24E-03,
	9.42E-03,1.13E-02,1.36E-02,1.67E-02,2.67E-02,4.67E-02,7.03E-02};
	double atlashRCP_sys[]={0.057072,0.058682,0.060379,0.062557,0.06492,0.06684,0.068425,0.070329,0.071995,0.071757,
	0.071281,0.068425,0.06307,0.05616,0.04788,0.0396,0.033,0.0288,0.02652,0.02652,0.0276,0.02904,0.0312,
	0.03396,0.036295,0.039508,0.042602,0.04824,0.05136,0.052756,0.064852,0.072,0.08671,0.091395,0.086658};
	TGraphErrors* gatlashRCP=new TGraphErrors(nATLAShRCPbins,atlashRCP_x,atlashRCP_y,0,atlashRCP_stat);
	TGraphErrors* gatlashRCP_sys=new TGraphErrors(nATLAShRCPbins,atlashRCP_x,atlashRCP_y,atlashRCP_xw,atlashRCP_sys);
	
	//Alice RCP (arXiv:1311.0633)
	TGraphErrors* galiceRCP[C_nR];
	TGraphAsymmErrors* galiceRCP_sys[C_nR];
	const int nAliceRCPbins=7;
	
	double alxRCP02[nAliceRCPbins]={25,35,45,55,65,75,85};
	double alxRCP02_err[nAliceRCPbins]={5,5,5,5,5,5,5};
	double alyRCP02[nAliceRCPbins]={0.325,0.3186,0.4175,0.5223,0.5334,0.4747,0.4032};
	double alyRCP02_err[nAliceRCPbins]={0.033,0.0441,0.06,0.0792,0.0972,0.1195,0.141};
	double alyRCP02_syserr_up[nAliceRCPbins]={0.034,0.041,0.041,0.042,0.079,0.124,0.136};
	double alyRCP02_syserr_down[nAliceRCPbins]={0.076,0.043,0.055,0.116,0.068,0.105,0.027};
	/*
	double alxRCP02[3]={25,35,45};
	double alxRCP02_err[3]={5,5,5};
	double alyRCP02[3]={0.3254,0.318,0.418};
	double alyRCP02_err[3]={0.033,0.034,0.066};
	double alyRCP02_syserr_up[3]={0.034,0.049,0.034};
	double alyRCP02_syserr_down[3]={0.076,0.041,0.034};*/
	
	galiceRCP[0]=new TGraphErrors(nAliceRCPbins,alxRCP02,alyRCP02,alxRCP02_err,alyRCP02_err);
	galiceRCP_sys[0]=new TGraphAsymmErrors(nAliceRCPbins,alxRCP02,alyRCP02,alxRCP02_err,alxRCP02_err,alyRCP02_syserr_up,alyRCP02_syserr_down);
	
	double alxRCP03[nAliceRCPbins-2]={35,45,55,65,75};
	double alxRCP03_err[nAliceRCPbins-2]={5,5,5,5,5};
	double alyRCP03[nAliceRCPbins-2]={0.352,0.418,0.5055,0.5133,0.4549};
	double alyRCP03_err[nAliceRCPbins-2]={0.072,0.073,0.1164,0.1091,0.1331};
	double alyRCP03_syserr_up[nAliceRCPbins-2]={0.057,0.053,0.086,0.054,0.054};
	double alyRCP03_syserr_down[nAliceRCPbins-2]={0.118,0.066,0.165,0.149,0.118};
	
	/*
	double alxRCP03[2]={35,45};
	double alxRCP03_err[2]={5,5};
	double alyRCP03[2]={0.352,0.418};
	double alyRCP03_err[2]={0.072,0.073};
	double alyRCP03_syserr_up[2]={0.0566,0.0533};
	double alyRCP03_syserr_down[2]={0.117,0.0658};*/
	
	galiceRCP[1]=new TGraphErrors(nAliceRCPbins-2,alxRCP03,alyRCP03,alxRCP03_err,alyRCP03_err);
	galiceRCP_sys[1]=new TGraphAsymmErrors(nAliceRCPbins-2,alxRCP03,alyRCP03,alxRCP02_err,alxRCP02_err,alyRCP03_syserr_up,alyRCP03_syserr_down);
	
	//STAR hadron RCP (arxiv:nucl-ex/0305015)
	const int nhRCPbins=35;
	double hRCPx[nhRCPbins]={0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,1.55,
		1.65,1.75,1.85,1.95,2.05,2.15,2.25,2.35,2.49,2.70,2.90,3.16,3.55,4.06,4.70,5.48,6.42,7.43,8.43,9.44,10.78};
	double hRCPy[nhRCPbins]={0.37,0.4,0.42,0.44,0.47,0.49,0.51,0.52,0.55,0.57,0.59,0.6,0.62,0.63,0.63,
		0.63,0.62,0.61,0.6,0.58,0.57,0.54,0.52,0.47,0.44,0.41,0.36,0.32,0.27,0.23,0.21,0.19,0.27,0.24,0.26};
	double hRCPerr[nhRCPbins]={0.02,0.02,0.02,0.02,0.02,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,
		0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.02,0.02,0.02,0.02,0.02,0.01,0.02,0.02,0.04,0.1,0.12,0.15};
	TGraphErrors* gSTARhRCP=new TGraphErrors(nhRCPbins,hRCPx,hRCPy,0,hRCPerr);
	
	
	//declare arrays
	TFile* fp[C_nR];

	TGraphErrors *gVitRAA_g22[C_nR];
	TGraphErrors *gVitRAA_g20[C_nR];
	TGraphErrors *gVitRAA_g22b[C_nR];
	TGraphErrors *gVitRAA_g20b[C_nR];
	TGraphErrors *gVitRAA_NLO_down[C_nR];
	TGraphErrors *gVitRAA_NLO_up[C_nR];
			
	TGraphErrors *gKrishnaRAA[C_nR];
	TGraphErrors *gKrishnaRAA1[C_nR];
	TGraphErrors *gKrishnaRAA2[C_nR];
	TGraphErrors *gKrishnaRAA3[C_nR];

	TGraphAsymmErrors *gjanRAA[C_nR];
	TGraphAsymmErrors *gjanRAA_sys[C_nR];
	TGraphAsymmErrors *gjanRAA_gray[C_nR]; //gray data points in biased region
	TGraphAsymmErrors *gjanRAA_sys_gray[C_nR]; //gray data points in biased region
	TGraphAsymmErrors *gjanRCP[C_nR];
	TGraphAsymmErrors *gjanRCP_sys[C_nR];
	TH1D* htmp1[C_nR];
	TH1D* htmp2[C_nR];
	
   TLatex *latex = new TLatex();
   latex->SetNDC();
   latex->SetTextSize(0.05);
	float Rxpos[C_nR]={0.75,0.15,0.15};
	float Rxpos2[C_nR]={0.7,0.4,0.5};
	float shrinkKrishna[C_nR]={40,40,40};
	
	float infoXmin,infoYmin,infoXmax,infoYmax;
	float posx1,posy1,posx2,posy2;
	
	TLine *one = new TLine(xminRAA, 1, xmaxRAA,1);
	one->SetLineWidth(lineOneWidth);
	one->SetLineStyle(lineOneStyle);
	one->SetLineColor(lineOneColor);
	
	//STAR hadron normalization error
	double boxX=xminRAA+9*(xmaxRAA-xminRAA)/10;
	TBox* ppboxStr=new TBox(boxX-2, 1-normErrStr[system], boxX-1, 1+normErrStr[system]);
	ppboxStr->SetFillColor(cstarH);
	ppboxStr->SetFillStyle(1000);
	
	//PHENIX hadron normalization error
	TBox* ppboxPhx=new TBox(boxX-4, 1-normErrPhx[system], boxX-3, 1+normErrPhx[system]);
	ppboxPhx->SetFillColor(cphx);
	ppboxPhx->SetFillStyle(1000);
	
	//STAR jet normalization error
	TBox* ppboxStrJet=new TBox(boxX, 1-normErrStrJet[system], boxX+1, 1+normErrStrJet[system]);
	ppboxStrJet->SetFillColor(cstarJ);
	ppboxStrJet->SetFillStyle(1000);
	
	//********************
	//RAA
	//********************
	if(showRAA){
	//loop over R
	for(int r=0; r<C_nR; r++)
	{

	float R=C_Rs[r];
	bool plotAlice=(R<0.25) ? true : false; //we have Alice results only for R=0.2

	//---------------------------------------------------------------------------
  //Vitev
	double cronin1=1.5;
	double eloss1=1.0;
	double cronin2=1.0;
	double eloss2=1.5;
	Color_t cvit1=ctheory_vit; //line color
	Color_t cvit2=ctheory_vit;
	Color_t cvit3=ctheory_vitNLO;
	int lsvit1=1; //line style
	int lsvit2=1;
	int lsvit3=3;
	int lwvit1=1; //line width
	int lwvit2=1;
	int lwvit3=1; 
	
	//SCET model
	gVitRAA_g22[r]=gReadData(Form("theory/vitev/AuAu_cron%.1lf_eloss%.1lf_200GeV_g2.2_R%0.lf.txt",cronin1,eloss1,R*10),1,0);
	gVitRAA_g20[r]=gReadData(Form("theory/vitev/AuAu_cron%.1lf_eloss%.1lf_200GeV_g2.0_R%0.lf.txt",cronin1,eloss1,R*10),1,0);
	
	gVitRAA_g22b[r]=gReadData(Form("theory/vitev/AuAu_cron%.1lf_eloss%.1lf_200GeV_g2.2_R%0.lf.txt",cronin2,eloss2,R*10),1,0);
	gVitRAA_g20b[r]=gReadData(Form("theory/vitev/AuAu_cron%.1lf_eloss%.1lf_200GeV_g2.0_R%0.lf.txt",cronin2,eloss2,R*10),1,0);
	
	gVitRAA_g20[r]->SetLineColor(cvit1);
	gVitRAA_g20[r]->SetLineStyle(lsvit1);
	gVitRAA_g20[r]->SetLineWidth(lwvit1);
	gVitRAA_g20b[r]->SetLineColor(cvit2);
	gVitRAA_g20b[r]->SetLineStyle(lsvit2);
	gVitRAA_g20b[r]->SetLineWidth(lwvit2);
	
	
	//NLO 
	if(r!=1)
	{
		gVitRAA_NLO_up[r]=gReadData(Form("theory/vitev/Raa_NLO-Au-Wmin0.R%.1lfcoldUP.dat",R),1,0);
		gVitRAA_NLO_down[r]=gReadData(Form("theory/vitev/Raa_NLO-Au-Wmin0.R%.1lfcoldDN.dat",R),1,0);
	
		gVitRAA_NLO_down[r]->SetLineColor(cvit3);
		gVitRAA_NLO_down[r]->SetLineStyle(lsvit3);
		gVitRAA_NLO_down[r]->SetLineWidth(lwvit3);
	}
	
	//Krishna 
	//gKrishnaRAA[r]=gReadData3(Form("theory/hybrid_results_2017/cent010jetRAA_R0%.0lf.dat",R*10),5);
	gKrishnaRAA[r]=gReadData4(Form("theory/hybrid_results_2019/010_RAA_R%.0lf_wake_%i_ignore_neg_%i.dat",R*10,wake,wake_positive_only),5);
	if(showKrishnaAll)
	{
		gKrishnaRAA1[r]=gReadData4(Form("theory/hybrid_results_2019/010_RAA_R%.0lf_wake_0_ignore_neg_0.dat",R*10),5);
		gKrishnaRAA2[r]=gReadData4(Form("theory/hybrid_results_2019/010_RAA_R%.0lf_wake_1_ignore_neg_1.dat",R*10),5);
		gKrishnaRAA3[r]=gReadData4(Form("theory/hybrid_results_2019/010_RAA_R%.0lf_wake_1_ignore_neg_0.dat",R*10),5);
	}
	
	//Jan
	TString indir=Form("%s%s",evo.Data(),sufPPbias[ppBias].Data());
	fp[r]=new TFile(Form("../../../plotting_out/systematics/%s/%s/bining%i/systematics_R%.1lf_pT%0.lf.root",C_system_name[system].Data(),indir.Data(),binning,R,pTlead),"read");
	TGraphAsymmErrors *gjanRAA_unf=(TGraphAsymmErrors*) fp[r]->Get("RAA_unfold_BSL");
	TGraphAsymmErrors *gjanRAA_corr=(TGraphAsymmErrors*) fp[r]->Get("RAA_corr_BSL");
	gjanRAA_sys[r]=gcombine_syserr(gjanRAA_unf,gjanRAA_corr,1,1);
	gjanRAA[r]=(TGraphAsymmErrors*) fp[r]->Get("RAA_BSL");

	gjanRAA_sys_gray[r]=(TGraphAsymmErrors*) gjanRAA_sys[r]->Clone();
	gjanRAA_gray[r]=(TGraphAsymmErrors*) gjanRAA[r]->Clone();
	
	//shrink graphs
	double shirnk_x=xbiasLine[r];
	if(ppBias==1) shirnk_x=pTlead;
	ShrinkGraph(gjanRAA[r],shirnk_x,30);
	ShrinkGraph(gjanRAA_sys[r],shirnk_x,30);
	ShrinkGraph(gjanRAA_gray[r],pTlead,shirnk_x);
	ShrinkGraph(gjanRAA_sys_gray[r],pTlead,shirnk_x);
	ShrinkGraph(gKrishnaRAA[r],pTlead,shrinkKrishna[r]);
	if(showKrishnaAll)
	{
		ShrinkGraph(gKrishnaRAA1[r],pTlead,shrinkKrishna[r]);
		ShrinkGraph(gKrishnaRAA2[r],pTlead,shrinkKrishna[r]);
		ShrinkGraph(gKrishnaRAA3[r],pTlead,shrinkKrishna[r]);
	}
	

	TLatex *latexL = new TLatex();
	latexL->SetNDC();
	latexL->SetTextSize(G_latex_sz);
	latexL->SetTextColor(G_latex_cl);
	
	//********************
	//Draw plots
	//********************
  TString xtitle="";
  if(r==C_nR-1) xtitle="p_{T, jet}^{ch}, p_{T, jet} (GeV/#it{c})";
  //if(!plotAlice) xtitle="p_{T, jet}^{ch} [GeV/#it{c}]";
  cfig1->cd(r+1);
  gPad->SetMargin(0.20,0.03,0.20,0.03);
  if(r==0) gPad->SetRightMargin(0);
  if(r==1){
	  gPad->SetLeftMargin(0);
	  gPad->SetRightMargin(0);
	  xminRAA=xminRAA+0.1;}
  if(r==2) gPad->SetLeftMargin(0);
  gPad->SetTicks(1);
  if(plotLogY_RAA) gPad->SetLogy();
  if(plotLogX_RAA) gPad->SetLogx();
  htmp1[r]=new TH1D(Form("htmp1_%i",r),"",100,xminRAA,xmaxRAA);
  htmp1[r]->SetXTitle(xtitle);htmp1[r]->SetYTitle("R_{AA}^{Pythia}");
  htmp1[r]->SetTitleOffset(1.2,"x");htmp1[r]->SetTitleOffset(1.1,"y");
  htmp1[r]->SetTitleSize(0.075,"x");htmp1[r]->SetTitleSize(0.08,"y");
  htmp1[r]->SetLabelSize(0.055,"x");htmp1[r]->SetLabelSize(0.055,"y");
  //htmp1[r]->SetNdivisions(505,"x");htmp1[r]->SetNdivisions(505,"y");
  htmp1[r]->SetMinimum(yminRAA);htmp1[r]->SetMaximum(ymaxRAA);
  htmp1[r]->DrawCopy();
  
    //Krishna's RAA
	gKrishnaRAA[r]->SetLineColor(ckrsh);
	gKrishnaRAA[r]->SetFillStyle(krshFill);
	gKrishnaRAA[r]->SetFillColor(ckrsh);
	if(showKrishna) gKrishnaRAA[r]->Draw("2");
	
	if(showKrishnaAll)
	{	
		gKrishnaRAA1[r]->SetLineColor(ckrsh1);
		gKrishnaRAA1[r]->SetFillStyle(krshFill1);
		gKrishnaRAA1[r]->SetFillColor(ckrsh1);
		gKrishnaRAA1[r]->Draw("2");
		
		gKrishnaRAA2[r]->SetLineColor(ckrsh2);
		gKrishnaRAA2[r]->SetFillStyle(krshFill2);
		gKrishnaRAA2[r]->SetFillColor(ckrsh2);
		gKrishnaRAA2[r]->Draw("2");
	
		gKrishnaRAA3[r]->SetLineColor(ckrsh3);
		gKrishnaRAA3[r]->SetFillStyle(krshFill3);
		gKrishnaRAA3[r]->SetFillColor(ckrsh3);
		gKrishnaRAA3[r]->Draw("2");
	}
 
   //Ivans's RAA
  double tsize=0.038;
  if(showVitev) plot_g2(gVitRAA_g20[r],gVitRAA_g22[r],cvit1,lsvit1,lwvit1,"",0,-1,tsize);
  //if(showVitev) plot_g2(gVitRAA_g20b[r],gVitRAA_g22b[r],cvit2,lsvit2,lwvit2,"",0,-1,tsize);
  if(showVitev && r!=1) plot_g2(gVitRAA_NLO_down[r],gVitRAA_NLO_up[r],cvit3,lsvit3,lwvit3,"",0,-1,tsize);
  
  //Phenix charged hadron RAA
	gphxRAAchhad[system]->SetLineColor(cphx);
   gphxRAAchhad[system]->SetLineWidth(1);
	gphxRAAchhad[system]->SetMarkerStyle(markphx);
	gphxRAAchhad[system]->SetMarkerSize(msizephx);
	gphxRAAchhad[system]->SetMarkerColor(cphx);
	gphxRAAchhad_sys[system]->SetLineColor(cphx);
   gphxRAAchhad_sys[system]->SetLineWidth(1);
	gphxRAAchhad_sys[system]->SetFillStyle(0);
	if(showPhenixhRAAchh) {
		gphxRAAchhad_sys[system]->Draw("[]");
		gphxRAAchhad[system]->Draw("P");
	}
	//Phenix pi0 RAA
	gphxRAApi0[system]->SetLineColor(cphx);
   gphxRAApi0[system]->SetLineWidth(1);
	gphxRAApi0[system]->SetMarkerStyle(markphx2);
	gphxRAApi0[system]->SetMarkerSize(msizephx);
	gphxRAApi0[system]->SetMarkerColor(cphx);
	gphxRAApi0_sys[system]->SetLineColor(cphx);
   gphxRAApi0_sys[system]->SetLineWidth(1);
	gphxRAApi0_sys[system]->SetFillStyle(0);
	if(showPhenixhRAApi0) {
		gphxRAApi0_sys[system]->Draw("[]");
		gphxRAApi0[system]->Draw("P");
	}
	
	//STAR hadron RAA
	gstrRAA[system]->SetLineWidth(1);
	gstrRAA[system]->SetLineColor(cstarH);
	gstrRAA[system]->SetMarkerStyle(markstarH);
	gstrRAA[system]->SetMarkerSize(msizestarH);
	gstrRAA[system]->SetMarkerColor(cstarH);
	if(showSTARRAAchh) gstrRAA[system]->Draw("P");
   
  //STAR RAA
  gjanRAA_sys[r]->SetLineWidth(1);
  gjanRAA_sys[r]->SetFillStyle(0);
  gjanRAA_sys[r]->SetLineColor(cstar);
  gjanRAA[r]->SetLineWidth(2);
  gjanRAA[r]->SetLineColor(cstar);
  gjanRAA[r]->SetMarkerStyle(markstar);
  gjanRAA[r]->SetMarkerSize(msizestar);
  gjanRAA[r]->SetMarkerColor(cstar);
  gjanRAA_sys_gray[r]->SetLineWidth(1);
  gjanRAA_sys_gray[r]->SetFillStyle(0);
  gjanRAA_sys_gray[r]->SetLineColor(cgray);
  gjanRAA_gray[r]->SetLineWidth(2);
  gjanRAA_gray[r]->SetLineColor(cgray);
  gjanRAA_gray[r]->SetMarkerStyle(markstar);
  gjanRAA_gray[r]->SetMarkerSize(msizestar);
  gjanRAA_gray[r]->SetMarkerColor(cgray);
  
  gjanRAA[r]->Draw("P");
  gjanRAA_sys[r]->Draw("2");
  if(drawBiasLine)
  {
		gjanRAA_gray[r]->Draw("P");
		gjanRAA_sys_gray[r]->Draw("2");
  }
  
  //ALICE RAA
  galiceRAA_sys->SetLineWidth(1);
  galiceRAA_sys->SetLineColor(calice);
  galiceRAA_sys->SetFillStyle(0);
  galiceRAA->SetLineWidth(2);
  galiceRAA->SetLineColor(calice);
  galiceRAA->SetMarkerStyle(markalice);
  galiceRAA->SetMarkerSize(msizealice);
  galiceRAA->SetMarkerColor(calice);
  if(plotAlice && showAliceRAA) galiceRAA->Draw("P");
  if(plotAlice && showAliceRAA) galiceRAA_sys->Draw("2");

 
	one->DrawClone("same");
	
	//normalization error
	if(showSTARRAAchh) ppboxStr->DrawClone("");
	if(showPhenixhRAAchh) ppboxPhx->DrawClone("");
	ppboxStrJet->DrawClone("");
  
  //LEGENDS
  if(r==0){
	  //general info
	posx1=0.25; posx2=0.6; posy1=0.65; posy2=0.9;
	if(system==peri){posx1=0.25;	posx2=0.6;	posy1=0.25;	posy2=0.45;}
	TLegend *model_info = new TLegend(posx1,posy1,posx2,posy2);
	model_info->SetTextSize(0.045);
	model_info->SetFillStyle(0);
	model_info->SetBorderSize(0);
	model_info->SetMargin(0.05);
	model_info->SetHeader("Au+Au #sqrt{s_{NN}}=200 GeV");
	//model_info->AddEntry("", Form("Run11, %s","MB"), "");
	//model_info->AddEntry("", "Charged jets", "");
	model_info->AddEntry("", "Central (0-10%)", "");
	model_info->AddEntry("","anti-k_{T}", "");
	//model_info->AddEntry("", Form("p_{T}^{leading} > %.1lf GeV/#it{c}",pTlead), "");
	model_info->DrawClone("same");
	}
 	else if(r==1)
	{
  //legend 1
	posx1=0.36;	posx2=0.95;	posy1=0.58;	posy2=0.92;
	if(compareRAA!=all) posy1=0.63;
	if(system==peri){posx1=0.05;	posx2=0.5;	posy1=0.25;	posy2=0.45;}
	TLegend *legraa = new TLegend(posx1,posy1,posx2,posy2);
	legraa->SetTextSize(0.038);
	legraa->SetFillStyle(0);
	legraa->SetBorderSize(0);
	//legraa->AddEntry(gjanRAA[r], "","");
	legraa->AddEntry(gjanRAA[r], "STAR charged jets","lp");
	legraa->AddEntry(gjanRAA[r], Form("  p_{T, lead}^{min} = %.0lf GeV/#it{c}",pTlead),"");
	legraa->AddEntry(ppboxStrJet, "      jet norm.unc.", "f");
	if(showAliceRAA)legraa->AddEntry(galiceRAA, "ALICE Pb+Pb (full j.)","lp");
	//legraa->AddEntry(gphxRAA, "PHENIX ch. hadrons ", "f");
	if(showSTARRAAchh)
	{
		legraa->AddEntry(gstrRAA[system], "STAR ch. had. (0-5%) ", "lp");
		legraa->AddEntry(ppboxStr, "      STAR ch.h. norm.unc.", "f");
	}
	if(compareRAA==all)
	{
		if(showPhenixhRAAchh) {
		legraa->AddEntry(gphxRAAchhad[system], "PHENIX ch. hadrons ", "lp");
		//legraa->AddEntry(ppboxPhx, "      hadron norm.unc.", "f");
		}	
		if(showPhenixhRAApi0) 
		legraa->AddEntry(gphxRAApi0[system], "PHENIX #pi^{0} ", "lp");
	}
	if(compareRAA==theory)
	{
		if(showVitev) legraa->AddEntry("","Full jets:","");
		if(showVitev) legraa->AddEntry(gVitRAA_g20[r], "  SCET", "l");
		//if(showVitev) legraa2->AddEntry(gVitRAA_g20b[r], "  SCET 2 ", "l");
		if(showVitev) legraa->AddEntry(gVitRAA_NLO_down[0], "  NLO pQCD ", "l");
	}
	legraa->DrawClone("same");
	
	latexL->DrawLatex(0.25, 0.7,label);
	}
	else if(r==2)
	{
		//legend 2
	posx1=0.35;	posx2=0.9;	posy1=0.58;	posy2=0.93;
	if(compareRAA==theory) posy1=0.65;
	else if(compareRAA==hadrons) posy1=0.75;
	if(system==peri){posx1=0.05;	posx2=0.55;	posy1=0.25;	posy2=0.45;}
	TLegend *legraa2 = new TLegend(posx1,posy1,posx2,posy2);
	legraa2->SetTextSize(0.04);
	legraa2->SetFillStyle(0);
	legraa2->SetBorderSize(0);
	if(compareRAA==hadrons)
	{
		if(showPhenixhRAAchh) 
			legraa2->AddEntry(gphxRAAchhad[system], "PHENIX ch. hadrons ", "lp");
		if(showPhenixhRAApi0) 
			legraa2->AddEntry(gphxRAApi0[system], "PHENIX #pi^{0} ", "lp");
		if(showPhenixhRAAchh || showPhenixhRAApi0)
			legraa2->AddEntry(ppboxPhx, "    PHENIX  had. norm.unc.", "f");

	}
	if(showKrishna) legraa2->AddEntry(gKrishnaRAA[r], "Hybrid model", "f");
	if(showKrishnaAll) 
	{
		legraa2->AddEntry("","Hybrid model, ch. jets:","");
		legraa2->AddEntry(gKrishnaRAA1[r], "   no medium resp.", "f");
		legraa2->AddEntry(gKrishnaRAA2[r], "   pos. resp. from wake ", "f");
		legraa2->AddEntry(gKrishnaRAA3[r], "   full medium resp.", "f");
	}
	if(compareRAA==all)
	{
		if(showVitev) legraa2->AddEntry("","Full jets:","");
		if(showVitev) legraa2->AddEntry(gVitRAA_g20[r], "  SCET", "l");
		//if(showVitev) legraa2->AddEntry(gVitRAA_g20b[r], "  SCET 2 ", "l");
		if(showVitev) legraa2->AddEntry(gVitRAA_NLO_down[r], "  NLO pQCD ", "l");
	}	
		
	legraa2->DrawClone("same");
	}
	

	latex->DrawLatex(Rxpos[r], 0.85,Form("R=%.1lf",R));
	
	TLine *bias = new TLine(xbiasLine[r], yminRAA, xbiasLine[r], ybiasLine);
	bias->SetLineWidth(lineBiasWidth);
	bias->SetLineStyle(2);
	bias->SetLineColor(lineBiasColor);
	if(ppBias==0 && drawBiasLine) bias->DrawClone("same");
	TLatex *latexbias = new TLatex();
	latexbias->SetNDC();
	latexbias->SetTextSize(0.04);
	latexbias->SetTextColor(lineBiasColor);
	if(ppBias==0 && drawBiasLine) latexbias->DrawLatex(xbiasDesc[r], ybiasDesc,"--> ~ UNBIASED");
	
	}//R loop
	cfig1->SaveAs(Form("%s/RAA_%s_pTl%.0lf_bin%i%s%s.%s",outDir.Data(),C_system_name[system].Data(),pTlead,binning,sufPPbias[ppBias].Data(),comp_suff[compareRAA].Data(),ext.Data()));
  }
  if(pTlead>6) return;
  if(!showRCP) return;
  
	//********************
	//RCP
	//********************

	int nR_RCP=2;
	TCanvas *cfig2=new TCanvas("cfig2","fig2",10,10,nR_RCP*860,700);
	cfig2->Divide(nR_RCP,1,0,0);
	
	for(int r=0; r<nR_RCP; r++)
	{

	float R=C_Rs[r];
	
	bool plotAliceRCP=(R<0.35) ? true : false; //we have Alice results only for R=0.2 and R=0.3
	float xmaxRCP_use=xmaxRCP;
	if(!plotAliceRCP) xmaxRCP_use=xmaxRCP_short;
	//Jan
	//fp[r]=new TFile(Form("AuAu/bin%i/systematics_R%.1lf_pT%0.lf.root",R,binning,pTlead),"read");
	TGraphAsymmErrors *gjanRCP_unf=(TGraphAsymmErrors*) fp[r]->Get("RCP_unfold_BSL");
	TGraphAsymmErrors *gjanRCP_corr=(TGraphAsymmErrors*) fp[r]->Get("RCP_corr_BSL");
	gjanRCP[r]=(TGraphAsymmErrors*) fp[r]->Get("RCP_BSL");

	ShrinkGraph(gjanRCP_unf,ptminRCP_STAR,ptmaxRCP_STAR);
	ShrinkGraph(gjanRCP_corr,ptminRCP_STAR,ptmaxRCP_STAR);
	ShrinkGraph(gjanRCP[r],ptminRCP_STAR,ptmaxRCP_STAR);
	gjanRCP_sys[r]=gcombine_syserr(gjanRCP_unf,gjanRCP_corr,1,0);

	//TGraphErrors* galiceRCP[r]=NULL;
	//TGraphAsymmErrors* galiceRCP[r]=NULL;

  //Draw plots
  //if(R>0.35) continue;
	cfig2->cd(r+1);
  gPad->SetMargin(0.15,0.03,0.20,0.03);
  if(r==0)gPad->SetRightMargin(0.0);
  if(r==1)
  {
	  gPad->SetLeftMargin(0.0);
	  //xminRCP=xminRCP+0.05;
  }
  gPad->SetTicks(1);
  if(showAliceRCP && plotLogX_RCP) gPad->SetLogx();
  if(plotLogY_RCP)gPad->SetLogy();
  htmp2[r]=new TH1D(Form("htmp2_%i",r),"",100,xminRCP,xmaxRCP_use);
  TString xtitle="";
  if(r==nR_RCP-1) xtitle="p_{T, jet}^{ch}, p_{T}^{ch}  (GeV/#it{c})  ";
  htmp2[r]->SetXTitle(xtitle);htmp2[r]->SetYTitle("R_{CP}");
  htmp2[r]->SetTitleOffset(1.05,"x");htmp2[r]->SetTitleOffset(0.95,"y");
  htmp2[r]->SetTitleSize(0.075,"x");htmp2[r]->SetTitleSize(0.080,"y");
  htmp2[r]->SetLabelSize(0.055,"x");htmp2[r]->SetLabelSize(0.055,"y");
  //htmp2->SetNdivisions(505,"x");htmp2->SetNdivisions(505,"y");
  htmp2[r]->SetMinimum(yminRCP);htmp2[r]->SetMaximum(ymaxRCP);
  htmp2[r]->DrawCopy();
  
	TLine *one = new TLine(xminRCP, 1, xmaxRCP,1);
	one->SetLineWidth(lineOneWidth);
	one->SetLineStyle(lineOneStyle);
	one->SetLineColor(lineOneColor);
	one->DrawClone("same");
  
    //ATLAS ch hadrons RCP
  gatlashRCP_sys->SetLineWidth(1);
  gatlashRCP_sys->SetLineColor(catlas);
  gatlashRCP_sys->SetFillStyle(0);
  gatlashRCP->SetLineWidth(1);
  gatlashRCP->SetLineColor(catlas);
  gatlashRCP->SetMarkerStyle(markatlasH);
  gatlashRCP->SetMarkerSize(msizeatlasH);
  gatlashRCP->SetMarkerColor(catlas);
  if(showAtlasRCP) gatlashRCP->Draw("P");
  if(showAtlasRCP) gatlashRCP_sys->Draw("2");
  
  //STAR ch. hadrons
  gSTARhRCP->SetLineWidth(1);
  gSTARhRCP->SetLineColor(cstarH);
  gSTARhRCP->SetMarkerStyle(markstarH);
  gSTARhRCP->SetMarkerSize(msizestarH);
  gSTARhRCP->SetMarkerColor(cstarH);
  if(showSTARhRCP) gSTARhRCP->Draw("P");
  
  //STAR recoil jets
  gIcp[r]->SetLineColor(calex);
  gIcp[r]->SetLineWidth(2);
  gIcp[r]->SetMarkerColor(calex);
  gIcp[r]->SetLineStyle(4);
  if(showAlex)gIcp[r]->Draw("l");
  gIcp_sys[r]->SetLineColor(calex);
  gIcp_sys[r]->SetLineWidth(2);
  gIcp_sys[r]->SetMarkerColor(calex);
  gIcp_sys[r]->SetLineStyle(3);
  if(showAlex)gIcp_sys[r]->Draw("l");
  
  //Alice inclusive jets
  galiceRCP_sys[r]->SetLineWidth(1);
  galiceRCP_sys[r]->SetLineColor(calice);
  galiceRCP_sys[r]->SetFillStyle(0);
  galiceRCP[r]->SetLineWidth(1);
  galiceRCP[r]->SetLineColor(calice);
  galiceRCP[r]->SetMarkerStyle(markalice);
  galiceRCP[r]->SetMarkerSize(msizealice);
  galiceRCP[r]->SetMarkerColor(calice);
  if(plotAliceRCP && showAliceRCP) galiceRCP[r]->Draw("P");
  if(plotAliceRCP && showAliceRCP) galiceRCP_sys[r]->Draw("2");
  
  //STAR inclusive jets  
    gjanRCP_sys[r]->SetLineWidth(1);
  gjanRCP_sys[r]->SetLineColor(cstar);
	gjanRCP_sys[r]->SetFillStyle(0);
  gjanRCP[r]->SetLineWidth(1);
  gjanRCP[r]->SetLineColor(cstar);
  gjanRCP[r]->SetMarkerStyle(markstar);
  gjanRCP[r]->SetMarkerSize(msizestar);
  gjanRCP[r]->SetMarkerColor(cstar);
	gjanRCP[r]->Draw("P");
  gjanRCP_sys[r]->Draw("2");
  
	if(r==0){
	infoXmin=0.7;
	infoXmax=0.9;
	infoYmin=(plotLogY_RCP)?0.25:0.75;
	infoYmax=(plotLogY_RCP)?0.4:0.9;
	TLegend *model_info = new TLegend(infoXmin,infoYmin,infoXmax,infoYmax);
	model_info->SetTextSize(0.05);
	model_info->SetFillStyle(0);
	model_info->SetBorderSize(0);
	model_info->SetMargin(0.05);
	//model_info->SetHeader("STAR Au+Au #sqrt{s_{NN}}=200 GeV");
	//model_info->AddEntry("", Form("Run11, %s","MB"), "");
	//model_info->AddEntry("", "Charged jets", "");
	//model_info->AddEntry("", "0-10% / 60-80% central collisions", "");
	model_info->AddEntry("",Form("anti-k_{T}, R=%.1lf",R), "");
	model_info->AddEntry("", Form("p_{T, lead}^{min} = %.0lf GeV/#it{c}",pTlead), "");
	model_info->Draw("same");
	
	latexL->DrawLatex(0.55, 0.85,label);
	
	//legend
	posx1=0.18;
	posx2=0.55;
	posy1=(plotLogY_RCP)?0.25:0.65;
	posy2=(plotLogY_RCP)?0.5:0.9;
	if(showAlex && !plotLogY_RCP) posy1=0.55;
	TLegend *legrcp = new TLegend(posx1,posy1,posx2,posy2);
	legrcp->SetTextSize(0.04);
	legrcp->SetFillStyle(0);
	legrcp->SetBorderSize(0);
	//legrcp->AddEntry(gjanRCP[r], "STAR Au+Au","lp");
	legrcp->AddEntry(gjanRCP[r], "Au+Au #sqrt{s_{NN}}=200 GeV","");
		//if(showAlex)legrcp->AddEntry(gjanRCP[r], " inclusive jets","");
	legrcp->AddEntry(gjanRCP[r], "  ch. jets 0-10% / 60-80%","lp");
		//legrcp->AddEntry(gjanRCP[r], "  0-10% / 60-80%","");
	if(showAlex)
	{
		//legrcp->AddEntry(gIcp[r], "STAR Au+Au","l");
		legrcp->AddEntry(gIcp[r], "  ch. recoil jets","l");
		legrcp->AddEntry(gIcp[r], "  0-10% / 60-80%","");
	}
	if(showSTARhRCP)
	{
		legrcp->AddEntry(gSTARhRCP,"  ch. hadrons 0-5% / 60-80%", "lp");
		//legrcp->AddEntry(gSTARhRCP,"  0-5% / 60-80%", "");
	}
	legrcp->Draw("same");

	}

	if(r==1)
	{
	infoXmin=0.75;
	infoXmax=0.9;
	infoYmin=(plotLogY_RCP)?0.25:0.8;
	infoYmax=(plotLogY_RCP)?0.4:0.9;
	TLegend *model_info2 = new TLegend(infoXmin,infoYmin,infoXmax,infoYmax);
	model_info2->SetTextSize(0.05);
	model_info2->SetFillStyle(0);
	model_info2->SetBorderSize(0);
	model_info2->SetMargin(0.05);
	model_info2->AddEntry("",Form("R=%.1lf",R), "");
	model_info2->Draw("same");
	
	posx1=0.05;
	posx2=0.5;
	posy1=(plotLogY_RCP)?0.25:0.65;
	posy2=(plotLogY_RCP)?0.5:0.9;
	TLegend *legrcp2 = new TLegend(posx1,posy1,posx2,posy2);
	legrcp2->SetTextSize(0.04);
	legrcp2->SetFillStyle(0);
	legrcp2->SetBorderSize(0);
		
	if(plotAliceRCP && showAliceRCP)	
	{
		legrcp2->AddEntry(galiceRCP[r], "Pb+Pb #sqrt{s_{NN}}=2.76 TeV","");
		legrcp2->AddEntry(galiceRCP[r], "  ch. jets 0-10% / 50-80%","lp");
		//legrcp2->AddEntry(galiceRCP[r], "  0-10% / 50-80% ","");
		legrcp2->AddEntry(gatlashRCP, "  ch. hadrons 0-5% / 60-80%","lp");
		//legrcp2->AddEntry(galiceRCP[r], "  0-10% / 50-80% ","");
	}
	legrcp2->Draw("same");
	}

	//latex->DrawLatex(Rxpos2[r], 0.75,Form("R=%.1lf",R));

	}//R loop
	
	cfig2->SaveAs(Form("%s/RCP_pTl%.0lf_bin%i%s.%s",outDir.Data(),pTlead,binning,sufPPbias[ppBias].Data(),ext.Data()));
  
}

//===============================================================
void ShrinkGraph(TGraphAsymmErrors* graph, double min, double max)
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
void ShrinkGraph(TGraphErrors* graph, double min, double max)
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

void plot_g2(TGraphErrors *g1,TGraphErrors *g2,int ic=2,int is=1,int iw=1,char *txt="",double xtxt=0.5,double ytxt=0.5,double tsiz=0.04){
  int n=g1->GetN()+g2->GetN()+1;
  //const int N=n;
  const int N=200;
  double xx[N],yy[N];
  double x,y;
  for(int i=0;i<g1->GetN();++i){
    g1->GetPoint(i,x,y);
    xx[i]=x;yy[i]=y;
  }
  for(int i=0;i<g2->GetN();++i){
    g2->GetPoint(i,x,y);
    xx[n-2-i]=x;yy[n-2-i]=y;
  }
  xx[n-1]=xx[0];yy[n-1]=yy[0];
  TGraph *g=new TGraph(n,xx,yy);
  g->SetLineColor(ic);g->SetLineStyle(is);g->SetLineWidth(iw);g->Draw("l");
  keyLine(xtxt,ytxt,txt,ic,is,tsiz,iw,1);
  return;
  /*
  g1->SetLineColor(ic);g1->SetLineStyle(is);g1->SetLineWidth(iw);g1->Draw("l");
  g2->SetLineColor(ic);g2->SetLineStyle(is);g2->SetLineWidth(iw);g2->Draw("l");
  TLine *line=new TLine();
  line->SetLineColor(ic);line->SetLineStyle(is);line->SetLineWidth(iw);
  double x1,y1,x2,y2;
  g1->GetPoint(0,x1,y1);g2->GetPoint(0,x2,y2);line->DrawLine(x1,y1,x2,y2);
  g1->GetPoint(g1->GetN()-1,x1,y1);g2->GetPoint(g2->GetN()-1,x2,y2);line->DrawLine(x1,y1,x2,y2);
  keyLine(xtxt,ytxt,txt,ic,is,tsiz,iw,1);*/
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
TGraphAsymmErrors* gcombine_syserr(TGraphAsymmErrors *gr1,TGraphAsymmErrors *gr2,bool setErrX=1,bool RAA=1/*RAA or RCP*/)
{
	int npoints=gr1->GetN();
	const int n=30;
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
	
		double norm=(RAA)? C_RAA_NormErr[0] : C_RCP_NormErr; //relative error of the normalization

		gx[i]=x1;
		gx_l[i]=x1_down;
		gx_h[i]=x1_up;
		gy[i]=y1;
		gy_l[i]=TMath::Sqrt(y1_down*y1_down+y2_down*y2_down+y1*norm*y1*norm);
		gy_h[i]=TMath::Sqrt(y1_up*y1_up+y2_up*y2_up+y1*norm*y1*norm);
		
	}
	TGraphAsymmErrors* graph=new TGraphAsymmErrors(npoints, gx, gy, gx_l, gx_h, gy_l, gy_h);
	return graph;
		
}
//===============================================================

TGraph *graph_band(int n,double *x,double *y_high,double *y_low,double xmax=10){
  for(int i=0;i<n;++i)
  //const int nn=2*n+1;double xx[nn],yy[nn];
	  const int nn=50;double xx[nn],yy[nn];
  int N=0;
  for(int i=0;i<n;++i){
    if(x[i]>xmax)continue;
    xx[N]=x[i];yy[N]=y_high[i];
    N++;
  }
  for(int i=n;i<2*n;++i){
    if(x[2*n-1-i]>xmax)continue;
    xx[N]=x[2*n-1-i];yy[N]=y_low[2*n-1-i];
    N++;
  }
  xx[N]=xx[0];yy[N]=yy[0];
  TGraph *g=new TGraphErrors(N+1,xx,yy);
  return(g);
}
