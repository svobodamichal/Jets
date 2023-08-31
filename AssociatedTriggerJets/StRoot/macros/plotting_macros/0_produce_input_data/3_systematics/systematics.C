#include "../../utils/common.h" //several constants
#include "../../utils/qa_cuts.h" //cuts on quality of the unfolded solution
#include "../../utils/utils.C" //contains sveral useful functions
#include "../../utils/iter_data.h" //optimal iterations
#include "../../utils/iter_toymodel.h" //optimal iterations
//#include "utils/fit_functions.C" //declarations of fitting functions

typedef unsigned int uint;
typedef TGraphAsymmErrors Graph;

//enum tests { chi2, camb, KS };
//enum unfoldings { Bayes, SVD };
//enum systems { cent, peri, pp };
enum ratios { closure, RAA, RCP, plpl, RR , sys };
enum syserrors { unfolding, correlated, normalization };
enum ppSets { PYTHIA, STAR };

	
//TString C_unfoldings_name[]={"Bayes","SVD"};
TString system_desc[]={"0-10% AuAu","60-80% AuAu", "pp"};
TString system_desc_short[]={"cent","peri", "pp"};

void systematics(syserrors syserr_type=correlated, ratios ratio_type=RAA, C_systems system=cent, float R=0.3, float pTlead=5.0, float pTlead_denom=2.0,
					  TString dir_suffix1="omicron",short bining=1, TString sysSource="normal", ppSets ppBase=PYTHIA, bool ppBaseWpTlead=0, bool doToy=0,bool saveSysFile=0,bool printTextFile=0)
{
	//syserr_type: type of systematic uncertainties to calculate:
	//unfolding - unfolding error
	//correlated - correlated errors (efficiency uncertainty,...)
	//normalization - normalization error (C_TAA uncertainty,...)
	TString syserr_name[]={"unfold","corr","norm"};
    
	//ratio_type: ratio to calculate
	//0 - do not calculate any ratio
	//RAA - RAA
	//RCP - RCP
	//plpl - pTlead1/pTlead2
	//RR - R1/R2
    TString ratio_name[]={"Closure","RAA","RCP","pTlpTl","RR"};
    
	
	//system:
	//cent - central AuAu
	//peri - peripheral AuAu
	//pp - pp collisions
   
	//***********************
	//initial setting
	//***********************
	bool printInfo=0; //verbose mode
	float pTlead_denom_min[]={5.0,2.0,4.0};  //minimal pTlead cut value for the donominator in the ratio of spectra with two different pTlead cuts
	bool ppHT=1; //for pp, use combined MB+HT pp data?: 0=MB only | 1=MB+HT
	bool showBBB=0; //show also bin-by-bin corrected solution (for pp data only)
	//bool ppBaseWpTlead=1; //0: unbiased pp reference | 1: apply pTlead cut also on pp reference
	bool useMPV=0; //use most probable value instead of mean for ratios obtained through MC propagation of two gaussians
	C_tests test=camb; //statistical test for finding optimal solutions: 0 chi2, 1 camberra, 2 Kolmogorov-Smirnov
	//ppSets ppBase=PYTHIA; //pp baseline for RAA: [STAR|PYTHIA]
	short unf=2; //which unfolding to use: 0: Bayes, 1: SVD, 2: both
	bool separate_solutions=1; //true: calculate unfolding uncertainty for RAA and RCP as a spread around average between various RAAs/RCPs for each prior
	TString dir_suffix2=dir_suffix1; //use the same evolution name for denominator as for numerator 
	bool constant_iter=0; //0: use optimal iteration, 1: use fixed iteration defiter
	short defiter[C_nunfoldings]={3,4}; //default iteration [Bayes,SVD] which will be used in case there is no good iteration at all or in case we want to use constant iteration instead of default
	unsigned short err_norm_spec=1; //normalization of systematic errors for spectra: 0: 1, 1: sqrt(n), 2: sqrt(n*(n-1))
	float RAA_toy[/*centrality*/]={0.5,0.8,0.5}; //RAA value for toymodel
	//priors
	int priors[]={5,2,4,6,7,8,9,10,11,12,13,14,15};
	int npriors=13;
	const int nshowsys=3; //for how many pT bins we want to print out the size of the correlated systematic errors
	float pTshowsys[nshowsys]={15,18.0,22.5}; //pT for which we want to print out the size of the (correlated) systematic errors
		//if(system==peri)=pTshowsys[1]=18.0;
	
	
	//array sizes
	const short nunfolds=(unf==2) ? C_nunfoldings : 1; //number of unfolding types
	const int nsolutions_max=npriors*2*nunfolds; //number of solutions to compare, N_priors*(i, i+1 iteration)*N_unfoldings
	const int nratios_max=20; //maximal number of averages to compare - needed for syserr_type==correlated and syserr_type==normalization
	
	//declare variables
	float R1, R2, pTl1, pTl2, yline;
	unsigned short nrats,err_norm;
	bool useaverages4sys; 
	C_systems system1, system2;
	TString str, y_title, sdesc1, sdesc2, rdesc;

	//set default values, DO NOT CHANGE! (desired values for each specific case are set in 'switch (ratio_type)' and 'switch (syserr_type)' )
	R1=R;
	R2=R;
	pTl1=pTlead;
	pTl2=pTlead;
	system1=system; //system for the numerator (central, peripheral, pp)
	system2=system; //system for the denominator (central, peripheral, pp)
	nrats=1; //use (nrats-1) average values to calculate the systematic error (nrats=1: calculate it as spread of different solutions, do not use average values)
	useaverages4sys=true; //0: calculate the systematic error from the spread of ratios of different solutions 
							//1, (nrats=1): use errors of average values of numerator and denominator to calculate error of the ratio
							//1, (nrats>1): use (nrats-1) average values to calculate the systematic error
	err_norm=0; //normalization of systematic errors for ratios: 0: 1, 1: sqrt(n), 2: sqrt(n*(n-1))

	//declare arrays
	TString insuf1[nratios_max]; //suffix for path
	TString insuf2[nratios_max]; 
	C_unfoldings unf1[nsolutions_max][nratios_max]; //unfolding type
	C_unfoldings unf2[nsolutions_max][nratios_max];
	uint priorNo1[nsolutions_max][nratios_max]; //prior type
	uint priorNo2[nsolutions_max][nratios_max]; 
	short iter1[nsolutions_max][nratios_max]; //number of iterations
	short iter2[nsolutions_max][nratios_max];

	//graphics
	Color_t colorList[]={13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,28,13,14,15,16,17,18,19,20,13,14,
		15,16,17,18,19,20,21,22,23,24,25,26,29,28,13,14,15,16,17,18,19,20,13,14,15,16,17,18,19,20,21,22,23,24,25,26}; //for spectra
	Color_t colorList2[]={kGray+3,kRed,kRed,kGreen+2,kGreen+2,kBlue,kBlue,kBlue,kYellow+2,kOrange-9,kMagenta-5,kCyan+3,kOrange}; //for ratio histograms
	Color_t colorList3[]={kMagenta,kViolet-6,kCyan+2}; //for total sys errors
	float marker_size=1.0;
	float line_width=2;
	Int_t markerS[]={22,23,29,20,21}; //marker style for unfolding
	Int_t markerSsys[]={1,22,26,21,25,31,29,30,20,33,34,5,2}; //marker style for systematic variations

	//canvas size
	int can_x=1200; //1600
	int can_y=680; //900
	//histogram ranges
	float spectraXmin=-20;
	float spectraXmax=50;
	float spectraYmin=1E-9;
	float spectraYmax=1E-2;
	float spec_pTmax=40.0; //maximal pT value to be saved in tgraphs
	float spec_pTmin; //calculated later
	float ratio_pTmax=40.0; //maximal pT value to be saved in tgraphs
	float ratio_pTmin; //calculated later
	float ratioXmin=0;
	float ratioXmax=(system==cent) ? 32 : 27;
	float ratioYmin=0.01;
	float ratioYmax=1.2;
		
	TH1I *frame = new TH1I("frame", "", 1000, -100, +100);
	//TH1I *frame2 = new TH1I("frame2", "", 1000, -100, +100);
	
	for(int rat=0; rat<nratios_max; rat++){
	//triple loop over all solutions
	int isol=0;
	for(int uf=0; uf<nunfolds; uf++){
	for(int itr=0; itr<2; itr++){
	for(int pr=0; pr<npriors; pr++){
		priorNo1[isol][rat]=priors[pr];
		priorNo2[isol][rat]=priors[pr];
		unf1[isol][rat]=uf;
		unf2[isol][rat]=uf;
		iter1[isol][rat]=itr; //here we only set the number of iteration to 0 or 1 , since we don't know the optimal value yet. Once we know it, it will be added to the current value (0,1 => N, N+1)
		iter2[isol][rat]=itr; 
		//cout<<isol<<" prior: "<<C_prior_type[priorNo1[isol]].Data()<<" iteration: "<<iter1[isol]<<" unfolding: "<<unf1[isol]<<endl;
		isol++;
	}}} //solutions
	} //ratios
	
	//***************************************************
	//which ratio of the unfolded spectra do we want to calculate
	//***************************************************
	
	switch (ratio_type)
	{
		case RAA : //RAA
			useaverages4sys=(separate_solutions) ? false : true; //this affects only unfolding errors, correlated errors have always useaverages4sys=true
			system1=system; 
			system2=pp;

			y_title="R_{AA}";
			rdesc=Form("p_{T}^{leading}>%.1lf",pTl1);
			yline=1.0;
			break;
			
		case RCP : //RCP
			useaverages4sys=(separate_solutions) ? false : true; //this affects only unfolding errors, correlated errors have always useaverages4sys=true
			system1=cent; 
			system2=peri;

			y_title="R_{CP}";
			rdesc=Form("p_{T}^{leading}>%.1lf",pTl1);
			yline=0.35;
			break;
			
		case plpl : //ratios of different pTleading cuts
			useaverages4sys=false; //this affects only unfolding errors, correlated errors have always useaverages4sys=true
			pTl2=(pTlead_denom<pTlead_denom_min[system]) ? pTlead_denom_min[system] : pTlead_denom; //pTlead value of the denominator
			y_title=Form("p_{T}^{leading}>%.1lf/p_{T}^{leading}>%.1lf",pTl1,pTl2);
			rdesc=Form("R=%.1lf",R1);
			yline=1.0;
			ratio_name[3]=Form("pTlpTl%.0lf",pTl2);
			break;

		case RR :// ratio R1/R2
			useaverages4sys=false; //this affects only unfolding errors, correlated errors have always useaverages4sys=true
			R1=0.2;
			R2=R;
			if(R<0.21 && R>0.19) 
			{
				R1=0.3;
				R2=0.4;
			}
			y_title=Form("ratio R=%.1lf/R=%.1lf",R2,R1);
			rdesc=Form("p_{T}^{leading}>%.1lf",pTl1);
			yline=1.0;
			break;
			
		case sys : //compare various systematic effects
			useaverages4sys=false; //this affects only unfolding errors, correlated errors have always useaverages4sys=true
			//syserr_type=correlated;

			y_title="variation/nominal";
			rdesc=Form("%s, R=%.1lf, p_{T}^{leading}>%.1lf", system_desc[system1].Data(),R1,pTl1);
			yline=1.0;
			ratioYmin=0.5;
			ratioYmax=1.6;
			marker_size=2;
			break;
			
		default : //no ratio at all
			useaverages4sys=false;
			yline=1.0;

	}
	
	float pTlmax=(pTl2>pTl1) ? pTl2 : pTl1; //maximum from (pTl1,pTl2)
	ratio_pTmin=pTlmax;
    spec_pTmin=pTl1;
	
	//**********************************************
	//systematic error which we want to calculate
	//**********************************************

	TString leg_name[nratios_max];
	TString sys_source[nratios_max];
				sys_source[0]=sysSource; //default systematic variation (usualy "normal") when we do not look at correlated uncertainties
				sys_source[1]=sysSource; //default systematic variation (usualy "normal") when we do not look at correlated uncertainties
				sys_source[2]=sysSource; //default systematic variation (usualy "normal") when we do not look at correlated uncertainties
	switch (syserr_type)
	{
		case unfolding : //unfolding error
			nrats=1; 
			err_norm=1; //normalization of systematic errors: 0: 1, 1: sqrt(n), 2: sqrt(n*(n-1))
			break;
		case correlated : //correlated errors
			useaverages4sys=true;
			if(system1==pp)nrats=7; //cannot be higher than nratios_max!!!
			else nrats=11; //13
			err_norm=0; //normalization of systematic errors: 0: 1, 1: sqrt(n), 2: sqrt(n*(n-1))
			
			sys_source[0]="normal";
			sys_source[1]="p5";
			sys_source[2]="m5";
			sys_source[3]="g";
			sys_source[4]="u";
			
			if(system1==pp)
			{
			sys_source[5]="tofp5";
			sys_source[6]="tofm5";
			}
			else
			{
			sys_source[5]="nrem-1";
			sys_source[6]="RRho02";
			sys_source[7]="RRho04";
			sys_source[8]="pp";
			//sys_source[9]="nfit14";
			//sys_source[10]="nfit18";
			sys_source[9]="v2";
			sys_source[10]="pythiaU";
			}
			
			leg_name[0]="nominal";
			leg_name[1]="#epsilon(TPC) +5";
			leg_name[2]="#epsilon(TPC) -5";
			leg_name[3]="g frag.";
			leg_name[4]="u frag.";
			if(system1==pp)
			{
			leg_name[5]="#epsilon(TOFvBEMC) +5";
			leg_name[6]="#epsilon(TOFvBEMC) -5";
			}
			else
			{
			leg_name[5]="#rho: nremove=x-1";
			leg_name[6]="#rho: R=0.2";
			leg_name[7]="#rho: R=0.4";
			leg_name[8]="pp had. rat.";
			//leg_name[9]="nfit>=14";
			//leg_name[10]="nfit>=18";
			leg_name[9]="v2 corr";
			leg_name[10]="pythia dpT";
			
			}
			
			break;
		case normalization : //normalization errors
			if(ratio_type!=RAA) 
			{
				cout<<"WARNING: Normalization errors can be calculated only for RAA!"<<endl;
				cout<<"TERMINATING..."<<endl;
				return;
			}
			useaverages4sys=true;
			nrats=3; //cannot be higher than nratios_max!!!
			err_norm=0; //normalization of systematic errors: 0: 1, 1: sqrt(n), 2: sqrt(n*(n-1))
			break;
	}
	
		
	//fill arrays
	for(int rat=0; rat<nratios_max; rat++)
	{
		insuf1[rat]=Form("_%s",dir_suffix1.Data());
		insuf2[rat]=Form("_%s",dir_suffix2.Data());
	
		if(syserr_type!=correlated)
		{
			insuf1[rat]+=Form("_%s",sysSource.Data());
			insuf2[rat]+=Form("_%s",sysSource.Data());
		}
		else
		{
			insuf1[rat]+=Form("_%s",sys_source[rat].Data());
			insuf2[rat]+=Form("_%s",sys_source[rat].Data());
		}
	}
	
		
	//declare additional variables
	const int nratios=nrats; //size of the ratio histogram arrays
	//calculate index of R12 and pTl12 in Rs[] and pTls[] arrays (see "utils/common.h" and "utils/utils.C")
	int R1_idx=calculate_index(R1,0);
	int R2_idx=calculate_index(R2,0);
	int pTl1_idx=calculate_index(pTl1,1);
	int pTl2_idx=calculate_index(pTl2,1);
	
	float pTl_pythia=(ppBaseWpTlead) ? pTlead : 0; //pTleading cut value for pp baseline 
	
	TString sdesc1=Form("%s, R=%.1lf, p_{T}^{lead}>%.0lf",system_desc[system1].Data(),R1,pTl1);
	TString sdesc2=Form("%s, R=%.1lf, p_{T}^{lead}>%.0lf",system_desc[system2].Data(),R2,pTl2);
	
	
	if(printInfo) cout<<"R1_idx, R2_idx:"<<R1_idx<<" "<<R2_idx<<" pTl1_idx, pTl2_idx: "<<pTl1_idx<<" "<<pTl2_idx<<endl;
	
	//fill optimal iteration value array
	for(int rat=0; rat<nratios; rat++){
	for(int isol=0; isol<nsolutions_max; isol++)
	{
		int *optiter;
		if(doToy) optiter=diter_toy_normal;
		else
		{
			if(syserr_type!=correlated)
			{
				optiter=(system==peri) ? diter_data_peripheral_normal : diter_data_normal;
			}
			else
			{
				if(rat==0) optiter=(system==peri) ? diter_data_peripheral_normal : diter_data_normal;
				else if (rat==1) optiter=(system==peri) ? diter_data_peripheral_p5 : diter_data_p5;
				else if (rat==2) optiter=(system==peri) ? diter_data_peripheral_m5 : diter_data_m5;
				else if (rat==3) optiter=(system==peri) ? diter_data_peripheral_g : diter_data_g;
				else if (rat==4) optiter=(system==peri) ? diter_data_peripheral_u : diter_data_u;
				else if (rat==5) optiter=(system==peri) ? diter_data_peripheral_nrem1 : diter_data_nrem1;
				else if (rat==6) optiter=(system==peri) ? diter_data_peripheral_RRho02 : diter_data_RRho02;
				else if (rat==7) optiter=(system==peri) ? diter_data_peripheral_RRho04 : diter_data_RRho04;
				else if (rat==8) optiter=(system==peri) ? diter_data_peripheral_pp : diter_data_pp;
				else if (rat==9) optiter=(system==peri) ? diter_data_peripheral_trcuts2 : diter_data_trcuts2;
				else if (rat==10) optiter=(system==peri) ? diter_data_peripheral_v2 : diter_data_v2;
				else if (rat==11) optiter=(system==peri) ? diter_data_peripheral_pythia : diter_data_pythia;
			}
		}

		
		short o1=optiter[priorNo1[isol][rat]][R1_idx][pTl1_idx][unf1[isol][rat]][system1][bining][test];
		//cout<<"diter_data_normal["<<priorNo1[isol][rat]<<"]["<<R1_idx<<"]["<<pTl1_idx<<"]["<<unf1[isol][rat]<<"]["<<system1<<"]["<<bining<<"]["<<test<<"]"<<endl;
		short o2=optiter[priorNo2[isol][rat]][R2_idx][pTl2_idx][unf2[isol][rat]][system2][bining][test];
		if(o1>0) iter1[isol][rat]+=o1;
		else iter1[isol][rat]=defiter[unf1[isol][rat]]; //in case we don't have the optimal iteration, set it to default value
		if(o2>0)iter2[isol][rat]+=o2;
		else iter2[isol][rat]=defiter[unf2[isol][rat]]; //in case we don't have the optimal iteration, set it to default value
		
		if(constant_iter>0)
		{
			iter1[isol][rat]=defiter[unf1[isol][rat]];
			iter2[isol][rat]=defiter[unf2[isol][rat]];
		}
		
		if(printInfo) 
			cout<<"iteration1,2:"<<iter1[isol][rat]<<", "<<iter2[isol][rat]<<endl;
	}}


		
	//************************************
	//Load histograms from input files
	//************************************
	TString sysname1,sysname2,RAAsuf1,RAAsuf2;
	TString unfType1="BGD";
	TString unfType2="BGD";
	if(doToy)
	{ 
		sysname1="toymodel";
		sysname2="toymodel";
		RAAsuf1=Form("_RAA%.1lf",RAA_toy[system1]);
		RAAsuf2=Form("_RAA%.1lf",RAA_toy[system2]);
	}
	else
	{
		sysname1="MB";
		sysname2="MB";
		RAAsuf1="";
		RAAsuf2="";
	}
	if(system1==peri) sysname1+="_peripheral";
	else if(system1==cent) sysname1+="_central";
	else if(system1==pp){
		sysname1+=(ppHT) ? "HT_pp" : "_pp";
		unfType1="dete";
	}
	if(system2==peri) sysname2+="_peripheral";
	else if(system2==cent) sysname2+="_central";
	else if(system2==pp){
		sysname2+=(ppHT) ? "HT_pp" : "_pp";
		unfType2="dete";
	}
	
	TString wrkdir1[nsolutions_max][nratios];
	TString wrkdir2[nsolutions_max][nratios];

	int nevents1[nratios];
	int nevents2[nratios];
	
	for(int rat=0; rat<nratios; rat++){
	for(int isol=0; isol<nsolutions_max; isol++)
	{
		wrkdir1[isol][rat] = Form("../../../plotting_out/root/%s/%s/Unfolded_R%.1lf_%s_VARbins_bining%i_%s%s%s",sysname1.Data(),dir_suffix1.Data(),R1,C_unfoldings_name[unf1[isol][rat]].Data(),bining,unfType1.Data(),RAAsuf1.Data(),insuf1[rat].Data());
		wrkdir2[isol][rat] = Form("../../../plotting_out/root/%s/%s/Unfolded_R%.1lf_%s_VARbins_bining%i_%s%s%s",sysname2.Data(),dir_suffix2.Data(),R2,C_unfoldings_name[unf2[isol][rat]].Data(),bining,unfType2.Data(),RAAsuf2.Data(),insuf2[rat].Data());
		
		if(ppBase==PYTHIA && ratio_type==RAA)wrkdir2[isol][rat]=wrkdir1[isol][rat];
	}//unf. solution loop
	
	//get total number of events
	str = Form("%s/../histos_inclusivejet_%s.root", wrkdir1[0][0].Data(),sys_source[rat].Data());
	if(doToy) str = Form("%s/../root_RAA%.1lf_%s/histos_jets_R%.1lf_pTcut0.2.root", wrkdir1[0][0].Data(),RAA_toy[system1],dir_suffix1.Data(),R1);
	TFile *f1 = new TFile(str.Data(), "OPEN");
	TH1D* hevents1= (TH1D*) f1->Get("hevents");
	str = Form("%s/../histos_inclusivejet_%s.root", wrkdir2[0][0].Data(),sys_source[rat].Data());
	if(doToy) str = Form("%s/../root_RAA%.1lf_%s/histos_jets_R%.1lf_pTcut0.2.root", wrkdir2[0][0].Data(),RAA_toy[system2],dir_suffix2.Data(),R2);
	TFile *f2 = new TFile(str.Data(), "OPEN");
	TH1D* hevents2= (TH1D*) f2->Get("hevents");
	
	nevents1[rat]=hevents1->GetBinContent(2);
	nevents2[rat]=hevents2->GetBinContent(2);
	if(doToy)
	{
		nevents1[rat]=hevents1->GetEntries();
		nevents2[rat]=hevents2->GetEntries();
	}
	
	cout<<sys_source[rat].Data()<<endl;
	cout<<"nevents1: "<<nevents1[rat]<<endl;
	cout<<"nevents2: "<<nevents2[rat]<<endl;
	
	f1->Close();
	f2->Close();
	}//systematic sources loop
//return;

	//show bin-by-bin correction as well
	if(system==pp && showBBB)
	{
		str = Form("%s/../histos_inclusivejet_BBBcorr.root", wrkdir1[0][0].Data());
		TFile *fbbb = new TFile(str.Data(), "OPEN");
		TH1D* hbbb=(TH1D*) fbbb->Get(Form("hpT_pTl%.0lf_R0%.0lf_corr",pTlead,R*10));
		TH1I* heventsbbb= (TH1I*) fbbb->Get("hevents");
		//int neventsbbb=heventsbbb->GetEntries();
		int neventsbbb=heventsbbb->GetBinContent(2);
	}
	
	//cout<<"nevents1: "<<nevents1<<endl;
	//cout<<"nevents2: "<<nevents2<<endl;
	
	//TOYMODEL true spectrum
	if(doToy && ratio_type==closure)
	{
		str = Form("%s/../root_RAA%.1lf_%s/jetonly_R%.1lf_pTcut0.2.root", wrkdir1[0][0].Data(),RAA_toy[system1],dir_suffix1.Data(),R1);
		TFile *ftrue=new TFile(str.Data(), "OPEN");
		TH2D* htrue2D=ftrue->Get("fhPtRecpTleading");
		int bint_first=htrue2D->GetXaxis()->FindBin(pTlead);
		int bint_last=htrue2D->GetNbinsX();
		TH1D* htrue=htrue2D->ProjectionY("htrue",bint_first,bint_last);
		TH1I* hevents_true=ftrue->Get("hevents");
		for(int rat=0; rat<nratios; rat++)
		{
			nevents2[rat]=hevents_true->GetEntries();
            cout<<"nevents2: "<<nevents2[rat]<<endl;
		}
		
	}
	
	//Normalization
	double scale_jets1[nratios];
	double scale_jets2[nratios];
	for(int rat=0; rat<nratios; rat++){
		scale_jets1[rat] = 1./(2.*(1.-R1)*2.*TMath::Pi()*nevents1[rat]);
		scale_jets2[rat] = 1./(2.*(1.-R2)*2.*TMath::Pi()*nevents2[rat]);
		if(dir_suffix1=="eta1-R03")
		{
			scale_jets1[rat]= 1./(2.*(1.-0.3)*2.*TMath::Pi()*nevents1[rat]);
			scale_jets2[rat]= 1./(2.*(1.-0.3)*2.*TMath::Pi()*nevents1[rat]);
		}	
		if(ratio_type==RAA && ppBase==PYTHIA)
			scale_jets2[rat] = 1.0; //PYTHIA spectrum is already normalized
	}//ratio (=syserr source) loop
	

	
	double scale_denum[/*nratios_max*/]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; //aditional normalization for the denominator
	if(ratio_type==RAA)
	{
		if(ppBase==PYTHIA) 
		{
			for(int rat=0; rat<nratios; rat++)
			{scale_denum[rat]=C_TAA[system1];}
			if(syserr_type==normalization)
			{
				scale_denum[1]=C_TAA[system1]+C_TAA_err[system1];
				scale_denum[2]=C_TAA[system1]-C_TAA_err[system1];
			}
		}
		else if(ppBase==STAR) 
		{
			for(int rat=0; rat<nratios; rat++)
			{scale_denum[rat]=C_NBIN[system1]/C_NBIN[system2];}
			if(syserr_type==normalization)
			{
				scale_denum[1]=C_NBIN[system1]/C_NBIN[system2]*(1+C_NBIN_rerr[system1]+C_NBIN_rerr[system2]);
				scale_denum[2]=C_NBIN[system1]/C_NBIN[system2]*(1-C_NBIN_rerr[system1]-C_NBIN_rerr[system2]);
			}
		}
	}
	else if(ratio_type==RCP)
	{
		for(int rat=0; rat<nratios; rat++)
		{scale_denum[rat]=C_NBIN[system1]/C_NBIN[system2];}
		if(syserr_type==normalization)
		{
			scale_denum[1]=C_NBIN[system1]/C_NBIN[system2]*(1+C_NBIN_rerr[system1]+C_NBIN_rerr[system2]);
			scale_denum[2]=C_NBIN[system1]/C_NBIN[system2]*(1-C_NBIN_rerr[system1]-C_NBIN_rerr[system2]);
		}
	}
	
	TFile *funfolding1[nsolutions_max][nratios];
	TFile *funfolding2[nsolutions_max][nratios];
	TFile* fqa1[nsolutions_max][nratios];
	TFile* fqa2[nsolutions_max][nratios];
    
	TH1D* hunfolded1[nsolutions_max];
	TH1D* hunfolded2[nsolutions_max];
	TH1D* hunfolded1_good12[nsolutions_max][nratios]; //array with solutions which are OK both for numerator and denominator
	TH1D* hunfolded2_good12[nsolutions_max][nratios];
	TH1D* hunfolded1_avg[nratios]; //average value
	TH1D* hunfolded2_avg[nratios]; //average value
	TH1D* hunfolded2_avg_scaled[nratios]; //average value
	TH1D* hratio_avg[nratios]; //ratio of averages
	TH1D* hratio_var[nsolutions_max][nsolutions_max]; //ratio of different solutions
	TH1D* hpythia[nratios];
    
	TF1* FitF1;// function for fitting the spectrum
	Graph* gspec_avg; //tgraph for spectrum average
	Graph* gspec1_BSF; //tgraph for spectrum average, with corrected bin center (fit)
	Graph* gspec1_BSL; //tgraph for spectrum average, with corrected bin center (linear interpolation)
	Graph* gspec2_BSL; //tgraph for denominator spectrum average, with corrected bin center (linear interpolation)
	Graph* gspec_sys; //tgraph for systematic error
	Graph* gratio_avg; //tgraph for ratio average
	Graph* gratio_sys; //tgraph for systematic error
	Graph* gratio_BSL; //tgraph for ratio average, with corrected bin center (linear interpolation)
	Graph* gratio_sys_BSL; //tgraph for systematic error on ratio, with corrected bin center (linear interpolation)
	Graph* gtoytrue; //toymodel true spectrum
	
	int nsolutions1[nratios]=0; //number of OK solutions for numerator
	int nsolutions2[nratios]=0; //number of OK solutions for denominator
	int nsolutions12[nratios]=0; //number of OK solutions for both numerator & denominator
	
	TCanvas *cspectrum1[nratios];
	TCanvas *cspectrum2[nratios];
		
	TString filen;
	TString fileqa;
    
    short firstbin_spec;
    short firstbin;

    
	for(int rat=0; rat<nratios; rat++){
		int igsol1=0;  //incremented for every good solution
		int igsol2=0;
		int igsol12=0;
	//loop over unfolded solutions with different priors / effects
	for(int isol=0;isol<nsolutions_max;isol++)
	{
		if(printInfo)cout<<"loading files"<<endl;
		filen = Form("%s/%s/unfolded_SignalSpectrum_R%.1lf_pTthresh%.1lf.root", wrkdir1[isol][rat].Data(),C_prior_type[priorNo1[isol][rat]].Data(), R1, pTl1);
		if(printInfo)cout<<filen.Data()<<endl;
		funfolding1[isol][rat] = new TFile(filen.Data(), "OPEN");
		fileqa = Form("%s/%s/unfoldingQA_R%.1lf_pTthresh%.1lf.root", wrkdir1[isol][rat].Data(),C_prior_type[priorNo1[isol][rat]].Data(), R1, pTl1);
		fqa1[isol][rat] = new TFile(fileqa.Data(), "OPEN");
		
		if(ratio_type!=RAA || ppBase!=PYTHIA) //not PYTHIA
		{
		filen = Form("%s/%s/unfolded_SignalSpectrum_R%.1lf_pTthresh%.1lf.root", wrkdir2[isol][rat].Data(),C_prior_type[priorNo2[isol][rat]].Data(), R2, pTl2);
		if(printInfo)cout<<filen.Data()<<endl;
		funfolding2[isol][rat] = new TFile(filen.Data(), "OPEN");
		fileqa = Form("%s/%s/unfoldingQA_R%.1lf_pTthresh%.1lf.root", wrkdir2[isol][rat].Data(),C_prior_type[priorNo2[isol][rat]].Data(), R2, pTl2);
		fqa2[isol][rat] = new TFile(fileqa.Data(), "OPEN");
		}
        
		bool good1=false;
		bool good2=false;
		//numerator
		if(iter1[isol][rat]>0)  //this prior has an optimal iteration and can be used
		{
			//perform QA checks
			if(printInfo) cout<<"performing QA check"<<endl;
			//load QA vectors
			TVectorD *test_change = (TVectorD*)fqa1[isol][rat]->Get(Form("%schange",C_test_name[test].Data())); //ratio of succeessive iterations
			TVectorD *test_backf_trunk=(TVectorD*)fqa1[isol][rat]->Get(Form("%sbackf_trunk",C_test_name[test].Data()));
			TVectorD *test_curvature=(TVectorD*)fqa1[isol][rat]->Get("Curvature");
			float test_ch=test_change(iter1[isol][rat]-1);
			float test_bf=test_backf_trunk(iter1[isol][rat]-1);
			float test_curv=test_curvature(iter1[isol][rat]-1);
			good1=qa_check(test_ch,test_bf,test_curv,system1,unf1[isol][rat],bining,test,R1_idx,pTl1_idx);
        
			if(good1)
			{
				str = Form("iter%d",iter1[isol][rat]-1);
				TDirectoryFile* dir1 = (TDirectoryFile*)funfolding1[isol][rat]->Get(str.Data());
				hunfolded1[igsol1] = (TH1D*) dir1->Get("hunfolded");
				//hunfolded1[igsol1]->Sumw2();
				hunfolded1[igsol1]->Scale(scale_jets1[rat],"width");
				histo_set_atributes(hunfolded1[igsol1], colorList[isol+1], colorList[isol+1], line_width, marker_size, markerS[unf1[isol][rat]]);
				igsol1++;
			}
			
			fqa1[isol][rat]->Close();
		}
	
		//denominator
		if(ratio_type==RAA && ppBase==PYTHIA) //PYTHIA
		{
			if(isol==0)
			{
				str = Form("iter%d",defiter[unf1[isol][rat]]-1);
				TDirectoryFile* dir2 = (TDirectoryFile*)funfolding1[0][rat]->Get(str.Data());
				TH1D* htemplate=(TH1D*)dir2->Get("hunfolded");
				//str=Form("../../../plotting_out/root/pp_Fuqiang/Fuqiang_fastjet2/NLO_convolute_starjet_pythia_R%.1lf_parton2chrg_new2_syst.root",R2); 
				//str="../../../plotting_out/root/pp_PYTHIA_Fuqiang/Jan_fastjet3/starjet_pythia_xsection_charged_jets.root";
				str="../../../plotting_out/root/pp_PYTHIA_Jana/jethistos.root";
				TFile* fpp = new TFile(str.Data(), "OPEN");
				//str="hNLO_conv";
				//str="hpt_chrg";
				//str=Form("hpT_R0%.0lf_pTl%.0lf",R2*10,pTl_pythia);
				str=Form("hchjet_pT_R0%.0lf_pTl%.0lf_sum_fudge",R2*10,pTl_pythia);
				hpythia[rat]=(TH1D*)fpp->Get(str.Data());
				hunfolded2[0] = rebinhisto(hpythia[rat],htemplate,Form("hpythia_%i",rat),true); //input histogram is normalized by the bin width -> set the last argument to TRUE
				hunfolded2[0]->Scale(scale_jets2[rat]); //we don't use Scale(..., "width"), because the output of the rebinhisto() function is already normalized by the binwidth
				histo_set_atributes(hunfolded2[0], colorList[isol+1], colorList[isol+1], line_width, marker_size, markerS[unf2[isol][rat]]);
				igsol2++;
				good2=true;
				
				if(system==pp && rat==0 && showBBB)
				{
					
					hbbb->Scale(scale_jets1[rat],"width");
					hbbb->Scale((float)(nevents1[rat])/neventsbbb);
					TH1D* hbbb_rebin = (TH1D*) rebinhisto(hbbb,htemplate,"hbbb_rebin",true);
					//hbbb_rebin->Scale(scale_jets1,"width");
				}
			}
			else
			{
				continue;
				/*
				hunfolded2[isol]=hunfolded2[0];
				igsol2++;
				good2=true;*/
			}
		}//PYTHIA
		else if(ratio_type==closure && doToy) //TOY CLOSURE
		{
			if(isol==0)
			{
				str = Form("iter%d",defiter[unf1[isol][rat]]-1);
				TDirectoryFile* dir2 = (TDirectoryFile*)funfolding1[0][0]->Get(str.Data());
				TH1D* htemplate=(TH1D*)dir2->Get("hunfolded");
				htrue->Scale(scale_jets2[rat],"width");
				hunfolded2[0] = rebinhisto(htrue,htemplate,Form("htrue_%i",rat),true);
				histo_set_atributes(hunfolded2[0], colorList[isol+1], colorList[isol+1], line_width, marker_size, markerS[unf2[isol][rat]]);
				igsol2++;
				good2=true;
				if(rat==0)gtoytrue=histo2graph(hunfolded2[0]);
				tgraph_set_atributes(gtoytrue, 0, spec_pTmax, spectraYmin, spectraYmax, kRed+1,kRed+1,kRed+1, 1, 2, 22, 0);
			}
			else
			{
				hunfolded2[isol]=hunfolded2[0];
				igsol2++;
				good2=true;
			}
		}//TOY CLOSURE
		
		else{
			if(iter2[isol][rat]>0){  //this prior has an optimal iteration and can be used
				//perform QA checks
				if(printInfo) cout<<"performing QA check"<<endl;
				//load QA vectors
				TVectorD *test_change = (TVectorD*)fqa2[isol][rat]->Get(Form("%schange",C_test_name[test].Data())); //ratio of succeessive iterations
				TVectorD *test_backf_trunk=(TVectorD*)fqa2[isol][rat]->Get(Form("%sbackf_trunk",C_test_name[test].Data()));
				TVectorD *test_curvature=(TVectorD*)fqa2[isol][rat]->Get("Curvature");
				float test_ch=test_change(iter2[isol][rat]-1);
				float test_bf=test_backf_trunk(iter2[isol][rat]-1);
				float test_curv=test_curvature(iter2[isol][rat]-1);
				good2=qa_check(test_ch,test_bf,test_curv,system2,unf2[isol][rat],bining,test,R2_idx,pTl2_idx);
        
				if(good2)
				{
					str = Form("iter%d",iter2[isol][rat]-1);
					TDirectoryFile* dir2 = (TDirectoryFile*)funfolding2[isol][rat]->Get(str.Data());
					hunfolded2[igsol2] = (TH1D*) dir2->Get("hunfolded");
					//hunfolded2[igsol2]->Sumw2();
					hunfolded2[igsol2]->Scale(scale_jets2[rat],"width");
					histo_set_atributes(hunfolded2[igsol2], colorList[isol+1], colorList[isol+1], line_width, marker_size, markerS[unf2[isol][rat]]);
					igsol2++;
				}
				
				fqa2[isol][rat]->Close();
			}//good iteration
		}//STAR data, not PYTHIA
		
		
		if(good1 && good2) //both numerator and denominator are OK
		{
			hunfolded1_good12[igsol12][rat]=hunfolded1[igsol1-1];
			hunfolded2_good12[igsol12][rat]=hunfolded2[igsol2-1];
			igsol12++;
		}

	}// loop over unfolded solutions
	
	
		nsolutions1[rat]=igsol1;
		nsolutions2[rat]=igsol2;
		//nsolutions12[rat]=igsol12;
		nsolutions12[rat]=igsol1*igsol2;
		cout<<"good1:"<<igsol1<<endl;
		cout<<"good2:"<<igsol2<<endl;
		cout<<"good12:"<<nsolutions12[rat]<<endl;
		
		//in case there is no good solution take at least one (not so) bad
		if(nsolutions1[rat]==0)
		{
			str = Form("iter%d",defiter[unf1[0][rat]]-1);
			cout<<"WARNING: no good solution => setting default iter to "<<defiter[unf1[0][rat]]<<endl;
			TDirectoryFile* dir1 = (TDirectoryFile*)funfolding1[0][rat]->Get(str.Data());
			hunfolded1[0] = (TH1D*) dir1->Get("hunfolded");
			hunfolded1[0]->Scale(scale_jets1[rat],"width");
			histo_set_atributes(hunfolded1[0], colorList[0+1], colorList[0+1], line_width, marker_size, markerS[unf1[0][rat]]);
			nsolutions1[rat]=1;
		}
		
		if(nsolutions2[rat]==0)
		{
			str = Form("iter%d",defiter[unf2[0][rat]]-1);
			cout<<"WARNING: no good solution => setting default iter to "<<defiter[unf2[0][rat]]<<endl;
			TDirectoryFile* dir2 = (TDirectoryFile*)funfolding2[0][rat]->Get(str.Data());
			hunfolded2[0] = (TH1D*) dir2->Get("hunfolded");
			hunfolded2[0]->Scale(scale_jets2[rat],"width");
			histo_set_atributes(hunfolded2[0], colorList[0+1], colorList[0+1], line_width, marker_size, markerS[unf2[0][rat]]);
			nsolutions2[rat]=1;
		}
		/* this part is not used anymore
		if(nsolutions12[rat]==0)
		{
			str = Form("iter%d",defiter[unf1[0][rat]]-1);
			cout<<"WARNING: no good solution => setting default iter to "<<defiter[unf1[0][rat]]<<endl;
			TDirectoryFile* dir1 = (TDirectoryFile*)funfolding1[0][rat]->Get(str.Data());
			hunfolded1_good12[0][rat] = (TH1D*) dir1->Get("hunfolded");
			hunfolded1_good12[0][rat]->Scale(scale_jets1[rat],"width");
			cout<<hunfolded1_good12[0][rat]->Integral()<<endl;
			histo_set_atributes(hunfolded1_good12[0][rat], colorList[0+1], colorList[0+1], line_width, marker_size, markerS[unf1[0][rat]]);

			if(ratio_type==RAA && ppBase==PYTHIA)
			{
				hunfolded2_good12[0][rat]=hunfolded2[0];
				nsolutions12[rat]=1;
			}
			else
			{
				str = Form("iter%d",defiter[unf2[0][rat]]-1);
				cout<<"WARNING: no good solution => setting default iter to "<<defiter[unf2[0][rat]]<<endl;
				TDirectoryFile* dir2 = (TDirectoryFile*)funfolding2[0][rat]->Get(str.Data());
				hunfolded2_good12[0][rat] = (TH1D*) dir2->Get("hunfolded");
				hunfolded2_good12[0][rat]->Scale(scale_jets2[rat],"width");
				histo_set_atributes(hunfolded2_good12[0], colorList[0+1], colorList[0+1], line_width, marker_size, markerS[unf2[0][rat]]);
				nsolutions12[rat]=1;
			}
		}*/
		
		//==============================
		//calculate average and spread
		hunfolded1_avg[rat]=(TH1D*) hunfolded1[0]->Clone(Form("haverage1_rat%i",rat));
		hunfolded1_avg[rat]->Reset("MICE");
		TH1D* hsigmap1=(TH1D*) hunfolded1_avg[rat]->Clone(Form("sigp1_rat%i",rat));
		hsigmap1->Reset("MICE");
		TH1D* hsigman1=(TH1D*) hunfolded1_avg[rat]->Clone(Form("sign1_rat%i",rat));
		hsigman1->Reset("MICE");
		calculate_avrg_spread(nsolutions1[rat], hunfolded1, hunfolded1_avg[rat], hsigmap1, hsigman1, err_norm_spec);
		hunfolded1_avg[rat]->SetLineColor(kRed);
		hsigmap1->SetLineColor(colorList3[syserr_type]);
		hsigman1->SetLineColor(colorList3[syserr_type]);
		
		hunfolded2_avg[rat]=(TH1D*) hunfolded2[0]->Clone(Form("haverage2_rat%i",rat));
		hunfolded2_avg[rat]->Reset("MICE");
		TH1D* hsigmap2=(TH1D*) hunfolded2_avg[rat]->Clone(Form("sigp2_rat%i",rat));
		hsigmap2->Reset("MICE");
		TH1D* hsigman2=(TH1D*) hunfolded2_avg[rat]->Clone(Form("sign2_rat%i",rat));
		hsigman2->Reset("MICE");
		calculate_avrg_spread(nsolutions2[rat], hunfolded2, hunfolded2_avg[rat], hsigmap2, hsigman2, err_norm_spec);
		hunfolded2_avg[rat]->SetMarkerColor(kRed);
		hunfolded2_avg[rat]->SetLineColor(kRed);
		hunfolded2_avg_scaled[rat]=(TH1D*) hunfolded2_avg[rat]->Clone(Form("hunfolded2_avg_scaled_%i",rat));
		hunfolded2_avg_scaled[rat]->Scale(scale_denum[rat]);
		hunfolded2_avg_scaled[rat]->SetLineColor(kBlue);
      hunfolded2_avg_scaled[rat]->SetMarkerColor(kBlue);
		
        //create TGraphAsymmErrors
        if(rat==0){
		short nbins_spec_tmp=calc_nbins(hunfolded1_avg[rat], spec_pTmin, spec_pTmax, firstbin_spec, 1);
		//if(printInfo)<<"# of bins to show in spectrum: "<<nbins_spec_tmp<<" firstbin:"<<firstbin_spec<<endl;
		const int nbins_spec=nbins_spec_tmp;
		double x_spec[nbins_spec];
		double x_specBSF[nbins_spec]; //for bin shift correction
		double x_specBSL[nbins_spec]; //for bin shift correction
		double x_spec_left[nbins_spec];
		double x_spec_right[nbins_spec];
		double x_spec_leftBSF[nbins_spec];
		double x_spec_rightBSF[nbins_spec];
		double x_spec_leftBSL[nbins_spec];
		double x_spec_rightBSL[nbins_spec];
		double y_spec[nbins_spec];
		double y_spec_up[nbins_spec];
		double y_spec_down[nbins_spec];
		double y_spec_sysunf_up[nbins_spec];
		double y_spec_sysunf_down[nbins_spec];
        
		//for denominator
		double x_spec2BSL[nbins_spec];
		double x_spec2_leftBSL[nbins_spec];
		double x_spec2_rightBSL[nbins_spec];
		double y_spec2[nbins_spec];
		double y_spec2_up[nbins_spec];
		double y_spec2_down[nbins_spec];
        
		for(int bn=0; bn<nbins_spec; bn++)
		{
			ushort hbin=bn+firstbin_spec;
			y_spec_up[bn]=hunfolded1_avg[rat]->GetBinError(hbin);
			y_spec_down[bn]=hunfolded1_avg[rat]->GetBinError(hbin);
			x_spec[bn]=hunfolded1_avg[rat]->GetBinCenter(hbin);
			y_spec[bn]=hunfolded1_avg[rat]->GetBinContent(hbin);
			x_spec_right[bn]=hunfolded1_avg[rat]->GetBinWidth(hbin)/2.0;
			x_spec_left[bn]=hunfolded1_avg[rat]->GetBinWidth(hbin)/2.0;
			
			//arrays for bin shift correction
			x_specBSF[bn]=hunfolded1_avg[rat]->GetBinCenter(hbin);
			x_specBSL[bn]=hunfolded1_avg[rat]->GetBinCenter(hbin);
			x_spec_rightBSF[bn]=hunfolded1_avg[rat]->GetBinWidth(hbin)/2.0;
			x_spec_leftBSF[bn]=hunfolded1_avg[rat]->GetBinWidth(hbin)/2.0;
			x_spec_rightBSL[bn]=hunfolded1_avg[rat]->GetBinWidth(hbin)/2.0;
			x_spec_leftBSL[bn]=hunfolded1_avg[rat]->GetBinWidth(hbin)/2.0;
			
			y_spec_sysunf_up[bn]=hsigmap1->GetBinContent(hbin);
			y_spec_sysunf_down[bn]=hsigman1->GetBinContent(hbin);
			
			//denominator
			x_spec2BSL[bn]=hunfolded2_avg_scaled[rat]->GetBinCenter(hbin);
			x_spec2_leftBSL[bn]=hunfolded2_avg_scaled[rat]->GetBinWidth(hbin)/2.0;
			x_spec2_rightBSL[bn]=hunfolded2_avg_scaled[rat]->GetBinWidth(hbin)/2.0;
			y_spec2[bn]=hunfolded2_avg_scaled[rat]->GetBinContent(hbin);
			y_spec2_up[bn]=hunfolded2_avg_scaled[rat]->GetBinError(hbin);
			y_spec2_down[bn]=hunfolded2_avg_scaled[rat]->GetBinError(hbin);
            
        }
        
		//fill average value into a tgraph
		gspec_avg=new TGraphAsymmErrors(nbins_spec, x_spec, y_spec, x_spec_left, x_spec_right, y_spec_down, y_spec_up);
		tgraph_set_atributes(gspec_avg, 0, spec_pTmax, spectraYmin, spectraYmax, kBlue,kBlue,kBlue, 1, 2, 29, 0);
		
		//create TGraph with unfolding errors
		gspec_sys=new TGraphAsymmErrors(nbins_spec, x_spec, y_spec, x_spec_left, x_spec_right, y_spec_sysunf_down, y_spec_sysunf_up);
        
		//BIN SHIFT CORRECTION (fit)
		FitF1=fit_for_binshift(R1_idx,C_TAA[system1]);
		bin_shift_LR(FitF1, nbins_spec, x_specBSF, y_spec, x_spec_leftBSF, x_spec_rightBSF, y_spec_down, y_spec_up, 3);
		gspec1_BSF=new TGraphAsymmErrors(nbins_spec, x_specBSF, y_spec, x_spec_leftBSF, x_spec_rightBSF, y_spec_down, y_spec_up);
		tgraph_set_atributes(gspec1_BSF, 0, spec_pTmax, spectraYmin, spectraYmax, kGreen+1,kGreen+1,kGreen+1, 1, 2, 29, 0);

		//BIN SHIFT CORRECTION (interpolation)
		bin_shift_LR(nbins_spec, x_specBSL, y_spec, x_spec_leftBSL, x_spec_rightBSL, y_spec_down, y_spec_up, 2);
		gspec1_BSL=new TGraphAsymmErrors(nbins_spec, x_specBSL, y_spec, x_spec_leftBSL, x_spec_rightBSL, y_spec_down, y_spec_up);
		tgraph_set_atributes(gspec1_BSL, 0, spec_pTmax, spectraYmin, spectraYmax, kBlue,kBlue,kBlue, 1, 2, 34, 0);
		//denominator
		bin_shift_LR(nbins_spec, x_spec2BSL, y_spec2, x_spec2_leftBSL, x_spec2_rightBSL, y_spec2_down, y_spec2_up, 2);
		gspec2_BSL=new TGraphAsymmErrors(nbins_spec, x_spec2BSL, y_spec2, x_spec2_leftBSL, x_spec2_rightBSL, y_spec2_down, y_spec2_up);
		tgraph_set_atributes(gspec2_BSL, 0, spec_pTmax, spectraYmin, spectraYmax, kBlue,kBlue,kBlue, 1, 2, 34, 0);
		
		//calculate ratio taking into account bin shift correction
		gratio_BSL=ratio_of_graphs(gspec1_BSL,gspec2_BSL);
		tgraph_set_atributes(gratio_BSL, 0, ratio_pTmax, ratioYmin, ratioYmax, kBlue,kBlue,kBlue, 1, 2, 29, 0);

        }//rat=0
        
        
		//==============================
		//calculate ratios
		
		hratio_avg[rat]=(TH1D*) hunfolded1_avg[rat]->Clone(Form("hratio_avg_%i",rat));
		//int bin=hratio_avg[rat]->GetXaxis()->FindBin(15);
		//cout<<"average:"<<endl;
		//cout<<"y1:"<<hratio_avg[rat]->GetBinContent(bin)<<endl;
		//hratio_avg[rat]->Scale(1.0/scale_denum[rat]);
		//cout<<"y2:"<<hratio_avg[rat]->GetBinContent(bin)<<endl;
		if(ratio_type==sys)
			hratio_avg[rat]->Divide(hunfolded2_avg[0]);
		else 
			hratio_avg[rat]->Divide(hunfolded2_avg_scaled[rat]);
		histo_set_atributes(hratio_avg[rat], colorList2[rat], colorList2[rat], line_width, marker_size, markerSsys[rat]);
		
		if(system==pp && rat==0 && showBBB)
		{
			TH1D* hbbb_ratio=(TH1D*) hbbb_rebin->Clone("hbbb_ratio");
			hbbb_ratio->Divide(hunfolded2_avg_scaled[0]);
			hbbb_ratio->SetLineColor(kGreen+1);
			hbbb_ratio->SetLineWidth(2);
		}

		
		//calculate assymetric statistical errors for ratio and save them into TGraphAsymmErrors
		if(rat==0){
		short nbins_tmp=calc_nbins(hratio_avg[rat], ratio_pTmin, ratio_pTmax, firstbin, 1);
		//if(printInfo)<<"# of bins to show: "<<nbins_tmp<<" firstbin:"<<firstbin<<endl;
		const int nbins=nbins_tmp;
		double x_rat[nbins];
		double x_rat_left[nbins];
		double x_rat_right[nbins];
		double y_rat[nbins];
		double y_MPVoverMU[nbins];
		double y_rat_up[nbins];
		double y_rat_down[nbins];
		
		TCanvas* cfit=new TCanvas("cfit","cfit",10,10,1600,1000);
		cfit->Divide(4,4);
		cfit->cd(1);
		for(int bn=0; bn<nbins; bn++)
		{
			cfit->cd(1+bn);
			ushort hbin=bn+firstbin;
			
			
			Double_t mu1,mu2,sigma1,sigma2,sigma_up,sigma_down;
			mu1=hunfolded1_avg[rat]->GetBinContent(hbin);
			mu2=hunfolded2_avg_scaled[rat]->GetBinContent(hbin);
			sigma1=hunfolded1_avg[rat]->GetBinError(hbin);
			sigma2=hunfolded2_avg_scaled[rat]->GetBinError(hbin);
			double mu=(mu2>0) ? (mu1/mu2) : 0;
			
			Double_t left=mu;
			Double_t right=mu;
			for(int i=0; i<100; i++)
			{
				Double_t rnd1=gRandom->Gaus(mu1,sigma1);
				Double_t rnd2=gRandom->Gaus(mu2,sigma2);
				Double_t xi=0;
				xi=(TMath::Abs(rnd2)>0) ? rnd1/rnd2 : 0;
				left=(xi<left) ? xi : left;
				right=(xi>right) ? xi : right;
			}

			if(right>1.5) right=1.5; //eliminate too large intervals
			//cout<<"left:"<<left<<" right:"<<right<<endl;
			//calculate asymmetric errors of the ratio
			TH1D* haseratio=new TH1D("haseratio","",40,left,right);
			TF1* gauss1=new TF1("gauss1","[0]*TMath::Exp(-((x-[1])*(x-[1]))/(2*[2]*[2]))",mu,right);
			TF1* gauss2=new TF1("gauss2","[0]*TMath::Exp(-((x-[1])*(x-[1]))/(2*[2]*[2]))",left,mu);
			
			if(printInfo)cout<<bn+1<<" mu1, mu2: "<<mu1<<", "<<mu2<<" sigma1,2: "<<sigma1<<", "<<sigma2<<endl;
			
			double mpv;
			calc_sigma_ratio(mu1,mu2,sigma1,sigma2,sigma_up,sigma_down,mpv,5000,1, haseratio, gauss1, gauss2,left,right); 
			
			if(printInfo)cout<<"sigma_up:"<<sigma_up<<" sigma_down:"<<sigma_down<<endl;
			
			
			double center=(useMPV) ? mpv : mu;
			y_rat_up[bn]=sigma_up;
			y_rat_down[bn]=sigma_down;
			x_rat[bn]=hratio_avg[rat]->GetBinCenter(hbin);
			y_rat[bn]=center; //hratio_avg[rat]->GetBinContent(hbin);
			y_MPVoverMU[bn]=(mu>0) ? (mpv/mu) : 0;
			hratio_avg[rat]->SetBinContent(hbin,center);
			x_rat_right[bn]=hratio_avg[rat]->GetBinWidth(hbin)/2.0;
			x_rat_left[bn]=hratio_avg[rat]->GetBinWidth(hbin)/2.0;
			
			
			haseratio->DrawCopy("");
			gauss1->DrawClone("same");
			gauss2->DrawClone("same");
			
			delete haseratio;
			delete gauss1;
			delete gauss2;
			
		}//bin loop
		
		//fill average value into a tgraph
		gratio_avg=new TGraphAsymmErrors(nbins, x_rat, y_rat, x_rat_left, x_rat_right, y_rat_down, y_rat_up);
		tgraph_set_atributes(gratio_avg, 0, ratio_pTmax, ratioYmin, ratioYmax, kBlue,kBlue,kBlue, 1, 2, 29, 0);

		//copy statistical error bars from ratio_avg to ratio_BSL
		copy_graph_errors(gratio_BSL, gratio_avg);
		//rescale the binshifted ratio graph from Mean to MPV
		if(useMPV) scale_graph(gratio_BSL,y_MPVoverMU);  
		
		//make ratio histograms for different unfolding solutions - these will be used to calculate systematic error from unfolding
		for(int isol=0; isol<nsolutions1[rat]; isol++){
		for(int jsol=0; jsol<nsolutions2[rat]; jsol++)
		{
			//cout<<"solution:"<<isol<<"/"<<nsolutions12[rat]<<endl;
			hratio_var[isol][jsol]=(TH1D*) hunfolded1[isol]->Clone(Form("hratio_var_%i_%i",isol,jsol));
			//cout<<"y1:"<<hratio_var[isol]->GetBinContent(bin)<<endl;
			hratio_var[isol][jsol]->Scale(1.0/scale_denum[rat]);
			//cout<<"y2:"<<hratio_var[isol]->GetBinContent(bin)<<endl;
			hratio_var[isol][jsol]->Divide(hunfolded2[jsol]);
			//cout<<"y3:"<<hratio_var[isol]->GetBinContent(bin)<<endl;
			
			hratio_var[isol][jsol]->SetLineColor(colorList[isol]);
			hratio_var[isol][jsol]->SetMarkerColor(colorList[jsol]);
			hratio_var[isol][jsol]->SetMarkerStyle(markerS[0]);
			hratio_var[isol][jsol]->SetMarkerSize(marker_size);
			hratio_var[isol][jsol]->SetLineWidth(line_width);
		}}
		}//rat==0
		//=================================
		//Draw histograms
		//=================================
		cout<<"DRAWING HISTOGRAMS"<<endl;
		
		//Spectrum numerator
		str=Form("cspectrum_%i",rat);
		cspectrum1[rat]= new TCanvas(str,str,10,10,2*can_x,can_y);
		cspectrum1[rat]->Divide(2,1);
		cspectrum1[rat]->cd(1);
		gPad->SetGrid();
		gPad->SetLogy();
		histo_set_atributes(frame, Form("Spectra, %s",sdesc1.Data()), "p_{T, corr}^{charged} (GeV/c)", "1/N_{events} dN_{jet}/(dp_{T}d#eta 2#pi)", 0, spec_pTmax, 1E-9, 1E0);
		frame->DrawCopy("");
		for(int isol=0;isol<nsolutions1[rat];isol++)
		{
			hunfolded1[isol]->DrawCopy("esame");
		}
		hunfolded1_avg[rat]->DrawCopy("esame");
		hunfolded2_avg_scaled[rat]->DrawCopy("esame");

		//draw systematic errors
		TH1D* hsigmap1_drw=(TH1D*)  hunfolded1_avg[rat]->Clone("hsigmap1_drw");
		TH1D* hsigman1_drw=(TH1D*)  hunfolded1_avg[rat]->Clone("hsigman1_drw");
		hsigmap1_drw->Add(hsigmap1,1);
		hsigman1_drw->Add(hsigman1,-1);
		hsigmap1_drw->SetLineColor(colorList3[0]);
		hsigman1_drw->SetLineColor(colorList3[0]);		
		hsigmap1_drw->SetMarkerColor(colorList3[0]);
		hsigman1_drw->SetMarkerColor(colorList3[0]);	
		hsigmap1_drw->DrawCopy("esame");
		hsigman1_drw->DrawCopy("esame");
		
		//Spectrum denominator
        /*
		str=Form("cspectrum2_%i",rat);
		cspectrum2[rat]= new TCanvas(str,str,10,10,can_x,can_y);
		cspectrum2[rat]->cd();
		cspectrum2[rat]->SetGrid();
		cspectrum2[rat]->SetLogy();*/
		cspectrum1[rat]->cd(2);
		gPad->SetGrid();
		gPad->SetLogy();
		histo_set_atributes(frame, Form("Spectra, %s",sdesc2.Data()), "p_{T, corr}^{charged} (GeV/c)", "1/N_{events} dN_{jet}/(dp_{T}d#eta 2#pi)", 0, spec_pTmax, 1E-9, 1E0);
		frame->DrawCopy("");
		for(int isol=0;isol<nsolutions2[rat];isol++)
		{
			hunfolded2[isol]->DrawCopy("esame");
		}
		hunfolded2_avg[rat]->DrawCopy("esame");
		if(ratio_type==RAA && ppBase==PYTHIA) hpythia[rat]->DrawCopy("same");
		//draw sys errors
		TH1D* hsigmap2_drw=(TH1D*)  hunfolded2_avg[rat]->Clone("hsigmap2_drw");
		TH1D* hsigman2_drw=(TH1D*)  hunfolded2_avg[rat]->Clone("hsigman2_drw");
		hsigmap2_drw->Add(hsigmap2,1);
		hsigman2_drw->Add(hsigman2,-1);
		hsigmap2_drw->SetLineColor(colorList3[0]);
		hsigman2_drw->SetLineColor(colorList3[0]);
		hsigmap2_drw->SetMarkerColor(colorList3[0]);
		hsigman2_drw->SetMarkerColor(colorList3[0]);
		hsigmap2_drw->DrawCopy("esame");
		hsigman2_drw->DrawCopy("esame");


		for(int isol=1;isol<nsolutions_max;isol++)
		{
			funfolding1[isol][rat]->Close();
			if(ratio_type!=RAA || ppBase!=PYTHIA)  funfolding2[isol][rat]->Close();
		}

	}// nratios loop
	

	
    
	//=========================================
	//calculate and draw systematic errors of spectra
	//=========================================
	//unfolding errors are already calculated, so here we calculate only correlated errors
	//here we also draw the errors as a TGraph - either unfolding errors or correlated errors, depending on syserr_type choice
    
	TH1D* hspecavg=(TH1D*) hunfolded1_avg[0]->Clone("hspecavg");
	hspecavg->Reset("MICE");
	TH1D* hspecsigmap=(TH1D*) hunfolded1_avg[0]->Clone("hspecsigmap");
	hspecsigmap->Reset("MICE");
	TH1D* hspecsigman=(TH1D*) hunfolded1_avg[0]->Clone("hspecsigman");
	hspecsigman->Reset("MICE");
    
	
	TCanvas *cspec=new TCanvas("cspec","spectrum_and_sys_errors",10,10,can_x,2*can_y);
	cspec->Divide(1,2);
	cspec->cd(1);
	//gPad->SetGrid();
	gPad->SetLogy();
	histo_set_atributes(frame, Form("Spectrum, %s",sdesc1.Data()), "p_{T, corr}^{charged} (GeV/c)", "1/N_{events} dN_{jet}/(dp_{T}d#eta 2#pi)", 0, spec_pTmax, 1E-9, 1E0);
	frame->DrawCopy("");
	
	
	if(useaverages4sys)//correlated errors
	{
		//draw averages
		for(int rat=0; rat<nratios; rat++)
		{
			hunfolded1_avg[rat]->DrawCopy("esame");
		}
		
		//calculate sys errors from spread between averages
		if(nratios>1) 
		{
			calculate_avrg_spread(nratios, hunfolded1_avg, hspecavg, hspecsigmap, hspecsigman, err_norm);
			
		}
		//draw sys errors
		TH1D* hspecsigmap_drw=(TH1D*)  hunfolded1_avg[0]->Clone("hspecsigmap_drw");
		TH1D* hspecsigman_drw=(TH1D*)  hunfolded1_avg[0]->Clone("hspecsigman_drw");
		hspecsigmap_drw->Add(hspecsigmap,1);
		hspecsigman_drw->Add(hspecsigman,-1);
		hspecsigmap_drw->SetLineColor(colorList3[syserr_type]);
		hspecsigman_drw->SetLineColor(colorList3[syserr_type]);
		hspecsigmap_drw->SetMarkerColor(colorList3[syserr_type]);
		hspecsigman_drw->SetMarkerColor(colorList3[syserr_type]);
		hspecsigmap_drw->DrawCopy("esame");
		hspecsigman_drw->DrawCopy("esame");
		
		//save the sys. error into arrays (to be used with TGraph)
		double y_spec_sys_up[nbins];
		double y_spec_sys_down[nbins];
		for(int bn=0; bn<nbins; bn++)
		{
			ushort hbin=bn+firstbin_spec;
			y_spec_sys_up[bn]=hspecsigmap->GetBinContent(hbin);
			y_spec_sys_down[bn]=hspecsigman->GetBinContent(hbin);
		}
	
		//create TGraph
		gspec_sys=new TGraphAsymmErrors(nbins, x_spec, y_spec, x_spec_left, x_spec_right, y_spec_sys_down, y_spec_sys_up);
		
	} //correlated errors 
	else //unfolding errors, normalization errors
	{
		hunfolded1_avg[0]->DrawCopy("esame");		 
	}//unfolding errors
	tgraph_set_atributes(gspec_sys, 0, spec_pTmax, spectraYmin, spectraYmax, colorList3[syserr_type],colorList3[syserr_type],colorList3[syserr_type], 1, 2, 29, 0);

	
	//Draw spectrum as a TGraph with assymetrical errors + systematic error also as TGraph
	cspec->cd(2);
	gPad->SetLogy();
	gspec_avg->DrawClone("AP");
	gspec1_BSF->DrawClone("P");
	gspec1_BSL->DrawClone("P");
	gspec_sys->DrawClone("2");
	if(ratio_type==closure && doToy)
	{
		gtoytrue->DrawClone("P");
		int nskipbins=0;
		int mlp=1E6;
		double intTrue=mlp*graph_integral(gtoytrue,none,5+nskipbins);
		double intUnf=mlp*graph_integral(gspec_sys,none,nskipbins);
		double intUnf_up=mlp*graph_integral(gspec_sys,up,nskipbins);
		double intUnf_down=mlp*graph_integral(gspec_sys,down,nskipbins);
		double errY_up=intUnf_up-intUnf;
		double errY_down=intUnf-intUnf_down;
		
		cout<<"Yield of unfolded spectra (x"<<mlp<<"): "<<intUnf<<"+"<<errY_up<<"-"<<errY_down<<endl;
		cout<<"Yiled of generated spectra (x"<<mlp<<"): "<<intTrue<<endl;
		
	}
		
	//draw also fit function
	FitF1->Draw("same");
   
    
    
	
	//=========================================
	//draw RATIO
	//=========================================
	float mplier=(ratio_type==sys) ? 1.2 : 2.0;
	TCanvas *cratio=new TCanvas("cratio","cratio",10,10,can_x,mplier*can_y);
	if(ratio_type!=sys)
	{
		cratio->Divide(1,2);
		cratio->cd(1);
	}
	histo_set_atributes(frame, "", "p_{T, corr}^{charged} (GeV/c)", y_title, ratioXmin, ratioXmax, ratioYmin,ratioYmax);
	frame->SetTitle(rdesc);
	//if(ratio_type==sys)frame->SetTitle("");
	frame->DrawCopy("");
    
    
	//=========================================
	//calculate and draw systematic errors of ratio 
	//=========================================
	//here we can calculate all types of systematic erroros - from unfolding, corelated errors, from normalization
    
	TH1D* hratavg=(TH1D*) hratio_avg[0]->Clone("hratavg");  //copy of average which can be replaced by the averge calculated by calculate_avrg_spread(....,calculate_avg=true)
	//hratavg->Reset("MICE");
	TH1D* hratsigmap=(TH1D*) hratio_avg[0]->Clone("hratsigmap");
	hratsigmap->Reset("MICE");
	TH1D* hratsigman=(TH1D*) hratio_avg[0]->Clone("hratsigman");
	hratsigman->Reset("MICE");
	
	TLegend *legend1 = new TLegend(0.2, 0.65, 0.35, 0.90);
	legend1->SetTextSize(0.03);
	legend1->SetFillStyle(0);
	legend1->SetBorderSize(0);

	TLegend *legend2 = new TLegend(0.2, 0.2, 0.35, 0.4);
	legend2->SetTextSize(0.03);
	legend2->SetFillStyle(0);
	legend2->SetBorderSize(0);

	ofstream sysErrFile;
	if(printTextFile)
	{
		sysErrFile.open(Form("sysErr/%s/%s/bin%i/sysErr_%s_R0%.0lf_pTl%.0lf.%s",dir_suffix1.Data(),syserr_name[syserr_type].Data(),bining,system_desc_short[system1].Data(),R1*10,pTl1,"txt"));
		sysErrFile<<"======================================"<<endl;
	}

	if(useaverages4sys) //correlated errors
	{
		//print out the size of the relative systematic error at pT=pTshowsys
		if(printTextFile)
		{
		for(int s=0;s<nshowsys;s++)
		{
			int binshow=hratio_avg[0]->FindBin(pTshowsys[s]);
			float val_base=hratio_avg[0]->GetBinContent(binshow);
			for(int rat=0; rat<nratios; rat++)
			{
				if(rat==0) 
				{
					sysErrFile<<"---------------------------------------"<<endl;
					sysErrFile<<"relative systematic error at "<<pTshowsys[s]<<" GeV/c:"<<endl;
				}
				float val_err=hratio_avg[rat]->GetBinContent(binshow);
				float rel_err=(val_base>0) ? ((val_err-val_base)/val_base) : 0;
				rel_err=rel_err*100; //in [%]
				sysErrFile<<leg_name[rat].Data()<<": "<<rel_err<<"%"<<endl;
			}
		}
		}
		
		//draw histograms create legend
		for(int rat=0; rat<nratios; rat++)
		{
			hratio_avg[rat]->DrawCopy("esame");
			if(rat<7)
				legend1->AddEntry(hratio_avg[rat], Form("%s",leg_name[rat].Data()), "lp");
			else
				legend2->AddEntry(hratio_avg[rat], Form("%s",leg_name[rat].Data()), "lp");
		}
		
		if(nratios>1) //from spread between averages
		{
			//cout<<"using normalization:"<<err_norm<<endl;
			calculate_avrg_spread(nratios-1, hratio_avg, hratavg, hratsigmap, hratsigman, err_norm,false,1);
			
		}
		else //from systematic error of numerator and denominator
		{
			//TODO
			for(int bn=1; bn<=hsigman1->GetNBinsX(); bn++)
			{
				float sp1=hsigmap1->GetBinContent(bn);
				float sn1=hsigman1->GetBinContent(bn);
				float sp2=hsigmap2->GetBinContent(bn);
				float sn2=hsigman2->GetBinContent(bn);
				float spr;
				float snr;
			}
		}
	}
	else //unfolding errors - calculate systematic error from spread of different solutions
	{
		calculate_avrg_spread2D(nsolutions1[0],nsolutions2[0], hratio_var, hratavg, hratsigmap, hratsigman, err_norm);
		//calculate_avrg_spread(nsolutions12[0],hratio_var, hratavg, hratsigmap, hratsigman, err_norm);
		//the main value we will not use hratavg (average calculated from spread of different ratios) but rather hratio_avg[0] (ratio of averages)
		hratio_avg[0]->DrawCopy("esame");
		legend1->AddEntry(hratio_avg[0], "ratio average", "lp");
		if(system==pp && showBBB)
		{
			hbbb_ratio->DrawCopy("same");
			legend1->AddEntry(hbbb_ratio, "bin-by-bin", "lp");
		}
		int cnt=0;
		for(int isol=0; isol<nsolutions1[0]; isol++){
		for(int jsol=0; jsol<nsolutions2[0]; jsol++)
		{
			cnt++;
			//if(cnt>50)continue;
			hratio_var[isol][jsol]->DrawCopy("esame");
			legend1->AddEntry(hratio_var[isol][jsol], Form("ratio unf_%i_%i",isol,jsol), "lp");
		}}
	}

	
	//draw sys errors
	TH1D* hratsigmap_drw=(TH1D*)  hratio_avg[0]->Clone("hratsigmap_drw");
	TH1D* hratsigman_drw=(TH1D*)  hratio_avg[0]->Clone("hratsigman_drw");
	hratsigmap_drw->Add(hratsigmap,1);
	hratsigman_drw->Add(hratsigman,-1);
	hratsigmap_drw->SetLineColor(colorList3[syserr_type]);
	hratsigman_drw->SetLineColor(colorList3[syserr_type]);
	hratsigmap_drw->SetMarkerColor(colorList3[syserr_type]);
	hratsigman_drw->SetMarkerColor(colorList3[syserr_type]);
	hratsigmap_drw->DrawCopy("esame");
	hratsigman_drw->DrawCopy("esame");
	legend2->AddEntry(hratsigmap_drw, "total sys. error", "lp");

	legend1->DrawClone("same");
	legend2->DrawClone("same");
    
	//printing out the relative size of the total sys errot at pT=pTshowsys
	if(printTextFile)
	{
		sysErrFile<<"======================================"<<endl;
		for(int s=0;s<nshowsys;s++)
		{
			int binshow=hratio_avg[0]->FindBin(pTshowsys[s]);
			float val_base=hratio_avg[0]->GetBinContent(binshow);
			float val_errp=hratsigmap_drw->GetBinContent(binshow);
			float rel_errp=(val_base>0) ? ((val_errp-val_base)/val_base) : 0;
			rel_errp=rel_errp*100; //in [%]
			float val_errn=hratsigman_drw->GetBinContent(binshow);
			float rel_errn=(val_base>0) ? ((val_errn-val_base)/val_base) : 0;
			rel_errn=rel_errn*100; //in [%]
			
			sysErrFile<<"relative systematic error at "<<pTshowsys[s]<<" GeV/c:"<<endl;
			sysErrFile<<"total sys. error:"<<rel_errp<<"%"<<endl;
			sysErrFile<<"total sys. error:"<<rel_errn<<"%"<<endl;
		}
		sysErrFile.close();
	}
	
    TLine *one = new TLine(0,yline, ratio_pTmax, yline);
    one->SetLineWidth(2);
    one->SetLineStyle(2);
    one->SetLineColor(kBlack);
    one->DrawClone("same");
	
	//save the sys. error into arrays (to be used with TGraph)
	double y_rat_sys_up[nbins];
	double y_rat_sys_down[nbins];
	for(int bn=0; bn<nbins; bn++)
	{
		ushort hbin=bn+firstbin;
		y_rat_sys_up[bn]=hratsigmap->GetBinContent(hbin);
		y_rat_sys_down[bn]=hratsigman->GetBinContent(hbin);
	}
	
	//create TGraph with systematic errors
	gratio_sys=new TGraphAsymmErrors(nbins, x_rat, y_rat, x_rat_left, x_rat_right, y_rat_sys_down, y_rat_sys_up);
	tgraph_set_atributes(gratio_sys, 0, ratio_pTmax, ratioYmin, ratioYmax, colorList3[syserr_type],colorList3[syserr_type],colorList3[syserr_type], 1, 2, 29, 0);
   
	//rescale the systematic error for the new position of points calculated using bincenter correction
	gratio_sys_BSL=new TGraphAsymmErrors(nbins, x_rat, y_rat, x_rat_left, x_rat_right, y_rat_sys_down, y_rat_sys_up);
	tgraph_set_atributes(gratio_sys_BSL, 0, ratio_pTmax, ratioYmin, ratioYmax, colorList3[syserr_type],colorList3[syserr_type],colorList3[syserr_type], 1, 2, 29, 0);
	move_graph_points(gratio_sys_BSL, gratio_BSL);


	
	//Draw ratio as a TGraph with assymetrical errors + systematic error also as TGraph
	if(ratio_type!=sys)
	{
	cratio->cd(2);
	gratio_avg->DrawClone("AP");
	gratio_BSL->DrawClone("P");
	gratio_sys->DrawClone("2");
	gratio_sys_BSL->DrawClone("2");
	one->DrawClone("same");
	if(system==pp && showBBB)
		{
			hbbb_ratio->DrawCopy("same");
		}
	}

	//save ratio plot as a figure
	if(ratio_type==sys)
	{
		TString sys_out_path=Form("../../../plotting_out/obr/intersteps/sys_errors/%s/bin%i",syserr_name[syserr_type].Data(),bining);
		TString ext[2]={"gif","pdf"};
		for(int e=0; e<2; e++)
		{
			cratio->SaveAs(Form("%s/sysErr_%s_R0%.0lf_pTl%.0lf.%s",sys_out_path.Data(),system_desc_short[system1].Data(),R1*10,pTl1,ext[e].Data()));
		}
	}
	
	//Save TGraphs to external file    
	if(saveSysFile)
	{
		TFile* fout=create_sys_file(system, R, pTlead,bining, dir_suffix1, doToy,ppBaseWpTlead);
		fout->cd();
		if(syserr_type==unfolding) //to make sure we write the graphs only once, not 3times (unfolding, correlated, normalization)
		{
			gratio_avg->Write(ratio_name[ratio_type]);
			gratio_BSL->Write(Form("%s_BSL",ratio_name[ratio_type].Data()));
			if(ratio_type==RAA) //to make sure we write the graphs only once, not 4times (RAA,RCP,RR,pTlpTl)
			{
				gspec_avg->Write("spectrum");
				gspec1_BSF->Write("spectrum_BSF");
				gspec1_BSL->Write("spectrum_BSL");
				FitF1->Write("spectrum_fit");
				hpythia[0]->Write("pythia");
			}
			else if(ratio_type==closure && doToy)
			{	
				gtoytrue->Write("spectrum_toy_true");
			}
		}
		gratio_sys->Write(Form("%s_%s",ratio_name[ratio_type].Data(),syserr_name[syserr_type].Data()));
		gratio_sys_BSL->Write(Form("%s_%s_BSL",ratio_name[ratio_type].Data(),syserr_name[syserr_type].Data()));
		

		if(ratio_type==RAA) //to make sure we write the given systematic error for the spectrum only once, not 4times (RAA,RCP,RR,pTlpTl)
			gspec_sys->Write(Form("spectrum_%s",syserr_name[syserr_type].Data())); 

		
		fout->Close();
	}
	
	
	
	return;
}

//========================================================================================================================================
//additional functions
//(see utils.C for other functions not declared here)
//========================================================================================================================================


TFile* create_sys_file(C_systems system=cent, float R=0.3, float pTlead=5.0,short bining=1, TString evolution="eta", bool doToy=0, /*tests test*/ bool pp_pTlead=0)
{
	TString sdir="../../../plotting_out/systematics";
	if(doToy)sdir+="/toymodel";
	if(system==peri) sdir+="/peripheral";
	else if(system==pp) sdir+="/pp";
	else sdir+="/central";
	sdir+=Form("/%s",evolution.Data());
	if(pp_pTlead) sdir+="_ppBiased";
	sdir+=Form("/bining%i",bining);
	//sdir+=Form("/%s",C_test_name[test].Data());
	TString str=Form("%s/systematics_R%.1lf_pT%0.lf.root",sdir.Data(),R,pTlead);
	TFile *foutput;
	foutput= new TFile(str.Data(), "UPDATE");
	return foutput;
}

TF1* fit_for_binshift(int ridx, float C_TAA)
{
		//TODO - system dependency
	
		/*
		Double_t B[3] = {2.41349e+01,1.75255e+01,1.27430e+01}; //R=0.2,0.3,0.4
      Double_t T[3] = {2.69190e+00,2.53251e+00,2.40521e+00};
      Double_t n[3] = {1.42495e+02,1.03738e+02,8.14043e+01};
      Double_t m0[3] = {-3.00000e+00,-7.83741e+00,-8.40985e+00};
      Double_t mu[3] = {-2.97634e+01,-2.01221e+01,-1.62948e+01};
      Double_t A[3] = {1.12732e+00,5.45759e-01,3.25500e-01};
      Double_t pwr[3] = {4.82216e+00,4.37727e+00,4.08808e+00};
      
      TF1 *fit_fnc=new TF1("fit_fnc",hardjet,3,40,9);
      fit_fnc->SetParameters(C_TAA,B[ridx],T[ridx],n[ridx],m0[ridx],mu[ridx],A[ridx],pwr[ridx],0.2);*/
		
		float A=0.1;
		float n=10;
		float T=1;
		TF1 *fit_fnc=new TF1("fit_fnc",TsalisFitFunc,3,40,3);
		fit_fnc->SetParameters(A,n,T);
      fit_fnc->SetNpx(10000);
		return fit_fnc;
}
 


