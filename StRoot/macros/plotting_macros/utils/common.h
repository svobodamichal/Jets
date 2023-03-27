#ifndef common_h
#define common_h

enum C_tests { chi2, camb, KS };
enum C_systems { cent, peri, pp };
enum C_unfoldings { Bayes, SVD };

const int C_nbinnings=5;
const int C_nsystems=3;
const TString C_system_name[C_nsystems]={"central","peripheral","pp"};
const int C_ntests=3;
const TString C_test_name[C_ntests]={"chi2","Camb","KS"};
const int C_nunfoldings=2;
const TString C_unfoldings_name[]={"Bayes","SVD"};
const int C_nR=3;
const float C_Rs[]={0.2,0.3,0.4,0.5}; //jet size R values
const int C_npTlead=8;
const float C_pTls[]={0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0}; //pTleading cut values
TString C_prior_type[]={/*0*/"truth",/*1*/"flat",/*2*/"pythia",/*3*/"powlaw3",/*4*/"powlaw45",/*5 */"powlaw5",/*6*/"powlaw55","tsalis_1","tsalis_2","tsalis_3","tsalis_4","tsalis_5","tsalis_6","tsalis_7","tsalis_8","tsalis_9"};
TString C_prior_type_name[]={"truth", "flat","biased Pythia","1/pT^{3}","1/pT^{4.5}","1/pT^{5}","1/(pT)^{5.5}","tsalis_1","tsalis_2","tsalis_3","tsalis_4","tsalis_5","tsalis_6","tsalis_7","tsalis_8","tsalis_9"};


//TAA and <Nbin>
const float C_ppXsection=30; //42 - zdc, 30 - vpd
const float C_TAA[/*systems*/]={22.75,0.49, 1.0/C_ppXsection}; //nuclear overlap function; 60-80%: TAA=0.49, 65-80%: TAA=0.37
const float C_TAA_err[/*systems*/]={1.56,0.14,0.01}; //error on TAA
const float C_NBIN[/*systems*/]={955.4,20.4,42.0/30.0}; //average number of binary collisions; 60-80%: 20.4, 65-80%: 15.77
const float C_NBIN_rerr[/*systems*/]={0.1,0.29,0.29}; //RELATIVE error on NBIN

//normalization relative errors on RAA and RCP
const float C_RAA_NormErr[]={/*central:*/0.07,/*peripheral:*/0.28};
const float C_RCP_NormErr=TMath::Sqrt(C_RAA_NormErr[0]*C_RAA_NormErr[0]+C_RAA_NormErr[1]*C_RAA_NormErr[1]);

//prior directory names
//TString prior_type[]={/*0*/"truth",/*1*/"flat",/*2*/"pythia",/*3*/"powlaw3",/*4*/"powlaw45",/*5 */"powlaw5",/*6*/"powlaw55","tsalis_1","tsalis_2","tsalis_3","tsalis_4","tsalis_5","tsalis_6","tsalis_7","tsalis_8","tsalis_9"};
//prior description
//TString prior_type_name[]={"truth", "flat","biased Pythia","1/pT^{3}","1/pT^{4.5}","1/pT^{5}","1/(pT)^{5.5}","tsalis_1","tsalis_2","tsalis_3","tsalis_4","tsalis_5","tsalis_6","tsalis_7","tsalis_8","tsalis_9"};

#endif
