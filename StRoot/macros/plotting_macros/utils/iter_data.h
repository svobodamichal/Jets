Int_t diter_data_normal[20][4][7][2][3][5][3]; //optimal iterations [prior][R][pTlead][Bayes/SVD][central/peripheral/pp][binning][test]
Int_t diter_data_peripheral_normal[20][4][7][2][3][5][3]; //optimal iterations [prior][R][pTlead][Bayes/SVD][central/peripheral/pp][binning][test]
Int_t diter_data_PP_normal[20][4][7][2][3][5][3]; //optimal iterations [prior][R][pTlead][Bayes/SVD][central/peripheral/pp][binning][test]

Int_t diter_data_p5[20][4][7][2][3][5][3]; //optimal iterations [prior][R][pTlead][Bayes/SVD][central/peripheral/pp][binning][test]
Int_t diter_data_m5[20][4][7][2][3][5][3]; //optimal iterations [prior][R][pTlead][Bayes/SVD][central/peripheral/pp][binning][test]
Int_t diter_data_u[20][4][7][2][3][5][3];
Int_t diter_data_g[20][4][7][2][3][5][3];
Int_t diter_data_nrem1[20][4][7][2][3][5][3];
Int_t diter_data_RRho02[20][4][7][2][3][5][3];
Int_t diter_data_RRho04[20][4][7][2][3][5][3];
Int_t diter_data_trcuts2[20][4][7][2][3][5][3];
Int_t diter_data_pp[20][4][7][2][3][5][3];
Int_t diter_data_v2[20][4][7][2][3][5][3];
Int_t diter_data_pythia[20][4][7][2][3][5][3];

Int_t diter_data_peripheral_p5[20][4][7][2][3][5][3]; //optimal iterations [prior][R][pTlead][Bayes/SVD][central/peripheral/pp][binning][test]
Int_t diter_data_peripheral_m5[20][4][7][2][3][5][3]; //optimal iterations [prior][R][pTlead][Bayes/SVD][central/peripheral/pp][binning][test]
Int_t diter_data_peripheral_u[20][4][7][2][3][5][3];
Int_t diter_data_peripheral_g[20][4][7][2][3][5][3];
Int_t diter_data_peripheral_nrem1[20][4][7][2][3][5][3];
Int_t diter_data_peripheral_RRho02[20][4][7][2][3][5][3];
Int_t diter_data_peripheral_RRho04[20][4][7][2][3][5][3];
Int_t diter_data_peripheral_trcuts2[20][4][7][2][3][5][3];
Int_t diter_data_peripheral_pp[20][4][7][2][3][5][3];
Int_t diter_data_peripheral_v2[20][4][7][2][3][5][3];
Int_t diter_data_peripheral_pythia[20][4][7][2][3][5][3];

/*
Int_t diter_data_p5[20][4][7][2][3][5][3]; //optimal iterations [prior][R][pTlead][Bayes/SVD][central/peripheral/pp][binning][test]
Int_t diter_data_m5[20][4][7][2][3][5][3]; //optimal iterations [prior][R][pTlead][Bayes/SVD][central/peripheral/pp][binning][test]

Int_t diter_data_AuAu[20][3][4][2][2][5][3]; //optimal iterations [prior][R][pTlead][Bayes/SVD][central/peripheral][binning]
Int_t diter_data_v2[20][3][4][2][2][5][3]; //optimal iterations [prior][R][pTlead][Bayes/SVD][central/peripheral][binning]
Int_t diter_data_2u1g[20][3][4][2][2][5][3]; //optimal iterations [prior][R][pTlead][Bayes/SVD][central/peripheral][binning]
Int_t diter_data_m5[20][3][4][2][2][5][3]; //optimal iterations [prior][R][pTlead][Bayes/SVD][central/peripheral][binning]
Int_t diter_data_p5[20][3][4][2][2][5][3]; //optimal iterations [prior][R][pTlead][Bayes/SVD][central/peripheral][binning]
*/
#include "optIter/optIter_data_normal_Camb_bin0.txt"
#include "optIter/optIter_data_normal_Camb_bin1.txt"
#include "optIter/optIter_data_normal_Camb_bin2.txt"
#include "optIter/optIter_data_normal_Camb_bin3.txt"
#include "optIter/optIter_data_normal_Camb_bin4.txt"
//#include "optIter/optIter_data_normal_chi2_bin1.txt"
//#include "optIter/optIter_data_normal_KS_bin1.txt"

#include "optIter/optIter_data_peripheral_normal_Camb_bin0.txt"
#include "optIter/optIter_data_peripheral_normal_Camb_bin1.txt"
#include "optIter/optIter_data_peripheral_normal_Camb_bin2.txt"
#include "optIter/optIter_data_peripheral_normal_Camb_bin3.txt"
#include "optIter/optIter_data_peripheral_normal_Camb_bin4.txt"
//#include "optIter/optIter_data_peripheral_normal_chi2_bin1.txt"
//#include "optIter/optIter_data_peripheral_normal_KS_bin1.txt"

/*
#include "optIter/optIter_data_PP_normal_Camb_bin1.txt"
#include "optIter/optIter_data_PP_normal_Camb_bin2.txt"
#include "optIter/optIter_data_PP_normal_Camb_bin3.txt"
#include "optIter/optIter_data_PP_normal_Camb_bin4.txt"
*/

#include "optIter/optIter_data_p5_Camb_bin0.txt"
#include "optIter/optIter_data_p5_Camb_bin1.txt"
#include "optIter/optIter_data_p5_Camb_bin4.txt"

#include "optIter/optIter_data_m5_Camb_bin0.txt"
#include "optIter/optIter_data_m5_Camb_bin1.txt"
#include "optIter/optIter_data_m5_Camb_bin4.txt"

#include "optIter/optIter_data_u_Camb_bin0.txt"
#include "optIter/optIter_data_u_Camb_bin1.txt"
#include "optIter/optIter_data_u_Camb_bin4.txt"

#include "optIter/optIter_data_g_Camb_bin0.txt"
#include "optIter/optIter_data_g_Camb_bin1.txt"
#include "optIter/optIter_data_g_Camb_bin4.txt"

#include "optIter/optIter_data_RRho02_Camb_bin0.txt"
#include "optIter/optIter_data_RRho02_Camb_bin1.txt"
#include "optIter/optIter_data_RRho02_Camb_bin4.txt"

#include "optIter/optIter_data_RRho04_Camb_bin0.txt"
#include "optIter/optIter_data_RRho04_Camb_bin1.txt"
#include "optIter/optIter_data_RRho04_Camb_bin4.txt"

#include "optIter/optIter_data_nrem1_Camb_bin0.txt"
#include "optIter/optIter_data_nrem1_Camb_bin1.txt"
#include "optIter/optIter_data_nrem1_Camb_bin4.txt"

#include "optIter/optIter_data_v2_Camb_bin0.txt"
#include "optIter/optIter_data_v2_Camb_bin1.txt"
#include "optIter/optIter_data_v2_Camb_bin4.txt"

#include "optIter/optIter_data_pythia_Camb_bin0.txt"
#include "optIter/optIter_data_pythia_Camb_bin1.txt"
#include "optIter/optIter_data_pythia_Camb_bin4.txt"

#include "optIter/optIter_data_trcuts2_Camb_bin0.txt"
#include "optIter/optIter_data_trcuts2_Camb_bin1.txt"
#include "optIter/optIter_data_trcuts2_Camb_bin4.txt"

#include "optIter/optIter_data_pp_Camb_bin0.txt"
#include "optIter/optIter_data_pp_Camb_bin1.txt"
#include "optIter/optIter_data_pp_Camb_bin4.txt"


#include "optIter/optIter_data_peripheral_p5_Camb_bin0.txt"
#include "optIter/optIter_data_peripheral_p5_Camb_bin1.txt"
#include "optIter/optIter_data_peripheral_p5_Camb_bin4.txt"

#include "optIter/optIter_data_peripheral_m5_Camb_bin0.txt"
#include "optIter/optIter_data_peripheral_m5_Camb_bin1.txt"
#include "optIter/optIter_data_peripheral_m5_Camb_bin4.txt"

#include "optIter/optIter_data_peripheral_u_Camb_bin0.txt"
#include "optIter/optIter_data_peripheral_u_Camb_bin1.txt"
#include "optIter/optIter_data_peripheral_u_Camb_bin4.txt"

#include "optIter/optIter_data_peripheral_g_Camb_bin0.txt"
#include "optIter/optIter_data_peripheral_g_Camb_bin1.txt"
#include "optIter/optIter_data_peripheral_g_Camb_bin4.txt"

#include "optIter/optIter_data_peripheral_RRho02_Camb_bin0.txt"
#include "optIter/optIter_data_peripheral_RRho02_Camb_bin1.txt"
#include "optIter/optIter_data_peripheral_RRho02_Camb_bin4.txt"

#include "optIter/optIter_data_peripheral_RRho04_Camb_bin0.txt"
#include "optIter/optIter_data_peripheral_RRho04_Camb_bin1.txt"
#include "optIter/optIter_data_peripheral_RRho04_Camb_bin4.txt"

#include "optIter/optIter_data_peripheral_nrem1_Camb_bin0.txt"
#include "optIter/optIter_data_peripheral_nrem1_Camb_bin1.txt"
#include "optIter/optIter_data_peripheral_nrem1_Camb_bin4.txt"

#include "optIter/optIter_data_peripheral_v2_Camb_bin0.txt"
#include "optIter/optIter_data_peripheral_v2_Camb_bin1.txt"
#include "optIter/optIter_data_peripheral_v2_Camb_bin4.txt"

#include "optIter/optIter_data_peripheral_pythia_Camb_bin0.txt"
#include "optIter/optIter_data_peripheral_pythia_Camb_bin1.txt"
#include "optIter/optIter_data_peripheral_pythia_Camb_bin4.txt"

#include "optIter/optIter_data_peripheral_trcuts2_Camb_bin0.txt"
#include "optIter/optIter_data_peripheral_trcuts2_Camb_bin1.txt"
#include "optIter/optIter_data_peripheral_trcuts2_Camb_bin4.txt"

#include "optIter/optIter_data_peripheral_pp_Camb_bin0.txt"
#include "optIter/optIter_data_peripheral_pp_Camb_bin1.txt"
#include "optIter/optIter_data_peripheral_pp_Camb_bin4.txt"
