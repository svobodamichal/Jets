Int_t diter_toy_normal[20][4][5][2][2][5][3];

#include "optIter/optIter_toy_normal_Camb_bin1.txt"
#include "optIter/optIter_toy_normal_chi2_bin1.txt"
#include "optIter/optIter_toy_normal_KS_bin1.txt"
#include "optIter/optIter_toy_normal_Camb_bin4.txt"
#include "optIter/optIter_toy_normal_chi2_bin4.txt"
#include "optIter/optIter_toy_normal_KS_bin4.txt"


#include "optIter/optIter_toy_normal_Camb_bin2.txt"
#include "optIter/optIter_toy_normal_chi2_bin2.txt"
#include "optIter/optIter_toy_normal_KS_bin2.txt"
#include "optIter/optIter_toy_normal_Camb_bin3.txt"
#include "optIter/optIter_toy_normal_chi2_bin3.txt"
#include "optIter/optIter_toy_normal_KS_bin3.txt"
/*

//prior|R|pTlead|unfolding|centrality|binning
//R=0.2
//pTlead=5
//Bayes
diter_toy[2][0][2][0][0][1]=3; //chi2bf=1.2, chi2ch=14.8
diter_toy[4][0][2][0][0][1]=4; //chi2bf=1.9, chi2ch=9.1
diter_toy[5][0][2][0][0][1]=4; //chi2bf=6.0, chi2ch=19.5
diter_toy[6][0][2][0][0][1]=6; //chi2bf=7.4, chi2ch=20.0
diter_toy[7][0][2][0][0][1]=-1;
diter_toy[8][0][2][0][0][1]=5; //chi2bf=0.8, chi2ch=18.2
diter_toy[9][0][2][0][0][1]=5; //chi2bf=4.7, chi2ch=18.6
diter_toy[10][0][2][0][0][1]=-1;
diter_toy[11][0][2][0][0][1]=6; //chi2bf=1.3, chi2ch=15.3
diter_toy[12][0][2][0][0][1]=4; //chi2bf=3.5, chi2ch=17.5
diter_toy[13][0][2][0][0][1]=-1;
diter_toy[14][0][2][0][0][1]=-1;
diter_toy[15][0][2][0][0][1]=7; //chi2bf=1.2, chi2ch=18.1
//SVD
diter_toy[2][0][2][1][0][1]=3; //chi2bf=2.6, chi2ch=296.1
diter_toy[4][0][2][1][0][1]=4; //chi2bf=4.1, chi2ch=140.2
diter_toy[5][0][2][1][0][1]=5; //chi2bf=1.3, chi2ch=163.1
diter_toy[6][0][2][1][0][1]=5; //chi2bf=1.8, chi2ch=221.5
diter_toy[7][0][2][1][0][1]=-1;
diter_toy[8][0][2][1][0][1]=4; //chi2bf=1.6, chi2ch=292.3
diter_toy[9][0][2][1][0][1]=-1;
diter_toy[10][0][2][1][0][1]=-1;
diter_toy[11][0][2][1][0][1]=-1;
diter_toy[12][0][2][1][0][1]=-1;
diter_toy[13][0][2][1][0][1]=-1;
diter_toy[14][0][2][1][0][1]=-1;
diter_toy[15][0][2][1][0][1]=-1;
//pTlead=7
//Bayes
diter_toy[2][0][3][0][0][1]=3; //chi2bf=2.9, chi2ch=17.5
diter_toy[4][0][3][0][0][1]=3; //chi2bf=0.5, chi2ch=11.2
diter_toy[5][0][3][0][0][1]=3; //chi2bf=1.3, chi2ch=8.6
diter_toy[6][0][3][0][0][1]=4; //chi2bf=2.2, chi2ch=9.5
diter_toy[7][0][3][0][0][1]=6; //chi2bf=2.0, chi2ch=13.5
diter_toy[8][0][3][0][0][1]=2; //chi2bf=2.0, chi2ch=16.7
diter_toy[9][0][3][0][0][1]=4; //chi2bf=2.9, chi2ch=17.3
diter_toy[10][0][3][0][0][1]=4; //chi2bf=1.9, chi2ch=14.6
diter_toy[11][0][3][0][0][1]=2; //chi2bf=2.6, chi2ch=13.0
diter_toy[12][0][3][0][0][1]=4; //chi2bf=1.7, chi2ch=10.3
diter_toy[13][0][3][0][0][1]=3; //chi2bf=2.2, chi2ch=19.3
diter_toy[14][0][3][0][0][1]=2; //chi2bf=3.1, chi2ch=12.1
diter_toy[15][0][3][0][0][1]=3; //chi2bf=2.3, chi2ch=13.8
//SVD
diter_toy[2][0][3][1][0][1]=3; //chi2bf=2.5, chi2ch=342.2
diter_toy[4][0][3][1][0][1]=3; //chi2bf=5.8, chi2ch=409.3
diter_toy[5][0][3][1][0][1]=3; //chi2bf=8.2, chi2ch=125.5
diter_toy[6][0][3][1][0][1]=3; //chi2bf=7.8, chi2ch=386.9
diter_toy[7][0][3][1][0][1]=6; //chi2bf=3.1, chi2ch=232.4
diter_toy[8][0][3][1][0][1]=3; //chi2bf=3.3, chi2ch=76.2
diter_toy[9][0][3][1][0][1]=4; //chi2bf=0.6, chi2ch=453.2
diter_toy[10][0][3][1][0][1]=4; //chi2bf=1.4, chi2ch=107.4
diter_toy[11][0][3][1][0][1]=3; //chi2bf=2.0, chi2ch=78.9
diter_toy[12][0][3][1][0][1]=4; //chi2bf=1.5, chi2ch=466.9
diter_toy[13][0][3][1][0][1]=3; //chi2bf=2.8, chi2ch=315.7
diter_toy[14][0][3][1][0][1]=3; //chi2bf=2.1, chi2ch=167.9
diter_toy[15][0][3][1][0][1]=3; //chi2bf=4.9, chi2ch=362.6
//R=0.4
//pTlead=5
//Bayes
diter_toy[2][2][2][0][0][1]=3; //chi2bf=2.4, chi2ch=18.6
diter_toy[4][2][2][0][0][1]=3; //chi2bf=7.5, chi2ch=19.0
diter_toy[5][2][2][0][0][1]=6; //chi2bf=14.7, chi2ch=16.3
diter_toy[6][2][2][0][0][1]=-1;
diter_toy[7][2][2][0][0][1]=-1;
diter_toy[8][2][2][0][0][1]=3; //chi2bf=1.6, chi2ch=11.7
diter_toy[9][2][2][0][0][1]=-1;
diter_toy[10][2][2][0][0][1]=-1;
diter_toy[11][2][2][0][0][1]=2; //chi2bf=5.5, chi2ch=6.8
diter_toy[12][2][2][0][0][1]=7; //chi2bf=9.3, chi2ch=13.8
diter_toy[13][2][2][0][0][1]=-1;
diter_toy[14][2][2][0][0][1]=3; //chi2bf=2.7, chi2ch=17.8
diter_toy[15][2][2][0][0][1]=2; //chi2bf=11.1, chi2ch=16.5
//SVD
diter_toy[2][2][2][1][0][1]=-1;
diter_toy[4][2][2][1][0][1]=4; //chi2bf=5.6, chi2ch=466.2
diter_toy[5][2][2][1][0][1]=-1;
diter_toy[6][2][2][1][0][1]=-1;
diter_toy[7][2][2][1][0][1]=-1;
diter_toy[8][2][2][1][0][1]=3; //chi2bf=4.7, chi2ch=150.2
diter_toy[9][2][2][1][0][1]=-1;
diter_toy[10][2][2][1][0][1]=-1;
diter_toy[11][2][2][1][0][1]=3; //chi2bf=3.3, chi2ch=379.6
diter_toy[12][2][2][1][0][1]=-1;
diter_toy[13][2][2][1][0][1]=-1;
diter_toy[14][2][2][1][0][1]=-1;
diter_toy[15][2][2][1][0][1]=-1;
//pTlead=7
//Bayes
diter_toy[2][2][3][0][0][1]=4; //chi2bf=4.1, chi2ch=16.4
diter_toy[4][2][3][0][0][1]=3; //chi2bf=2.0, chi2ch=8.2
diter_toy[5][2][3][0][0][1]=2; //chi2bf=6.7, chi2ch=13.7
diter_toy[6][2][3][0][0][1]=4; //chi2bf=5.7, chi2ch=14.6
diter_toy[7][2][3][0][0][1]=7; //chi2bf=2.3, chi2ch=13.8
diter_toy[8][2][3][0][0][1]=2; //chi2bf=2.0, chi2ch=3.7
diter_toy[9][2][3][0][0][1]=6; //chi2bf=8.4, chi2ch=15.0
diter_toy[10][2][3][0][0][1]=4; //chi2bf=1.2, chi2ch=9.9
diter_toy[11][2][3][0][0][1]=2; //chi2bf=8.5, chi2ch=17.6
diter_toy[12][2][3][0][0][1]=5; //chi2bf=7.5, chi2ch=18.3
diter_toy[13][2][3][0][0][1]=2; //chi2bf=1.4, chi2ch=5.7
diter_toy[14][2][3][0][0][1]=2; //chi2bf=6.7, chi2ch=12.0
diter_toy[15][2][3][0][0][1]=4; //chi2bf=5.5, chi2ch=13.6
//SVD
diter_toy[2][2][3][1][0][1]=-1;
diter_toy[4][2][3][1][0][1]=3; //chi2bf=5.5, chi2ch=185.0
diter_toy[5][2][3][1][0][1]=3; //chi2bf=6.3, chi2ch=184.9
diter_toy[6][2][3][1][0][1]=4; //chi2bf=1.5, chi2ch=245.7
diter_toy[7][2][3][1][0][1]=5; //chi2bf=2.1, chi2ch=453.6
diter_toy[8][2][3][1][0][1]=3; //chi2bf=2.1, chi2ch=15.4
diter_toy[9][2][3][1][0][1]=-1;
diter_toy[10][2][3][1][0][1]=3; //chi2bf=1.0, chi2ch=367.4
diter_toy[11][2][3][1][0][1]=3; //chi2bf=0.9, chi2ch=189.4
diter_toy[12][2][3][1][0][1]=-1;
diter_toy[13][2][3][1][0][1]=3; //chi2bf=1.1, chi2ch=92.3
diter_toy[14][2][3][1][0][1]=3; //chi2bf=0.7, chi2ch=302.9
diter_toy[15][2][3][1][0][1]=-1;
*/
