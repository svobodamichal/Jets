#include "common.h"


//cut values used to determine optimal unfolded solutions
//values were determined using toymodel simulation, using macro ../3_additional_studies/optIter/show_iter.C
float backfold_cut[C_ntests][C_nsystems][C_nbinnings][C_nunfoldings][C_nR][C_npTlead];
float change_cut[C_ntests][C_nsystems][C_nbinnings][C_nunfoldings][C_nR][C_npTlead];
float curv_cut[C_ntests][C_nsystems][C_nbinnings][C_nunfoldings][C_nR][C_npTlead];

//default values
float backfold_cut_default[C_nunfoldings]={0.3,0.3};
float change_cut_default[C_nunfoldings]={0.25,1};
float curv_cut_default[C_nunfoldings]={50,50};

//test:Camb
//system:central
//bin:1
//unfolding:Bayes
//backfolded vs measured
//R=0.2, distance cut=0.25
backfold_cut[1][0][1][0][0][5]=0.577826;
backfold_cut[1][0][1][0][0][6]=0.219348;
backfold_cut[1][0][1][0][0][7]=0.175784;
//R=0.3, distance cut=0.25
backfold_cut[1][0][1][0][1][5]=0.132524;
backfold_cut[1][0][1][0][1][6]=0.163216;
backfold_cut[1][0][1][0][1][7]=0.118776;
//R=0.4, distance cut=0.5
backfold_cut[1][0][1][0][2][5]=0.264208;
backfold_cut[1][0][1][0][2][6]=0.259395;
backfold_cut[1][0][1][0][2][7]=0.218236;
//successive iterations
//R=0.2, distance cut=0.25
change_cut[1][0][1][0][0][5]=0.296762;
change_cut[1][0][1][0][0][6]=0.413902;
change_cut[1][0][1][0][0][7]=0.238172;
//R=0.3, distance cut=0.25
change_cut[1][0][1][0][1][5]=0.0703048;
change_cut[1][0][1][0][1][6]=0.101614;
change_cut[1][0][1][0][1][7]=0.106275;
//R=0.4, distance cut=0.5
change_cut[1][0][1][0][2][5]=0.144111;
change_cut[1][0][1][0][2][6]=0.167343;
change_cut[1][0][1][0][2][7]=0.181964;

//unfolding:SVD
//backfolded vs measured
//R=0.2, distance cut=0.25
backfold_cut[1][0][1][1][0][5]=0.196961;
backfold_cut[1][0][1][1][0][6]=0.206371;
backfold_cut[1][0][1][1][0][7]=0.16953;
//R=0.3, distance cut=0.25
backfold_cut[1][0][1][1][1][5]=0.121077;
backfold_cut[1][0][1][1][1][6]=0.150532;
backfold_cut[1][0][1][1][1][7]=0.117649;
//R=0.4, distance cut=0.5
backfold_cut[1][0][1][1][2][5]=0.240535;
backfold_cut[1][0][1][1][2][6]=0.205029;
backfold_cut[1][0][1][1][2][7]=0.180703;
//successive iterations
//R=0.2, distance cut=0.25
change_cut[1][0][1][1][0][5]=1;
change_cut[1][0][1][1][0][6]=1;
change_cut[1][0][1][1][0][7]=1;
//R=0.3, distance cut=0.25
change_cut[1][0][1][1][1][5]=1;
change_cut[1][0][1][1][1][6]=1;
change_cut[1][0][1][1][1][7]=1;
//R=0.4, distance cut=0.5
change_cut[1][0][1][1][2][5]=1;
change_cut[1][0][1][1][2][6]=1;
change_cut[1][0][1][1][2][7]=1;



//bin:2
//unfolding:Bayes
//backfolded vs measured
//R=0.2, distance cut=0.25
backfold_cut[1][0][2][0][0][5]=0.272512;
backfold_cut[1][0][2][0][0][6]=0.26334;
backfold_cut[1][0][2][0][0][7]=0.202324;
//R=0.3, distance cut=0.25
backfold_cut[1][0][2][0][1][5]=0.0983119;
backfold_cut[1][0][2][0][1][6]=0.0913385;
backfold_cut[1][0][2][0][1][7]=0.14218;
//R=0.4, distance cut=0.5
backfold_cut[1][0][2][0][2][5]=0.288217;
backfold_cut[1][0][2][0][2][6]=0.2868;
backfold_cut[1][0][2][0][2][7]=0.264442;
//successive iterations
//R=0.2, distance cut=0.25
change_cut[1][0][2][0][0][5]=0.351693;
change_cut[1][0][2][0][0][6]=0.622305;
change_cut[1][0][2][0][0][7]=3.88456;
//R=0.3, distance cut=0.25
change_cut[1][0][2][0][1][5]=0.00519761;
change_cut[1][0][2][0][1][6]=0.0743515;
change_cut[1][0][2][0][1][7]=0.196739;
//R=0.4, distance cut=0.5
change_cut[1][0][2][0][2][5]=0.179546;
change_cut[1][0][2][0][2][6]=0.188311;
change_cut[1][0][2][0][2][7]=0.235979;

//unfolding:SVD
//backfolded vs measured
//R=0.2, distance cut=0.25
backfold_cut[1][0][2][1][0][5]=0.102912;
backfold_cut[1][0][2][1][0][6]=0.144537;
backfold_cut[1][0][2][1][0][7]=0.143509;
//R=0.3, distance cut=0.25
backfold_cut[1][0][2][1][1][5]=0.0844647;
backfold_cut[1][0][2][1][1][6]=0.100416;
backfold_cut[1][0][2][1][1][7]=0.0992652;
//R=0.4, distance cut=0.5
backfold_cut[1][0][2][1][2][5]=0.20598;
backfold_cut[1][0][2][1][2][6]=0.284592;
backfold_cut[1][0][2][1][2][7]=0.183484;
//successive iterations
//R=0.2, distance cut=0.25
change_cut[1][0][2][1][0][5]=1;
change_cut[1][0][2][1][0][6]=1;
change_cut[1][0][2][1][0][7]=1;
//R=0.3, distance cut=0.25
change_cut[1][0][2][1][1][5]=1;
change_cut[1][0][2][1][1][6]=1;
change_cut[1][0][2][1][1][7]=1;
//R=0.4, distance cut=0.5
change_cut[1][0][2][1][2][5]=1;
change_cut[1][0][2][1][2][6]=1;
change_cut[1][0][2][1][2][7]=1;


//bin:3
//unfolding:Bayes
//backfolded vs measured
//R=0.2, distance cut=0.25
backfold_cut[1][0][3][0][0][5]=0.2;
backfold_cut[1][0][3][0][0][6]=0.2;
backfold_cut[1][0][3][0][0][7]=0.2;
//R=0.3, distance cut=0.25
backfold_cut[1][0][3][0][1][5]=0.261301;
backfold_cut[1][0][3][0][1][6]=0.198134;
backfold_cut[1][0][3][0][1][7]=0.133765;
//R=0.4, distance cut=0.5
backfold_cut[1][0][3][0][2][5]=0.843968;
backfold_cut[1][0][3][0][2][6]=0.457536;
backfold_cut[1][0][3][0][2][7]=0.274018;
//successive iterations
//R=0.2, distance cut=0.25
change_cut[1][0][3][0][0][5]=0.20963;
change_cut[1][0][3][0][0][6]=0.204724;
change_cut[1][0][3][0][0][7]=0.247308;
//R=0.3, distance cut=0.25
change_cut[1][0][3][0][1][5]=0.123322;
change_cut[1][0][3][0][1][6]=0.129661;
change_cut[1][0][3][0][1][7]=0.083252;
//R=0.4, distance cut=0.5
change_cut[1][0][3][0][2][5]=0.187705;
change_cut[1][0][3][0][2][6]=0.188731;
change_cut[1][0][3][0][2][7]=0.176654;

//unfolding:SVD
//backfolded vs measured
//R=0.2, distance cut=0.25
backfold_cut[1][0][3][1][0][5]=0.11327;
backfold_cut[1][0][3][1][0][6]=0.271654;
backfold_cut[1][0][3][1][0][7]=0.207152;
//R=0.3, distance cut=0.25
backfold_cut[1][0][3][1][1][5]=0.1546;
backfold_cut[1][0][3][1][1][6]=0.187535;
backfold_cut[1][0][3][1][1][7]=0.160961;
//R=0.4, distance cut=0.5
backfold_cut[1][0][3][1][2][5]=0.668753;
backfold_cut[1][0][3][1][2][6]=0.669283;
backfold_cut[1][0][3][1][2][7]=0.490263;
//successive iterations
//R=0.2, distance cut=0.25
change_cut[1][0][3][1][0][5]=1;
change_cut[1][0][3][1][0][6]=1;
change_cut[1][0][3][1][0][7]=1;
//R=0.3, distance cut=0.25
change_cut[1][0][3][1][1][5]=1;
change_cut[1][0][3][1][1][6]=1;
change_cut[1][0][3][1][1][7]=1;
//R=0.4, distance cut=0.5
change_cut[1][0][3][1][2][5]=1;
change_cut[1][0][3][1][2][6]=1;
change_cut[1][0][3][1][2][7]=1;


//bin:4
//unfolding:Bayes
//backfolded vs measured
//R=0.2, distance cut=0.25
backfold_cut[1][0][4][0][0][5]=0.381628;
backfold_cut[1][0][4][0][0][6]=0.373785;
backfold_cut[1][0][4][0][0][7]=0.33431;
//R=0.3, distance cut=0.25
backfold_cut[1][0][4][0][1][5]=0.0985076;
backfold_cut[1][0][4][0][1][6]=0.136801;
backfold_cut[1][0][4][0][1][7]=0.182452;
//R=0.4, distance cut=0.5
backfold_cut[1][0][4][0][2][5]=0.522688;
backfold_cut[1][0][4][0][2][6]=0.412902;
backfold_cut[1][0][4][0][2][7]=0.316811;
//successive iterations
//R=0.2, distance cut=0.25
change_cut[1][0][4][0][0][5]=0.36095;
change_cut[1][0][4][0][0][6]=0.432813;
change_cut[1][0][4][0][0][7]=0.257282;
//R=0.3, distance cut=0.25
change_cut[1][0][4][0][1][5]=0.0301113;
change_cut[1][0][4][0][1][6]=0.0546406;
change_cut[1][0][4][0][1][7]=0.138043;
//R=0.4, distance cut=0.5
change_cut[1][0][4][0][2][5]=0.145249;
change_cut[1][0][4][0][2][6]=0.140149;
change_cut[1][0][4][0][2][7]=0.186433;

//unfolding:SVD
//backfolded vs measured
//R=0.2, distance cut=0.25
backfold_cut[1][0][4][1][0][5]=0.308717;
backfold_cut[1][0][4][1][0][6]=0.312946;
backfold_cut[1][0][4][1][0][7]=0.244476;
//R=0.3, distance cut=0.25
backfold_cut[1][0][4][1][1][5]=0.128603;
backfold_cut[1][0][4][1][1][6]=0.160723;
backfold_cut[1][0][4][1][1][7]=0.149356;
//R=0.4, distance cut=0.5
backfold_cut[1][0][4][1][2][5]=0.374931;
backfold_cut[1][0][4][1][2][6]=0.350433;
backfold_cut[1][0][4][1][2][7]=0.169147;
//successive iterations
//R=0.2, distance cut=0.25
change_cut[1][0][4][1][0][5]=1;
change_cut[1][0][4][1][0][6]=1;
change_cut[1][0][4][1][0][7]=1;
//R=0.3, distance cut=0.25
change_cut[1][0][4][1][1][5]=1;
change_cut[1][0][4][1][1][6]=1;
change_cut[1][0][4][1][1][7]=1;
//R=0.4, distance cut=0.5
change_cut[1][0][4][1][2][5]=1;
change_cut[1][0][4][1][2][6]=1;
change_cut[1][0][4][1][2][7]=1;
