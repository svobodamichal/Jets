#include <TSystem.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TRandom.h>
#include <TMath.h>

#include <iostream>

using namespace std;

/*
Double_t levy(Double_t *x, Double_t *par)
{
	Double_t pT=x[0];
	Double_t A=par[0];
	Double_t mu=par[1];
	Double_t c=par[2];
	Double_t pow=par[3];
	Double_t x0=par[4];
	Double_t y=(TMath::Exp(-(c/(2*(pT-mu))))/TMath::Power((pT-mu),pow));
   y = A*(1.0/(1.0+TMath::Exp(-(pT-x0)/0.8)))*y; //suppress low-pT part
	return y;
}*/

Double_t LevyFitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for d2N/(2pi*pT dpT dy)
    // taken from here: http://sampa.if.usp.br/~suaide/blog/files/papers/PLB6372006.pdf
    Double_t pT, y, A, T, n, m0, mu, x0;
    A    = par[0];
    T    = par[1];
    n    = par[2];
    m0   = par[3];
    mu   = par[4];
	 x0   = par[5];
    pT   = x_val[0];
    Double_t mT = TMath::Sqrt((pT-mu)*(pT-mu)+m0*m0);
    y = 1/TMath::Power(1.0+(mT-m0)/(n*T),n);
    y =(1.0/(1.0+TMath::Exp(-(pT-x0)/1.0)))*y; //suppress low-pT part

    return A*y;
}

Double_t powlaw(Double_t* x_val, Double_t* par)
{
	Double_t y, A, pT, pow, x0;
	pT=x_val[0];
	A=par[0];
	pow=par[1];
   x0=par[2];
	y=(1/TMath::Power(pT,pow));
	y=y*(1.0/(1.0+TMath::Exp(-(pT-x0)/(1.1-pow/10.)))); //suppress low-pT part
   return A*y;
}

Double_t flat(Double_t* x_val, Double_t* par)
{ 
	Double_t y=par[0];
	return y;
}

Double_t TsalisFitFunc(Double_t* x_val, Double_t* par)
{
   double A=par[0];
   double n=par[1];
   double T=par[2];
   double pT=x_val[0];
   
   double y=A*pT*TMath::Power((1.+pT/(n*T)),-n);
   return y;
}
 /*	
Double_t PythiaFitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for d2N/(2pi*pT dpT dy)
    // taken from here: http://sampa.if.usp.br/~suaide/blog/files/papers/PLB6372006.pdf
    Double_t pT, y,y1,y2, B1, T1, n1, m1, mu1,B2, T2, n2, m2, mu2;
    B1    = par[0];
    T1    = par[1];
    n1    = par[2];
    m1   = par[3];
	 mu1  = par[4];
	 B2    = par[5];
    T2    = par[6];
    n2    = par[7];
    m2   = par[8];
	 mu2   = par[9];

	 pT=x_val[0];
	 
    Double_t mT1 = TMath::Sqrt((pT-mu1)*(pT-mu1)+m1*m1);
	 Double_t mT2 = TMath::Sqrt((pT-mu2)*(pT-mu2)+m2*m2);
	 
	y1 =B1/TMath::Power(1.0+(mT1-m1)/(n1*T1),n1);
	y2=B2/TMath::Power(1.0+(mT2-m2)/(n2*T2),n2);
	
	
	
	if(pT<15)
		y=y1;
	else if(pT<25)
	{
		Double_t c1=(pT-15.0)/10.0;
		y=((1-c1)*y1+c1*y2);
	}
	 else y=y2;

	 return y;
}
*/


//___________________________________________________________________
prior_maker(const TString output="./")
{
	static const Int_t npriors_py=3;
	static const int nr=4;
	static const Int_t nLeadCuts=8;
	static const Int_t nptbins=800;
	static const Float_t ptminbin=-100;
	static const Float_t ptmaxbin=100;
	static const Int_t npows=4;
	static const Int_t firstpow=3;


	Int_t lastpow;
	Float_t pTleadcuts[nLeadCuts];
	Float_t R[nr];
 	TFile *fpyt[nr];
 	TFile *foutput;
	TH1D* hprior_py[nr][nLeadCuts];
	TH1D* hprior_pyde[nr][nLeadCuts];
	TH1D* hdiv_py[nr][nLeadCuts];
	TH1D* hdiv_pyde[nr][nLeadCuts];

   TH1D *hprior_pythia[nr][nLeadCuts];
   TH1D *hprior_pythiadete[nr][nLeadCuts];
	TH1D *hprior_powlaw[npows][nLeadCuts];
	TH1D *hprior_flat[nLeadCuts];
   TH1D *hprior_gauss[nLeadCuts];
   TH1D *hprior_levy[nLeadCuts];
   TH1D *hprior_levy_alex[nLeadCuts];

	
	TString str;

	R[0]=0.2;
	R[1]=0.3;
	R[2]=0.4;
	R[3]=0.5;
	pTleadcuts[0]=0;
	pTleadcuts[1]=1;
	pTleadcuts[2]=2;
	pTleadcuts[3]=3;
	pTleadcuts[4]=4;
	pTleadcuts[5]=5;
	pTleadcuts[6]=6;
	pTleadcuts[7]=7;

  str = Form("%s/priors.root", output.Data());
  foutput = new TFile(str.Data(), "RECREATE");

  TH2::SetDefaultSumw2();
  
    Double_t pTmax=ptmaxbin;

//Various parametrizations of Tsalis function
const int ntsalis=9; //number of parametrizations
float n_ts[ntsalis]={4,8,12,8,12,16,12,16,20}; //n
float T_ts[ntsalis]={0.6,0.6,0.6,0.9,0.9,0.9,1.2,1.2,1.2};//T

//charged Pythia jets fit parameter values

float A_ts_pyt[nr][nLeadCuts]={{2.83360e+01,2.83360e+01,2.83360e+01,8.81968e-03,4.02710e-04,1.43305e-04,5.42625e-05,2.18447e-05},
{2.53645e+01,2.53645e+01,2.53645e+01,2.00188e-02,6.41934e-04,1.60507e-04,4.93945e-05,1.72143e-05},
{1.12320e+01,1.12320e+01,1.12320e+01,1.01828e-02,8.80889e-04,1.93556e-04,5.61829e-05,1.89797e-05},
{5.53639e+00,5.53639e+00,5.53639e+00,8.77900e-03,9.89249e-04,2.04819e-04,5.62038e-05,1.79458e-05}};
float n_ts_pyt[nr][nLeadCuts]={{9.47989e+00,9.47989e+00,9.47989e+00,1.46500e+01,2.10481e+01,2.53750e+01,3.18281e+01,4.23864e+01},
{9.70878e+00,9.70878e+00,9.70878e+00,1.45632e+01,2.10545e+01,2.68948e+01,3.61376e+01,5.43851e+01},
{9.96063e+00,9.96063e+00,9.96063e+00,1.60387e+01,2.12647e+01,2.77678e+01,3.82714e+01,5.99171e+01},
{1.06730e+01,1.06730e+01,1.06730e+01,1.68246e+01,2.16281e+01,2.84937e+01,4.02036e+01,6.74942e+01}};
float T_ts_pyt[nr][nLeadCuts]={{1.87938e-01,1.87938e-01,1.87938e-01,8.01893e-01,1.33712e+00,1.57866e+00,1.83741e+00,2.11036e+00},
{2.08829e-01,2.08829e-01,2.08829e-01,7.63354e-01,1.32177e+00,1.64012e+00,1.96193e+00,2.29803e+00},
{2.51388e-01,2.51388e-01,2.51388e-01,9.10532e-01,1.33535e+00,1.68376e+00,2.02680e+00,2.37875e+00},
{3.16090e-01,3.16090e-01,3.16090e-01,9.72541e-01,1.35553e+00,1.71752e+00,2.08017e+00,2.46061e+00}};

TF1* ftsalis = new TF1("tasalis",TsalisFitFunc,0.2,pTmax,3);
		ftsalis->SetNpx(10000);
		ftsalis->SetParameter(0,1); // A, changes the amplitude, 1

TF1* levy1 = new TF1("levy1",LevyFitFunc,0.2,pTmax,6);
	  levy1->SetNpx(10000);
 	  levy1->SetParameter(0,1); // A, changes the amplitude, 0.1
	  levy1->SetParameter(3,0.0001); // m0, changes the width, 0.0001
	  levy1->SetParameter(4,0.0); // mu, changes the x-axis shift, 0.0
	  levy1->SetParameter(5,4.0); // x0, position of the low-pT knee
	  levy1->SetParameter(1,0.4); // T, changes slope, 0.4
	  levy1->SetParameter(2,5.8); // n, changes how fast spectrum drops, 5.8

TF1 *fplaw = new TF1("fplaw",powlaw,1,pTmax,3);
     fplaw->SetParameters(1.0,5,4);
     fplaw->SetParNames("amplitude","power","x0");
	  fplaw->SetNpx(10000);

TF1 *fflat = new TF1("fflat",flat,0,pTmax,1);
		fflat->SetParameter(0,1.0);

TF1 *fpythia = new TF1("fpythia",TsalisFitFunc,0.2,pTmax,3);
		fpythia->SetNpx(10000);

//TF1 *fpythia = new TF1("fpythia",PythiaFitFunc,0.2,pTmax,10);

//PYTHIA jet distribution parameters - obtained from fit to charged jet PYTHIA spectrum for R=0.2/0.3/0.4 and pTleading>5GeV/c
/*
Double_t B1[5][3]= //[pTleading cut: 3,4,5,7][R: 0.2,0.3,0.4]
{
	{2.84774e+00,2.38659e-01,2.83075e-01}, //pTlead>3 instead of 0 (I don't have a good fit for pTlead>0)
	{2.84774e+00,2.38659e-01,2.83075e-01}, //pTlead>3
	{1.73173e-01,1.41869e-01,2.17150e-01}, //pTlead>4
	{3.57068e-02,2.61248e-01,1.39465e+00}, //pTlead>5
	{7.48557e-02,3.59904e+00,7.22683e+00}  //pTlead>7
};

Double_t T1[5][3]= //[pTleading cut: 0,3,4,5,7][R: 0.2,0.3,0.4]
{
	{2.49993e-01,2.76819e-01,3.00461e-01},
	{2.49993e-01,2.76819e-01,3.00461e-01},
	{3.28990e-01,6.37244e-01,7.60106e-01},
	{7.43245e-01,7.87954e-01,1.02352e+00},
	{1.18772e+00,1.45543e+00,1.60960e+00}
};

Double_t n1[5][3]= //[pTleading cut: 0,3,4,5,7][R: 0.2,0.3,0.4]
{
	{6.76771e+00,6.14446e+00,6.92238e+00},
	{6.76771e+00,6.14446e+00,6.92238e+00},
	{6.23665e+00,8.87298e+00,1.06110e+01},
	{9.25722e+00,1.05483e+01,1.47328e+01},
	{1.34568e+01,2.26235e+01,2.57591e+01}
};

Double_t m01[5][3]= //[pTleading cut: 0,3,4,5,7][R: 0.2,0.3,0.4]
{
	{-1.05394e+00,9.83217e+00,1.62457e+01},
	{-1.05394e+00,9.83217e+00,1.62457e+01},
	{3.94223e+00,2.35110e+00,5.25541e+00},
	{4.45090e+00,2.78336e+00,-1.48311e+00},
	{-1.25106e+00,-8.02311e+00,-1.02241e+01}
};

Double_t mu1[5][3]= //[pTleading cut: 0,3,4,5,7][R: 0.2,0.3,0.4]
{
	{1.59985e+00,-3.16050e+00,-5.41446e+00},
	{1.59985e+00,-3.16050e+00,-5.41446e+00},
	{-1.63418e+00,-2.87416e+00,-6.28665e+00},
	{-4.55443e+00,-6.50919e+00,-7.17847e+00},
	{-6.78114e+00,-8.30787e+00,-8.73165e+00}
};

Double_t B2[5][3]= //[pTleading cut: 0,3,4,5,7][R: 0.2,0.3,0.4]
{
	{1.62851e+01,1.56709e+01,1.15048e+01},
	{1.62851e+01,1.56709e+01,1.15048e+01},
	{1.74172e+01,2.24953e+01,1.80006e+01},
	{2.04519e+01,1.73324e+01,2.81221e+01},
	{2.15284e+01,3.28406e+01,2.97309e+01}
};

Double_t T2[5][3]= //[pTleading cut: 0,3,4,5,7][R: 0.2,0.3,0.4]
{
	{2.35837e+00,2.21367e+00,2.28350e+00},
	{2.35837e+00,2.21367e+00,2.28350e+00},
	{2.27170e+00,2.47576e+00,2.48260e+00},
	{2.44757e+00,2.35363e+00,2.65813e+00},
	{2.41168e+00,2.64063e+00,2.61912e+00}
};

Double_t n2[5][3]= //[pTleading cut: 0,3,4,5,7][R: 0.2,0.3,0.4]
{
	{9.88238e+01,7.65644e+01,7.15251e+01},
	{9.88238e+01,7.65644e+01,7.15251e+01},
	{9.62992e+01,1.06243e+02,9.38929e+01},
	{1.12497e+02,9.04842e+01,1.30241e+02},
	{1.16419e+02,1.51931e+02,1.47563e+02}
};

Double_t m02[5][3]= //[pTleading cut: 0,3,4,5,7][R: 0.2,0.3,0.4]
{
	{-1.58889e+01,-1.49943e+01,1.20716e+01},
	{-1.58889e+01,-1.49943e+01,1.20716e+01},
	{-1.78200e+01,-1.27342e+01,-5.71891e+00},
	{-1.51199e+01,-1.42919e+01,-6.10978e+00},
	{-1.76363e+01,-1.52363e+01,-1.65992e+01}
};

Double_t mu2[5][3]= //[pTleading cut: 0,3,4,5,7][R: 0.2,0.3,0.4]
{
	{-1.01711e+01,-8.19611e+00,-3.76573e+01},
	{-1.01711e+01,-8.19611e+00,-3.76573e+01},
	{-5.21055e+00,-1.67921e+01,-2.47554e+01},
	{-1.34549e+01,-1.18774e+01,-2.79141e+01},
	{-8.76751e+00,-1.68410e+01,-1.34104e+01}
};

*/


foutput->cd();

fflat->Write("flat");

levy1->Write("levy_1");
/*
levy1->SetParameter(1,0.4); // T, changes slope, 0.4
levy1->SetParameter(2,4.4); // n, changes how fast spectrum drops, 5.8
levy1->Write("levy_2");
*/

for(int i=0; i<ntsalis;i++)
{
	ftsalis->SetParameter(1,n_ts[i]);
	ftsalis->SetParameter(2,T_ts[i]);
	ftsalis->Write(Form("tsalis_%i",i+1));
}

for(int power=3; power<7; power++){
	fplaw->SetParameter(1,power);
	fplaw->Write(Form("powlaw%i",power));
}
//pT^4.5
float fpwr=4.5;
fplaw->SetParameter(1,fpwr);
fplaw->Write(Form("powlaw%.0lf",fpwr*10));
//pT^5.5
fpwr=5.5;
fplaw->SetParameter(1,fpwr);
fplaw->Write(Form("powlaw%.0lf",fpwr*10));
//pT^5.5


//fit PYTHIA distributions
for(int pTlidx=0; pTlidx<nLeadCuts;pTlidx++){
for(int ridx=0;ridx<nr;ridx++){ 
	float Rpar=R[ridx];
   float pTlead=pTleadcuts[pTlidx];
	cout<<"PYTHIA: R="<<Rpar<<endl;
   fpythia->SetParameters(A_ts_pyt[ridx][pTlidx],n_ts_pyt[ridx][pTlidx],T_ts_pyt[ridx][pTlidx]);
   fpythia->Write(Form("pythia_R%.0lf_pTlead%.0lf",Rpar*10,pTlead));

}//R loop
}//pTleading cut loop
/*
foutput->cd();
foutput->Write();*/
foutput->Close();
delete foutput;

}
