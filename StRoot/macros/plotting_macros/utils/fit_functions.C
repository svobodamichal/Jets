Double_t TsalisFitFunc(Double_t* x_val, Double_t* par)
{
   double A=par[0];
   double n=par[1];
   double T=par[2];
   double pT=x_val[0];
   
   double y=A*pT*TMath::Power((1.+pT/(n*T)),-n);
   return y;
}

Double_t hardjet(double *x, double *par)
{
  double pT = *x;

  Double_t C = 1/(1 + TMath::Exp(-(pT - 2.8)/0.3)); //cutoff function
Double_t J,y1,y2, B, T, n, m0, mu,pT0,R,a,b,A,pwr;
    B    = par[1];
    T    = par[2];
    n    = par[3];
    m0   = par[4];
    mu   = par[5];
    A=par[6];
    pwr=par[7];
    Double_t mT = TMath::Sqrt((pT-mu)*(pT-mu)+m0*m0);

   y2 = B/TMath::Power(1.0+(mT-m0)/(n*T),n);
   y1=A*TMath::Power(pT,-pwr);

   if(pT<10)
      J=y1;
   else if(pT<20)
   {
      Double_t c1=(pT-10.0)/10.0;
      J=((1-c1)*y1+c1*y2);
   }
    else J=y2;
if(J<0)J=0;

//pT-dependent RAA
Double_t RAA;
/*
Double_t RAA_lowpT=0.5;
if(par[8]>RAA_lowpT)RAA_lowpT=par[8];
if(pT<5)RAA=RAA_lowpT;
//else if(pT<10)RAA=1.0;
else if(pT<15)RAA=(RAA_lowpT-par[8])*(15-pT)/10+par[8];
else RAA=par[8];*/

RAA=par[8];
return par[0]*J*RAA/*C*/; //hard jet spectrum (with suppressed low pT part)
}


Double_t LevyFitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for d2N/(2pi*pT dpT dy)
    // taken from here: http://sampa.if.usp.br/~suaide/blog/files/papers/PLB6372006.pdf
    Double_t pT, y, B, T, n, m0, mu,pT0,R,a,b;
    B    = par[0];
    T    = par[1];
    n    = par[2];
    m0   = par[3];
    mu   = par[4];
    /*
    pT0=par[5];
    R=par[6];
    a=par[7];
    b=par[8];*/
    pT   = x_val[0];
    Double_t mT = TMath::Sqrt((pT-mu)*(pT-mu)+m0*m0);
    y = B/TMath::Power(1.0+(mT-m0)/(n*T),n);
    //y=y*(a+b*(1.0/(1.0+TMath::Exp((pT-pT0)/R))));
    return y;
} 
