Double_t TsalisFitFunc(Double_t* x_val, Double_t* par)
{
   double A=par[0];
   double n=par[1];
   double T=par[2];
   double pT=x_val[0];
   
   double y=A*pT*TMath::Power((1.+pT/(n*T)),-n);
   return y;
} 
