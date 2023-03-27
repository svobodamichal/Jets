Double_t Eff_track_rec_function(Double_t* x,Double_t* par)
{   
    // Track reconstruction efficiency parametrization
    Double_t pt,y;
    Double_t A,B,C;A=par[0];B=par[1];C=par[2];
    pt=x[0];
    y=A*(exp(-pow(B/pt,C)));
    return y;
}

TF1* efficiency11(bool kcentral, TString ratio="AuAu")
{
   TF1* f_Efficiency = new TF1("f_EfficiencyCent",Eff_track_rec_function,0,10.0,3);
   if(kcentral)
   {
      if(ratio=="pp") //pp like ratio of hadrons (pi/K/p)
         f_Efficiency->SetParameters(7.45643e-01,1.43725e-01,2.02904e+00);
      else //AuAu like ratio of hadrons (pi/K/p)
         f_Efficiency->SetParameters(7.44207e-01,1.39140e-01,1.88680e+00);
   }
   else
   {
      if(ratio=="pp")
         f_Efficiency->SetParameters(9.06946e-01,1.45242e-01,2.87409e+00);
      else
         f_Efficiency->SetParameters(9.05783e-01,1.25129e-01,2.19900e+00);
   }
   //f_Efficiency->Draw();
return f_Efficiency;
}

/*
Double_t efficiency11_Stephen(Double_t pt, TF1* effLow, TF1* effHigh)
{
   Double_t eff;
   if(pt<=1.2)eff = effLow->Eval(pt);
   else eff = effHigh->Eval(pt);
   eff=eff*(24.0/22.0); //efficiency was calculated as N_reco(in 22 TPC sectors)/N_emb(24 sectors), but we need it as  N_reco(22)/N_emb(22)
   if(eff<=0)eff=0.0001;
return eff;
*/
