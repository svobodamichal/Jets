Double_t ppratio(Int_t particle, Double_t pt) {
  //this function returns abundancy of a given particle with transverse momentum pt in "hadrons"
  // particle labeling code: 1 piplus, 2 piminus, 3 kplus, 4 kminus, 5 proton, 6 antiproton
  //for pt<=1.2 we use low pt parametrization, above the high pt one
 
  Double_t p1l[7]={0,11.2496,13.3369,0.423814,0.397589,0.231914,0.191385};
  Double_t p2l[7]={0,1.48493,1.19841,3.09346,4.00881,14.526, 21.2117};
  Double_t p3l[7]={0,-11.4475,-10.0152,-13.762,-17.0636,-51.1358,-74.3381};

  Double_t p1h[7]={0,11.4916,10.7937,3.77487,2.56875,0.978768,0.533545};
  Double_t p2h[7]={0,1.22731,1.27492,1.27857,1.60754,2.31051,2.99409};
  Double_t p3h[7]={0,-10.0128,-10.1654,-10.1423,-11.3835,-13.1034,-15.161};

  Double_t fpart=0.0;
  Double_t fhadr=0.0;
  for(int k=1;k<7;k++) {
    if(pt<=1.2) {
      fhadr+=p1l[k]*TMath::Power(1.0+pt/p2l[k],p3l[k]);
    } else {
      fhadr+=p1h[k]*TMath::Power(1.0+pt/p2h[k],p3h[k]);
    }
  }

  if(pt<=1.2) {
    //for particles with pt<=1.2 we use low pt parametrization, above the high pt parametrization
    fpart=p1l[particle]*TMath::Power(1.0+pt/p2l[particle],p3l[particle])/fhadr;
  } else {
    fpart=p1h[particle]*TMath::Power(1.0+pt/p2h[particle],p3h[particle])/fhadr;
  }
  
  // cout << pt << "   " << fpart << endl;
  return fpart; 
}
  

