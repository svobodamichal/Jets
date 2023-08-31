Double_t auauratio(Int_t particle, Double_t pt) {
  //this function returns abundancy of a given particle with transverse momentum pt in "hadrons"
  // particle labeling code: 1 piplus, 2 piminus, 3 kplus, 4 kminus, 5 proton, 6 antiproton
  // in case of protons parametrization of data is divided into a lower  and a higher pt interval, for pions and kaons one function seems over the studied range is OK
 
  Double_t p1l[7]={0,32117,30229.8,2336.8,2202.89,53.1456,44.743};
  Double_t p2l[7]={0,1.07999,1.07965,1.07435,1.07334,0.989034,0.99638};
  Double_t p3l[7]={0,0.0963853,0.0973201,0.126937,0.127424,0.286439,0.274714};

  Double_t p1h[7]={0,0,0,0,0,311.627,227.458};
  Double_t p2h[7]={0,0,0,0,0,3.77473,3.64532};
  Double_t p3h[7]={0,0,0,0,0,-17.7066,-17.2519};

  Double_t fpart=0.0;
  Double_t fhadr=0.0;
  for(int k=1;k<7;k++) {
    if(pt<=2.0) {
      fhadr+=p1l[k]*pt*pt*TMath::Power(1+(p2l[k]-1)*pt/p3l[k],-p2l[k]/(p2l[k]-1.0));
    } else {
      if(k<5) {
	fhadr+=p1l[k]*pt*pt*TMath::Power(1+(p2l[k]-1)*pt/p3l[k],-p2l[k]/(p2l[k]-1.0));
      } else {
	fhadr+=p1h[k]*TMath::Power(1.0+pt/p2h[k],p3h[k]);
      }
    }
  }
  
  if(pt<=2.0) {
    fpart=p1l[particle]*pt*pt*TMath::Power(1+(p2l[particle]-1)*pt/p3l[particle],-p2l[particle]/(p2l[particle]-1.0))/fhadr;
  } else {
    if(particle<5) {
      fpart=p1l[particle]*pt*pt*TMath::Power(1+(p2l[particle]-1)*pt/p3l[particle],-p2l[particle]/(p2l[particle]-1.0))/fhadr;
    } else {
      fpart=p1h[particle]*TMath::Power(1.0+pt/p2h[particle],p3h[particle])/fhadr;
    }
  }
  // cout << pt << "   " << fpart << endl;
  return fpart; 
}
  
    
    

 
    