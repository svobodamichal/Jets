Double_t ppratio(Int_t PID, Double_t pt) {
  //this function returns abundancy of a given particle with transverse momentum pt in "hadrons"
	 // PID: p, anti-p, pi+, pi-, K+, K-
  //for pt<=1.2 we use low pt parametrization, above the high pt one
 
  Int_t SwitchArray[6]={5,6,1,2,3,4}; //change PID in turn
  // old particle labeling code: 1 piplus, 2 piminus, 3 kplus, 4 kminus, 5 proton, 6 antiproton
  Int_t particle=SwitchArray[PID];
  
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
Double_t auauratio(Int_t PID, Double_t pt) {
   //this function returns abundancy of a given particle with transverse momentum pt in "hadrons"
   // PID: p, anti-p, pi+, pi-, K+, K-
   // in case of protons parametrization of data is divided into a lower  and a higher pt interval, for pions and kaons one function seems over the studied range is OK
 
   Int_t SwitchArray[6]={5,6,1,2,3,4}; //change PID in turn
  // old particle labeling code: 1 piplus, 2 piminus, 3 kplus, 4 kminus, 5 proton, 6 antiproton
  Int_t particle=SwitchArray[PID];
  
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

//=========================================================
Double_t Eff_track_rec_function(Double_t* x,Double_t* par)
{
    // Track reconstruction efficiency parametrization
    Double_t pt,y;
    Double_t A,B,C;A=par[0];B=par[1];C=par[2];
    pt=x[0];
    y=A*(exp(-pow(B/pt,C)));
    return y;
}

//=========================================================

Double_t eff_rat_fit_function(Double_t* x,Double_t* par)
{
    // Track reconstruction efficiency parametrization
    Double_t pt,y;
    Double_t A,B,C,D; A=par[0];B=par[1];C=par[2];D=par[3];E=par[4];
    pt=x[0];
	 if(pt<1)
    y=A*(exp(-B*pt+C))+D*pt+E;
	 else
		 y=D*pt+E;
    return y;
}

//=========================================================
void make_histo(bool kCentrality=0/*0-central,1-peripheral*/,TString ratio="pp"/*particle ratio: pp or AuAu*/,Int_t ncount=1E3)
{
	TString fname[2]={"central","peripheral"};
	const Float_t pTmax=5.2;
	const Float_t pTmin=0.2;
	const Int_t nbins=500;
	
	 Double_t Global_ABC_Parameters_A[9][7][6];  //centrality(9),runID(7),PID(6)
    Double_t Global_ABC_Parameters_B[9][7][6];  //centrality(9),runID(7),PID(6)
    Double_t Global_ABC_Parameters_C[9][7][6];  //centrality(9),runID(7),PID(6)
	 Int_t Switch_Array[6]={2,3,4,5,0,1}; //change PID in turn
	 
    TString eEff_file="EffParameters2011_New.root";
    cout << "Read track reconstruction efficiency parameter file: " << eEff_file.Data() << endl;
    TFile *Eff_file=new TFile(eEff_file.Data());
    TH3D* hA  =(TH3D*)Eff_file->Get("hA");
    TH3D* hB  =(TH3D*)Eff_file->Get("hB");
    TH3D* hC  =(TH3D*)Eff_file->Get("hC");
    TH2D* hsA =(TH2D*)Eff_file->Get("hsA");
    TH2D* hsB =(TH2D*)Eff_file->Get("hsB");
    TH2D* hsC =(TH2D*)Eff_file->Get("hsC");

    for(Int_t i = 0; i < 6; i++) // PID: p, anti-p, pi+, pi-, K+, K-
    {
        for(Int_t j = 0; j < 9; j++) // Centrality
        {
            for(Int_t k = 0; k < 7; k++) // runID
            {
                Global_ABC_Parameters_A[j][k][Switch_Array[i]] = hA->GetBinContent(hA->FindBin(i,j,k));
                Global_ABC_Parameters_B[j][k][Switch_Array[i]] = hB->GetBinContent(hB->FindBin(i,j,k));
                Global_ABC_Parameters_C[j][k][Switch_Array[i]] = hC->GetBinContent(hC->FindBin(i,j,k));
            }
        }
    }
	TString outname=Form("./effi_%s_%s.root",fname[kCentrality].Data(),ratio.Data());
	cout<<"outname: "<<outname.Data()<<endl;
   TFile* fout= new TFile(outname.Data(),"RECREATE");
   fout->cd();
	TF1* f_EfficiencyVsPt = new TF1("f_EfficiencyVsPt",Eff_track_rec_function,0,50.0,3);
   TH1D* heffi=new TH1D("heffi","heffi",nbins,0,pTmax);	 
	TH1D* hpart=new TH1D("hpart","hpart",6,-0.5,5.5);

	for(int i=0; i<ncount; i++)
	{
		if(i%1000==0)cout<<"running event "<<i<<endl;
		Double_t pT=gRandom->Uniform(pTmin,pTmax);
		Double_t rnd=gRandom->Uniform(0,1);
		Double_t sum=0;
		Double_t particle=0;
			for(int part=0; part<6; part++)
			{
			if(rnd>sum)particle=part;
			if(ratio=="pp")sum+=ppratio(part, pT);
         else sum+=auauratio(part, pT);
			//cout<<"sum:"<<sum<<" rnd:"<<rnd<<endl;
			}
		//cout<<"pT: "<<pT<<" particle: "<<particle<<endl;
		hpart->Fill(particle);
	
	
	
	 // Centrality: 0-5, 5-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80
    Int_t cent9_eff = 8;
	 Double_t rnd2=gRandom->Uniform(0,1);
    if(kCentrality == 0) // 0-10% -> randomly choose between efficiency for 0-5% and 5%-10%
    {
        cent9_eff = 0;
        if(rnd2 < 0.5) cent9_eff = 1;
    }
    if(kCentrality == 1) // 60-80% -> randomly choose between efficiency for 60%-70% and 70%-80%
    {
        cent9_eff = 7;
        if(rnd2 < 0.5) cent9_eff = 8;
    }
    
    Int_t RunId     = 0; //all runIds together
	
	//cout << "Initialize track reconstruction efficiency functions" << endl;
    f_EfficiencyVsPt->SetParameter(0,Global_ABC_Parameters_A[cent9_eff][RunId][particle]);
    f_EfficiencyVsPt->SetParameter(1,Global_ABC_Parameters_B[cent9_eff][RunId][particle]);
    f_EfficiencyVsPt->SetParameter(2,Global_ABC_Parameters_C[cent9_eff][RunId][particle]);
    
	
	//TH1D* heffi=new TH1D*("heffi","heffi",(Int_t) 10*pTmax,0,pTmax);

    Double_t epsilon = f_EfficiencyVsPt->Eval(pT);
	 heffi->Fill(pT,epsilon*nbins/ncount);
    }//event loop
    heffi->Write();
	 hpart->Write();
	 fout->Close();
	 delete fout;
}
//==========================================================

void make_ratio(TString ratio="pp")
{
	TString eEff_file_per=Form("effi_peripheral_%s.root",ratio.Data());
	TString eEff_file_cent=Form("effi_central_%s.root",ratio.Data());
	TFile *Eff_file_per=new TFile(eEff_file_per.Data(),"OPEN");
	TFile *Eff_file_cent=new TFile(eEff_file_cent.Data(),"OPEN");
	
	TH1D* heff_per=(TH1D*)Eff_file_per->Get("heffi");
	TH1D* heff_cent=(TH1D*)Eff_file_cent->Get("heffi");
	
	heff_per->Divide(heff_cent);
	
	TF1* ffit=new TF1("ffit",eff_rat_fit_function,0.2,5.0,5);
	heff_per->Fit(ffit,"R");
	//heff_per->Fit(ffit,"R");
}

void fit_efficiency(TString ratio="pp")
{

	TF1* f_EfficiencyCent = new TF1("f_EfficiencyCent",Eff_track_rec_function,0,10.0,3);
	TF1* f_EfficiencyPer = new TF1("f_EfficiencyPer",Eff_track_rec_function,0,10.0,3);
	f_EfficiencyCent->SetParameters(1,1,1);
	f_EfficiencyPer->SetParameters(1,1,1);
	
   TString eEff_file_per=Form("effi_peripheral_%s.root",ratio.Data());
   TString eEff_file_cent=Form("effi_central_%s.root",ratio.Data());
	TFile *Eff_file_per=new TFile(eEff_file_per.Data(),"OPEN");
	TFile *Eff_file_cent=new TFile(eEff_file_cent.Data(),"OPEN");
	
	TH1D* heff_per=(TH1D*)Eff_file_per->Get("heffi");
	TH1D* heff_cent=(TH1D*)Eff_file_cent->Get("heffi");
	
	cout<<"fitting central..."<<endl;
	heff_cent->Fit(f_EfficiencyCent);
	cout<<"fitting peripheral..."<<endl;
	heff_per->Fit(f_EfficiencyPer);
	
	heff_per->Draw("");
   heff_per->SetLineColor(kBlack);
	heff_cent->Draw("same");
	f_EfficiencyCent->Draw("same");
	f_EfficiencyPer->Draw("same");
	
	Double_t pTval[6]={0.2,0.5,1.0,1.2,1.5,2.0};
	for(int i=0;i<6;i++)
	{
		cout<<"pT: "<<pTval[i]<<" central: "<<f_EfficiencyCent->Eval(pTval[i])<<" peripheral: "<<f_EfficiencyPer->Eval(pTval[i])<<endl;
	}
	
	
}