Double_t Eff_track_rec_function(Double_t* x,Double_t* par)
{
    // Track reconstruction efficiency parametrization
    Double_t pt,y;
    Double_t A,B,C;A=par[0];B=par[1];C=par[2];
    pt=x[0];
    y=A*(exp(-pow(B/pt,C)));
    return y;
}

//============================================================

Int_t Apply_mom_smearing_and_efficiency(Int_t kflab_prim_glob, Float_t qp_In, Int_t kCentrality, Int_t i_Particle_use_In, Float_t m2_In,
                                        Float_t ePYTHIA_eff_factor_In, TLorentzVector TLV_Particle_In, PseudoJet& Fill_PseudoJet_smear_out,
                                        TF1* f_EfficiencyVsPt_In[][7][6])
{
    Int_t PID_eff = 2; // PID: p, anti-p, pi+, pi-, K+, K-
    if(qp_In > 0.0)
    {
        if(m2_In < 1.0 && m2_In > 0.8 ) PID_eff = 0; // p+
        if(m2_In < 0.3 && m2_In > 0.2 ) PID_eff = 4; // K+
        if(m2_In < 0.1 && m2_In > -0.1) PID_eff = 2; // pi+
    }
    if(qp_In < 0.0)
    {
        if(m2_In < 1.0 && m2_In > 0.8 ) PID_eff = 1; // p-
        if(m2_In < 0.3 && m2_In > 0.2 ) PID_eff = 5; // K-
        if(m2_In < 0.1 && m2_In > -0.1) PID_eff = 3; // pi-
    }


    // Centrality: 0-5, 5-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80
    Int_t cent9_eff = 8;
    if(kCentrality == 0) // 0-10% -> randomly choose between efficiency for 0-5% and 5%-10%
    {
        cent9_eff = 0;
        if(ran_gen.Rndm() < 0.5) cent9_eff = 1;
    }
    if(kCentrality == 1) // 60-80% -> randomly choose between efficiency for 60%-70% and 70%-80%
    {
        cent9_eff = 7;
        if(ran_gen.Rndm() < 0.5) cent9_eff = 8;
    }


    Int_t RunId_In     = 0; //all runIds together
    Double_t epsilon = ePYTHIA_eff_factor_In*f_EfficiencyVsPt_In[cent9_eff][RunId_In][PID_eff]->Eval(TLV_Particle_smear.Pt());
    Double_t rnd = gRandom->Uniform(0,1);
    // cout << "pT = " << TLV_Particle_smear.Pt() << ", epsilon = " << epsilon << ", rnd = " << rnd << endl;

    if(rnd <= epsilon)
    {
        Fill_PseudoJet_smear_out = Fill_PseudoJet_smear;
        //cout << "Particle was accepted, efficiency: " << epsilon << ", pT: " << TLV_Particle_smear.Pt() << ", PID: " << PID_eff << endl;
        return 1;
    }
    else
    {
        //cout << "Particle was rejected due to track reconstruction efficiency" << endl;
        return 0;
    }
}

//==========================================================
    //----------------------------------------------------------------------
    // Parameters for track reconstruction efficiency
    // Centrality: 0-5,5-10,10-20,20-30,30-40,40-50,50-60,60-70,70-80
    // runID: Index 0 is for all run ranges, then starting from 1: 12138025,12145021,12152017,12154022,12165032,12171017
    // PID: p, anti-p", pi+, pi-, K+, K-
    Double_t Global_ABC_Parameters_A[9][7][6];  //centrality(9),runID(7),PID(6)
    Double_t Global_ABC_Parameters_B[9][7][6];  //centrality(9),runID(7),PID(6)
    Double_t Global_ABC_Parameters_C[9][7][6];  //centrality(9),runID(7),PID(6)
    TString eEff_file="EffParameters2011_New.root";

//==========================================================

Double_t* StJetAnalysis::Parameter_eff_track_rec_function(Int_t Centrality, Int_t runID, Int_t PID, Int_t All)
{
    Int_t runID_range = -1;
    Double_t *Parameter = new Double_t[4];
    for(int i = 0; i < 4; i++)
    {
        Parameter[i] = 0;
    }

    runID_range = Get_runID_range_Index(runID);
    if(runID_range < 0)
    {
        cout << "WARNING: runID is out of range" << endl;
        return 0;
    }
    if( PID < 0 || PID > 5)
    {
        cout << "PID is out of range" << endl;
        return 0;
    }

    runID_range += 1; // index 0 is for all runIDs together
    if(All == 1){runID_range = 0;}

    Parameter[0] = Global_ABC_Parameters_A[Centrality][runID_range][PID];
    Parameter[1] = Global_ABC_Parameters_B[Centrality][runID_range][PID];
    Parameter[2] = Global_ABC_Parameters_C[Centrality][runID_range][PID];
    Parameter[3] = runID_range;

    return Parameter;
}

//==========================================================

void StJetAnalysis::Read_ABC_params_eff_track_rec_function()  //read ABC Parameter from a root file
{
    Int_t Switch_Array[6]={2,3,4,5,0,1}; //change PID in turn

    cout << "Read track reconstruction efficiency parameter file: " << eEff_file.Data() << endl;
    TFile *Eff_file=new TFile(eEff_file.Data());
    TH3D* hA  =(TH3D*)Eff_file->Get("hA");
    TH3D* hB  =(TH3D*)Eff_file->Get("hB");
    TH3D* hC  =(TH3D*)Eff_file->Get("hC");
    TH2D* hsA =(TH2D*)Eff_file->Get("hsA");
    TH2D* hsB =(TH2D*)Eff_file->Get("hsB");
    TH2D* hsC =(TH2D*)Eff_file->Get("hsC");

    for(Int_t i = 0; i < 6; i++) // PID
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
}

//==========================================================

cout << "Initialize track reconstruction efficiency functions" << endl;
    TString FuncName;
    for(Int_t loop_Centr = 0; loop_Centr < 9; loop_Centr++) // Centrality
    {
        for(Int_t loop_runID = 0; loop_runID < 7; loop_runID++) // runID
        {
            for(Int_t loop_PID = 0; loop_PID < 6; loop_PID++) // PID
            {
                Int_t loop_runID_use = loop_runID - 1; // runIDs start from loop_runID = 1, for loop_runID = 0 all runIDs are used
                Int_t use_all_runIDs = 0;
                if(loop_runID == 0) // for first index use all runIDs
                {
                    loop_runID_use = 0;
                    use_all_runIDs = 1;
                }
                Double_t* Parameter = Parameter_eff_track_rec_function(loop_Centr,Array_runID_eff_track_rec[loop_runID_use+1],loop_PID,use_all_runIDs);

                FuncName =  "f_EfficiencyVsPt_runID_Centr_";
                FuncName += loop_Centr;
                FuncName += "_runID_";
                FuncName += loop_runID;
                FuncName += "_PID_";
                FuncName += Array_PID_eff_track_rec[loop_PID];

                f_EfficiencyVsPt[loop_Centr][loop_runID][loop_PID] = new TF1(FuncName.Data(),Eff_track_rec_function,0,50.0,3);
                f_EfficiencyVsPt[loop_Centr][loop_runID][loop_PID]->SetParameter(0,Parameter[0]);
                f_EfficiencyVsPt[loop_Centr][loop_runID][loop_PID]->SetParameter(1,Parameter[1]);
                f_EfficiencyVsPt[loop_Centr][loop_runID][loop_PID]->SetParameter(2,Parameter[2]);

                delete[] Parameter;
            }
        }


//==========================================================
//particle ratio in pp
//==========================================================
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
