void multiply_matrix()
{
Float_t pTlead = atof(gSystem->Getenv("PTTHRESH"));
Float_t R=atof(gSystem->Getenv("RPARAM"));
Int_t reverse=atoi(gSystem->Getenv("REVERSE"));
//TString jettype=gSystem->Getenv("JETTYPE");
//TString suffix=gSystem->Getenv("CSUFFIX");
//TString sys_suf=gSystem->Getenv("SYSSUF");
//TString type_suf=gSystem->Getenv("TSUFFIX");
//TString type_suf2=gSystem->Getenv("TSUFFIX2");
TString input_path1=gSystem->Getenv("DIR_IN1");
TString input_path2=gSystem->Getenv("DIR_IN2");
TString output_path=gSystem->Getenv("DIR_OUT");

cout<<"R: "<<R<<endl;
cout<<"pTlead: "<<pTlead<<endl;

TString str;
//TString data_path = Form("/global/homes/r/rusnak/jet_analysis/STARJet/out/MB/embedding%s%s",suffix.Data(),sys_suf.Data());
//TString outdir=data_path;
str = Form("%s/response_matrix_BG_sp_R%.1lf_pTlead%.0lf.root", input_path1.Data(),R,pTlead);
cout<<"input1: "<<str<<endl;
TFile *finput1 = new TFile(str.Data(), "OPEN");
TH2D* hmatrix1=(TH2D*) finput1->Get("hResponse_1E9");
//TH2D* hmatrix1=(TH2D*) finput->Get(Form("hdpT2D_%s",type.Data()));
//str = Form("%s/pythia_emb_R%.1lf.root", input_path2.Data(),R);
str = Form("%s/pythia6_normalized_combined_response.root", input_path2.Data());
cout<<"input2: "<<str<<endl;
//str = Form("%s/histos_deltapT_R%.1lf.root", data_path.Data(),R);
TFile *finput2 = new TFile(str.Data(), "OPEN");
//cout << "matrix name " << Form("hResponseMatrix_pTl%.0lf_R0%.0lf",pTlead, R*10) << endl;
TH2D* hmatrix2=(TH2D*)finput2->Get(Form("hResponseMatrix_pTl%.0lf_R0%.0lf",pTlead, R*10));
//TH2D* hmatrix2=(TH2D*) finput2->Get("hdpT2D_dete");


TString outfile=Form("%s/response_matrix_BGD_R%.1lf_pTlead%.0lf.root",output_path.Data(),R,pTlead);
cout<<"output: "<<outfile<<endl;
TFile* fout = new TFile(outfile.Data(), "recreate");

TH2D *hmatrix_out=(TH2D*) hmatrix2->Clone("hmatrix_out");
hmatrix_out->Reset("MICE");
Int_t nbins=hmatrix1->GetNbinsX();
//cout << nbins << " " << hmatrix2->GetNbinsX() << endl;
TMatrixD matrix1(nbins,nbins);
TMatrixD matrix2(nbins,nbins);
TMatrixD matrixo(nbins,nbins);

cout<<"filling vectors"<<endl;
for(Int_t i=0; i<nbins; i++){
for(Int_t j=0; j<nbins; j++){
 matrix1(i,j)=hmatrix1->GetBinContent(i + 1, j + 1);
 matrix2(i,j)=hmatrix2->GetBinContent(i + 1, j + 1);
}
}

if(!reverse)
{
	cout<<"multiplying BGxdete"<<endl;
	matrixo.Mult(matrix1,matrix2);
}
else 
{
	cout<<"multiplying detexBG"<<endl;
	matrixo.Mult(matrix2,matrix1);
}
for(Int_t i=0; i<nbins; i++){
for(Int_t j=0; j<nbins; j++){
hmatrix_out->SetBinContent(i+1, j+1, matrixo(i,j));
}
}

cout<<"writing output"<<endl;
hmatrix_out->Write("hResponse_1E9");
//hmatrix_out->Write("hdpT2D_BG_dete");

finput1->Close();
finput2->Close();
fout->Close();
}

