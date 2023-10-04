void plot_JER_AllinOne(float R=0.3, float pTlead=5.0,TString type="normal",TString evo="omicron",TString ext="pdf")
{
	
	int can_x=300; //canvas panel size
	int can_y=250;
	
	float xmin=0;
	float xmax=70;
	float ymin=0.005;
	float ymax=0.15;
	bool plotLogY=1; //plot log y-axis
	
	if(type=="u" || type=="g" )
	{
		ymax=1.0;
		ymin=0.03;
	}
	
	const int nprojections=6;
	const float pTtrue[nprojections]={7,10,15,20,30,40}; //pTtrue for which we want to show the pTdete distribution
	
	Color_t colors[nprojections]={kBlack, kGray+2,kGray+1,kBlack, kGray+2,kGray+1};
	
	const float fit_xmin[nprojections]={6,8,13,17.5,27,36};
	const float fit_xmax[nprojections]={8.5,12,17,22.5,33,44};
	
	//latex coordinates
	float latex_x[nprojections]={0.25,0.29,0.36,0.44,0.56,0.7};
	float latex_y[nprojections]={0.8,0.6,0.47,0.36,0.27,0.22};
	float latex_y_log[nprojections]={0.82,0.72,0.6,0.53,0.45,0.36}; //in case of log y-axis
	for(int pr=0;pr<nprojections;pr++)
	{
		latex_x[pr]=latex_x[pr]*60/xmax;
		latex_y[pr]=latex_y[pr];
	}
	
	TString frag="u-quark";
	if(type=="g")frag="gluon";
	else if(type=="normal")frag="u-quark + gluon (2:1)";
	TString figpath=Form("../../../plotting_out/obr/%s/effi/",evo.Data());
	TString path=Form("../../../plotting_out/intersteps/epsilon/%s/central/%s",evo.Data(),type.Data());
	
	//define histograms
	TH2D* hpTdete2D;
	TH1D* hpTdete[nprojections];
	
	TString name=Form("%s/pythia_emb_R%.1lf.root",path.Data(),R);
	infile=new TFile(name,"OPEN");
	name=Form("hresponse_pTl%.0lf",pTlead);
	hpTdete2D=(TH2D*) infile->Get(name);
	
	TCanvas* cprojects_all=new TCanvas("cprojects_all","projections",10,10,2*can_x,2*can_y);
	cprojects_all->cd();
	if(plotLogY) cprojects_all->SetLogy();
	TH1D* frame=new TH1D("hframe","hframe",2,xmin,xmax);
	frame->SetTitle("");
	frame->GetXaxis()->SetTitle("p_{T,jet}^{det} (GeV/c)");
	frame->GetYaxis()->SetTitle("probability");
	frame->GetYaxis()->SetRangeUser(ymin,ymax);
	frame->Draw(""); 
	
	TLatex *latex = new TLatex();
	latex->SetNDC();
	latex->SetTextSize(0.03);
			
	for(int i=0; i<nprojections;i++)
	{
		int bin=hpTdete2D->GetYaxis()->FindBin(pTtrue[i]);
		name=Form("projection_%.0lf",pTtrue[i]);
		hpTdete[i]=hpTdete2D->ProjectionX(name,bin,bin);
		hpTdete[i]->SetTitle("");
		hpTdete[i]->GetXaxis()->SetTitleSize(0.06);
		hpTdete[i]->SetLineColor(colors[i]);
		hpTdete[i]->Rebin(4);
		hpTdete[i]->Scale(1.0/hpTdete[i]->Integral());
		hpTdete[i]->Draw("same");
		
		//gaussian fit
		TF1 * fgaus = new TF1("gauss", "[0] / sqrt(2.0 * TMath::Pi()) / [2] * exp(-(x-[1])*(x-[1])/2./[2]/[2])", fit_xmin[i], fit_xmax[i]);
		fgaus->SetParNames("Constant","Mean","Sigma");
		fgaus->SetParameters(hpTdete[i]->Integral(),pTtrue[i],hpTdete[i]->GetRMS());
		//fgaus->SetParLimits(0,  1E2, 1E4);
		//fgaus->SetParLimits(1, -0.2,  0.1);
		//fgaus->SetParLimits(2, 1E-2,  5E-1);
		hpTdete[i]->Fit(fgaus,"R","",fit_xmin[i],fit_xmax[i]);
		float sigma=fgaus->GetParameter(2);
		float res=TMath::Abs(100*sigma/pTtrue[i]);
		
		latex->DrawLatex(latex_x[i],(plotLogY)?latex_y_log[i]:latex_y[i],Form("p_{T,jet}^{part}=%.0lf GeV/c, #delta p_{T}/p_{T}=%.0lf%%",pTtrue[i],res));
		
		TLegend *linfo = new TLegend(0.6, 0.55, 0.8, 0.85);
		linfo->SetFillStyle(0);
		linfo->SetBorderSize(0);
		linfo->SetTextSize(0.04);
		linfo->SetMargin(0.05);
		linfo->AddEntry("","PYTHIA 6","");
		linfo->AddEntry("",Form("%s",frag.Data()),"");
		linfo->AddEntry("", Form("anti-k_{T}, R = %.1lf",R), "");
		//	linfo->AddEntry("", Form("p_{T}^{lead} > %.1lf GeV/c",pTlead), ""); //Jana
		linfo->AddEntry("", Form("p_{T,lead}^{min} = %.0lf GeV/c",pTlead), ""); //Jana
		linfo->AddEntry("","det. effects for","");
		linfo->AddEntry("","  central Au+Au","");
		//linfo->AddEntry("", "p_{T}^{const} > 0.2 GeV/c", "");
		linfo->DrawClone("same");
	}//projections
	cprojects_all->SaveAs(Form("%s/JER_allinone_R0%.0lf_pTlead%.0lf_%s.%s",figpath.Data(),R*10,pTlead,type.Data(),ext.Data()));
	
}

//=======================================================
//plot JET ENERGY RESOLUTION pT dependence
//=======================================================
void plot_JER(float pTlead=5.0,TString type="normal",bool doSys=1)
{
	
	const int nprojections=9;
	const float low_x[nprojections]={pTlead,10,15,20,25,30,35,40,45}; //projection bin boundaries
	const float high_x[nprojections]={10,15,20,25,30,35,40,45,50};
	
	
	//number of datasets - 1st is central value, 2 others for systematic error calculation
	int nsets=1;
	if(doSys)nsets=3; //if we want systematic errors
	
	const int nr=2; //number of R parameters
	float R[nr]={0.2,0.4};
	const int ndatasets=nsets; 
	float mean[nr][nprojections]; //mean
	float mean_down[nr][nprojections]; //sys. error up
	float mean_up[nr][nprojections]; //sys. error down
	float mpv[nr][nprojections]; //most probabale value
	float mpv_down[nr][nprojections]; 
	float mpv_up[nr][nprojections]; 
	float q1[nr][nprojections]; //1st quartile
	float q1_up[nr][nprojections];
	float q1_down[nr][nprojections];
	float q2[nr][nprojections]; //median
	float q2_up[nr][nprojections];
	float q2_down[nr][nprojections];
	float q3[nr][nprojections]; //3rd quartile
	float q3_up[nr][nprojections];
	float q3_down[nr][nprojections];
	
	int can_x=300; //canvas panel size
	int can_y=250;
	
	//axis range
	float xmin=0;
	float xmax=50;
	float ymin=-0.5;
	float ymax=0.5;
	float xshift=1.0; //shift in x axis in order to display clerly different R
	
	TString frag="u-quark";
	if(type=="g")frag="gluon";
	else if(type=="normal")frag="u-quark + gluon (2:1)";
	TString figpath="obr";
	
	float med_x[nr][nprojections]; //x axis bin centers
	for(int i=0;i<nprojections;i++)
	{
		med_x[0][i]=(low_x[i]+high_x[i])/2.0;
		med_x[1][i]=(low_x[i]+high_x[i])/2.0+xshift;
	}
	
	//define histograms
	TH2D* hJER2D[ndatasets];
	TH1D* hJER[nprojections][ndatasets];
	
	TCanvas* cprojects[nr][ndatasets];
	
	TFile* infile[ndatasets];//input fiels
	TString path[]={Form("%s/",type.Data()),"m5/","p5/"};
	for(int r=0; r<nr;r++)
	{
		for (int j=0;j<ndatasets;j++)
		{
			TString name=Form("./%spythia_emb_R%.1lf.root",path[j].Data(),R[r]);
			infile[j]=new TFile(name,"OPEN");
			name=Form("hJER_pTl%.0lf",pTlead);
			hJER2D[j]=(TH2D*) infile[j]->Get(name);
			
			cprojects[r][j]=new TCanvas(Form("cprojects%i_R%i",j,r),Form("projections_%i_R0%.0lf",j,R[r]*10),10,10,3*can_x,3*can_y);
			cprojects[r][j]->Divide(3,3);
			TLatex *latex = new TLatex();
			latex->SetNDC();
			latex->SetTextSize(0.05);
			
			for(int i=0; i<nprojections;i++)
			{
				
				int bin_low=hJER2D[j]->GetXaxis()->FindBin(low_x[i]);
				int bin_high=hJER2D[j]->GetXaxis()->FindBin(high_x[i]);
				name=Form("projection_%.0lf-%.0lf",low_x[i],high_x[i]);
				hJER[i][j]=hJER2D[j]->ProjectionY(name,bin_low,bin_high);
				hJER[i][j]->SetTitle("");
				hJER[i][j]->GetXaxis()->SetTitleSize(0.06);
				//if(i==1) hJER[i][j]->SetTitle(Form("R=%.1lf",R[r]));
				
				cprojects[r][j]->cd(i+1);
				hJER[i][j]->Draw("e histo");
				
				//calculate most probable value as the mean of the gaussian fit
				TF1 * fgaus = new TF1("gauss", "[0] / sqrt(2.0 * TMath::Pi()) / [2] * exp(-(x-[1])*(x-[1])/2./[2]/[2])", -1, 1);
				fgaus->SetParNames("Constant","Mean","Sigma");
				fgaus->SetParameters(hJER[i][j]->Integral(),hJER[i][j]->GetMean(),hJER[i][j]->GetRMS());
				fgaus->SetParLimits(0,  1E2, 1E4);
				fgaus->SetParLimits(1, -0.2,  0.1);
				fgaus->SetParLimits(2, 1E-2,  5E-1);
				hJER[i][j]->Fit(fgaus,"","",hJER[i][j]->GetMean()/4-0.08,hJER[i][j]->GetMean()/4+0.08);
				
				if(i==0){
					TLegend *linfo = new TLegend(0.15, 0.3, 0.6, 0.9);
					linfo->SetFillStyle(0);
					linfo->SetBorderSize(0);
					linfo->SetTextSize(0.05);
					linfo->SetMargin(0.05);
					linfo->SetHeader("PYTHIA 6");
					linfo->AddEntry("",Form("frag.: %s",frag.Data()),"");
					linfo->AddEntry("", Form("anti-k_{T}, R = %.1lf",R[r]), "");
					linfo->AddEntry("", Form("p_{T}^{lead} > %.1lf GeV/c",pTlead), "");
					linfo->AddEntry("", "p_{T}^{const} > 0.2 GeV/c", "");
					linfo->DrawClone("same");
				}
				
				
				//calculate median and 1st and 3rd quartile
				float sum=0;
				float total=hJER[i][j]->Integral();
				//cout<<"integral: "<<hJER[i][j]->Integral()<<" nentries: "<<hJER[i][j]->GetEntries()<<endl;
				bool showq3=1;
				float q1_tmp=0;
				float q3_tmp=0;
				float q2_tmp=0;
				
				for(int bin=1;bin<=hJER[i][j]->GetNbinsX();bin++)
				{
					float x=hJER[i][j]->GetBinCenter(bin);
					sum+=hJER[i][j]->GetBinContent(bin);
					float ratio=(sum/total);
					//cout<<"bin:"<<bin<<" ratio:"<<ratio<<endl;
					if(ratio<=.25) q1_tmp=x;
					if(ratio<=.50) q2_tmp=x;
					if(ratio<=.75) q3_tmp=x;
					
				}//bin loop
				
				//save mean value or upper or lower error band
				switch(j)
				{
					case 0:
						mean[r][i]=hJER[i][j]->GetMean();
						mpv[r][i]=fgaus->GetParameter(1);
						q1[r][i]=q1_tmp;
						q2[r][i]=q2_tmp;
						q3[r][i]=q3_tmp;
						break;
					case 1:
						mean_up[r][i]=hJER[i][j]->GetMean();
						mpv_up[r][i]=fgaus->GetParameter(1);
						q1_up[r][i]=q1_tmp;
						q2_up[r][i]=q2_tmp;
						q3_up[r][i]=q3_tmp;
						break;
					case 2:
						mean_down[r][i]=hJER[i][j]->GetMean();
						mpv_down[r][i]=fgaus->GetParameter(1);
						q1_down[r][i]=q1_tmp;
						q2_down[r][i]=q2_tmp;
						q3_down[r][i]=q3_tmp;
						break;
				}
				
				latex->DrawLatex(0.65, 0.8,Form("%.0lf-%.0lf GeV/c",low_x[i],high_x[i]));
			}//projection loop
			cprojects[r][j]->SaveAs(Form("%s/projection_%s_set%i_R0%.0lf_pTlead%.0lf.pdf",figpath.Data(),type.Data(),j,R[r]*10,pTlead));
			cprojects[r][j]->SaveAs(Form("%s/projection_%s_set%i_R0%.0lf_pTlead%.0lf.gif",figpath.Data(),type.Data(),j,R[r]*10,pTlead));
		}//dataset loop
	}//R loop
	
	//Declare TGraphs
	TGraph *gr_mean[nr];
	TGraph *gr_mpv[nr];
	TGraph *gr_median[nr];
	TGraph *gr_q1[nr];
	TGraph *gr_q3[nr];
	
	Color_t colorList[]={kBlue,kRed,kBlack,kMagenta};
	
	for(int r=0;r<nr;r++)
	{
		
		float *med_xr=&med_x[r][0];
		float *meanr=&mean[r][0];
		float *mpvr=&mpv[r][0];
		float *q2r=&q2[r][0];
		float *q1r=&q1[r][0];
		float *q3r=&q3[r][0];
		
		
		gr_mean[r] = new TGraph(nprojections,med_xr,meanr);
		gr_mean[r]->SetMarkerStyle(20);
		gr_mean[r]->SetMarkerColor(colorList[r]);
		
		gr_mpv[r] = new TGraph(nprojections,med_xr,mpvr);
		gr_mpv[r]->SetMarkerStyle(21);
		//gr_mpv[r]->SetMarkerSize(2);
		gr_mpv[r]->SetMarkerColor(colorList[r]);
		
		gr_median[r] = new TGraph(nprojections,med_xr,q2r);
		gr_median[r]->SetMarkerStyle(34);
		gr_median[r]->SetMarkerColor(colorList[r]);
		
		gr_q1[r] = new TGraph(nprojections,med_xr,q1r);
		gr_q1[r]->SetMarkerStyle(22);
		gr_q1[r]->SetMarkerColor(colorList[r]);
		
		gr_q3[r] = new TGraph(nprojections,med_xr,q3r);
		gr_q3[r]->SetMarkerStyle(23);
		gr_q3[r]->SetMarkerColor(colorList[r]);
		
	}//R loop
	
	gr_mean[0]->GetXaxis()->SetLimits(xmin,xmax);
	gr_mean[0]->GetHistogram()->SetMaximum(ymax);   // along          
	gr_mean[0]->GetHistogram()->SetMinimum(ymin);
	gr_mean[0]->SetTitle("");
	gr_mean[0]->GetYaxis()->SetTitle("(p_{T}^{dete}-p_{T}^{part})/p_{T}^{part}");
	gr_mean[0]->GetXaxis()->SetTitle("p_{T}^{part} [GeV/c]");
	
	TCanvas* cJER=new TCanvas("cJER","JER",10,10,3*can_x,2*can_y);
	TString drawopt[]={"ap","p"};
	for(int r=0;r<nr;r++)
	{
		gr_mean[r]->Draw(drawopt[r]);
		gr_median[r]->Draw("p");
		gr_q1[r]->Draw("p");
		gr_q3[r]->Draw("p");
		gr_mpv[r]->Draw("p");
		
	}//R loop
	TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.90);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.035);
	leg->AddEntry(gr_mean[0], "mean", "p");
	leg->AddEntry(gr_mpv[0], "m.p.v.", "p");
	leg->AddEntry(gr_q1[0], "prob(ratio<y)=25%", "p");
	leg->AddEntry(gr_median[0], "prob(ratio<y)=50% (median)", "p");
	leg->AddEntry(gr_q3[0], "prob(ratio<y)=75%", "p");
	leg->AddEntry("", "", "");
	for(int r=0;r<nr;r++)
	{
		leg->AddEntry(gr_mean[r], Form("R=%.1lf",R[r]), "p");
	}
	leg->DrawClone("same");
	
	TLegend *model_info = new TLegend(0.2, 0.65, 0.4, 0.9);
	model_info->SetFillStyle(0);
	model_info->SetBorderSize(0);
	model_info->SetTextSize(0.035);
	model_info->SetMargin(0.05);
	model_info->SetHeader("PYTHIA 6");
	model_info->AddEntry("",Form("fragmentation: %s",frag.Data()),"");
	model_info->AddEntry("","anti-k_{T}","");
	//model_info->AddEntry("", Form("Anti-k_{T} \t \t R = %.1lf",R[0]), "");
	model_info->AddEntry("", Form("p_{T}^{lead} > %.1lf GeV/c",pTlead), "");
	model_info->AddEntry("", "p_{T}^{const} > 0.2 GeV/c", "");
	model_info->DrawClone("same");
	cJER->SaveAs(Form("%s/JER_%s_pTlead%.0lf.pdf",figpath.Data(),type.Data(),pTlead));
	cJER->SaveAs(Form("%s/JER_%s_pTlead%.0lf.gif",figpath.Data(),type.Data(),pTlead));
	
	
	//print out JER values
	cout<<"================================"<<endl;
	cout<<"JER"<<endl;
	cout<<"================================"<<endl;
	for(int r=0;r<nr;r++)
	{
		cout<<"--------------"<<endl;
		cout<<"R="<<R[r]<<endl;
		cout<<"--------------"<<endl;
		for(int prj=0; prj<nprojections; prj++)
		{
			float pt=med_x[r][prj];
			if(pt>xmax)continue;
			float JER_down=100*(mpv[r][prj]-q1[r][prj]);
			float JER_up=-100*(mpv[r][prj]-q3[r][prj]);
			cout<<"pT: "<<pt-r*xshift<<" MPV: "<<100*mpv[r][prj]<<"%, JER up: "<<JER_up<<"% JER down: "<<JER_down<<"%"<<endl;	
			if(doSys)
			{
				float JER_down_e1=100*(mpv_up[r][prj]-q1_up[r][prj]);
				float JER_down_e2=100*(mpv_down[r][prj]-q1_down[r][prj]);
				float JER_up_e1=-100*(mpv_up[r][prj]-q3_up[r][prj]);
				float JER_up_e2=-100*(mpv_up[r][prj]-q3_up[r][prj]);
				
				float err_down_1=TMath::Abs(JER_down-JER_down_e1);
				float err_down_2=TMath::Abs(JER_down-JER_down_e2);
				float err_up_1=TMath::Abs(JER_up-JER_up_e1);
				float err_up_2=TMath::Abs(JER_up-JER_up_e2);
				
				float err_down=(err_down_1>err_down_2) ? err_down_1 : err_down_2;
				float err_up=(err_up_1>err_up_2) ? err_up_1 : err_up_2;
				cout<< "                                                    JER_up (abs.) err.: "<<err_up<<"%, JER_down (abs.) err.: "<<err_down<<"%"<<endl;
			}
		}//pT bin loop
	}//r loop
	
	
	
}//end            
