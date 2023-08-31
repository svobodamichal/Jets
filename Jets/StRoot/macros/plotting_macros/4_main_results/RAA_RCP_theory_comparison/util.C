TH1D *getPythiaHisto(char *file="pythia_fullEta1/histos_pythiajet_R0.4.root",float norm=1,float binsize=1,float ptlead=0){
  TFile *fp=new TFile(file,"read");
  TH1D *hpt=(TH1D*)fp->Get(Form("hpT_pTl%g",ptlead));
  hpt->Scale(norm);
  int nrebin=binsize/hpt->GetBinWidth(1);
  printf("rebin full jet pt spectra by %d.\n",nrebin);
  hpt->Rebin(nrebin);hpt->Scale(1./nrebin);
  return(hpt);
}

TGraphErrors *gReadData(char *file="star_6m_200_0208_sc1.dat",float scale=1,float err=0){
  FILE *fp=fopen(file,"read");
  float x[1000],y[1000],e[1000];
  int n=0;
  while(fscanf(fp,"%e %e",x+n,y+n)==2)++n;
  //for(int i=0;i<n;++i)printf("%e %e\n",x[i],y[i]);
  for(int i=0;i<n;++i){y[i]*=scale;e[i]=err*y[i];}
  //for(int i=0;i<n;++i){double y0=13.861*exp(-0.2929*x[i]);y[i]/=y0;e[i]/=y0;}
  TGraphErrors *g=new TGraphErrors(n,x,y,0,e);
  return(g);
}

TGraphErrors *gReadData2(char *file="star_6m_200_0208_sc1.dat",float scale=1,float err=0){
  FILE *fp=fopen(file,"read");
  float x[1000],y[1000],ex[1000],ey[1000];
  int n=0;
  while(fscanf(fp,"%e %e",x+n,y+n)==2)++n;
  for(int i=0;i<n-1;++i){ex[i]=(x[i+1]-x[i])/2;}ex[n-1]=(x[n-1]-x[n-2])/2;
  for(int i=0;i<n;++i){y[i]*=scale;ey[i]=err*y[i];}
  TGraphErrors *g=new TGraphErrors(n,x,y,ex,ey);
  return(g);
}

TGraphErrors *gReadData3(char *file="file.dat",float xbin_width=0){
  FILE *fp=fopen(file,"read");
  float x[1000],y[1000],y_low[1000],y_high[1000],ex[1000],ey[1000];
  int n=0;
  while(fscanf(fp,"%e %e %e",x+n,y_high+n,y_low+n)==3)++n;
  for(int i=0;i<n;++i)
  {
	  ex[i]=xbin_width/2;
	  ey[i]=(y_high[i]-y_low[i])/2;
	  y[i]=y_low[i]+ey[i];
	}
  TGraphErrors *g=new TGraphErrors(n,x,y,ex,ey);
  return(g);
}

TGraphErrors *gReadData4(char *file="file.dat",float xbin_width=0){
  FILE *fp=fopen(file,"read");
  float x[1000],y[1000],y_low[1000],y_high[1000],ex[1000],ey[1000],skip;
  int n=0;
  while(fscanf(fp,"%e %e %e %e %e %e %e",x+n,&skip, &skip, &skip, &skip,y_high+n,y_low+n) != EOF)++n;
  for(int i=0;i<n;++i)
  {
	  ex[i]=xbin_width/2;
	  ey[i]=(y_high[i]-y_low[i])/2;
	  y[i]=y_low[i]+ey[i];
	}
  TGraphErrors *g=new TGraphErrors(n,x,y,ex,ey);
  return(g);
}

TH1D *g2h(TGraphErrors *g=0,TH1D *hcp=0){
  float x[1000],y[1000];
  double xx,yy;
  int n=g->GetN();
  for(int i=0;i<n;++i){
    g->GetPoint(i,xx,yy);
    //printf("%e %e\n",xx,yy);
    x[i]=xx;y[i]=yy;
  }
  TH1D *h=(TH1D*)hcp->Clone();h->Reset();
  for(int k=1;k<=h->GetNbinsX();++k){
    float xh=h->GetBinCenter(k);
    for(int i=0;i<n-1;++i){
      if(xh<x[i]||xh>x[i+1])continue;
      //float yh=exp(log(y[i])+(xh-x[i])*(log(y[i+1])-log(y[i]))/(x[i+1]-x[i]));
      float yh=y[i]*pow(y[i+1]/y[i],(xh-x[i])/(x[i+1]-x[i]));
      h->SetBinContent(k,yh);h->SetBinError(k,0);
    }
  }
  //h->Draw();g->SetMarkerStyle(24);g->Draw("lp");
  return(h);
}

TGraphErrors *draw_hist_graph(TH1D *h,char *option="l",float xmin=0,float xmax=100){
  //const int N=h->GetNbinsX();
  //float x[N],y[N],e[N];
  int N=h->GetNbinsX();
  float x[1000],y[1000],e[1000];
  int n=0;
  for(int i=0;i<N;++i){
    float xtmp=h->GetBinCenter(i+1);
    if(xtmp<xmin)continue;
    if(xtmp>xmax)continue;
    x[n]=xtmp;
    y[n]=h->GetBinContent(i+1);
    e[n]=h->GetBinError(i+1);
    n++;
  }
  if(strstr(option,"p"))TGraphErrors *g=new TGraphErrors(n,x,y,0,e);
  else TGraphErrors *g=new TGraphErrors(n,x,y,0,0);
  g->SetMarkerStyle(h->GetMarkerStyle());
  g->SetMarkerColor(h->GetMarkerColor());
  g->SetLineStyle(h->GetLineStyle());
  g->SetLineWidth(h->GetLineWidth());
  g->SetLineColor(h->GetLineColor());
  g->Draw(option);
  return(g);
}
