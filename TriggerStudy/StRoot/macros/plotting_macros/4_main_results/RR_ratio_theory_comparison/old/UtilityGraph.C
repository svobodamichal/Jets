/*
TGraphErrors *Gjoin(TGraphErrors *g1,TGraphErrors *g2)
{
  int n1=g1->GetN();
  int n2=g2->GetN();
  int n=n1+n2;
  //printf("%d %d %d\n",n1,n2,n);
  float x[1000],y[1000],ex[1000],ey[1000];
  if(n>1000){printf("Gjoin: n=%d error!\n",n);return(0);}
  double x1,y1;
  for(int k=0;k<n1;++k){
    g1->GetPoint(k,x1,y1);
    x[k]=x1;
    y[k]=y1;
    ex[k]=g1->GetErrorX(k);
    ey[k]=g1->GetErrorY(k);
  }
  for(int k=0;k<n2;++k){
    g2->GetPoint(k,x1,y1);
    x[n1+k]=x1;
    y[n1+k]=y1;
    ex[n1+k]=g2->GetErrorX(k);
    ey[n1+k]=g2->GetErrorY(k);
  }
  TGraphErrors *g=new TGraphErrors(n,x,y,ex,ey);
  return(g);
}
*/

TGraphErrors *Gjoin(TGraphErrors *g1,TGraphErrors *g2,int n1=0,int n2=0)
{
  int n1a=0,n1b=g1->GetN();
  int n2a=0,n2b=g2->GetN();
  if(n1>0){n1b=n1;}else if(n1<0){n1a=n1b+n1;}
  if(n2>0){n2b=n2;}else if(n2<0){n2a=n2b+n2;}
  int n=n1b-n1a+n2b-n2a;
  //printf("%d %d %d %d\n",n1a,n1b,n2a,n2b);
  float x[1000],y[1000],ex[1000],ey[1000];
  if(n>1000){printf("Gjoin: n=%d error!\n",n);return(0);}
  double x1,y1;
  int kk=0;
  for(int k=n1a;k<n1b;++k){
    g1->GetPoint(k,x1,y1);
    x[kk]=x1;
    y[kk]=y1;
    ex[kk]=g1->GetErrorX(k);
    ey[kk]=g1->GetErrorY(k);
    kk++;
  }
  for(int k=n2a;k<n2b;++k){
    g2->GetPoint(k,x1,y1);
    x[kk]=x1;
    y[kk]=y1;
    ex[kk]=g2->GetErrorX(k);
    ey[kk]=g2->GetErrorY(k);
    kk++;
  }
  TGraphErrors *g=new TGraphErrors(n,x,y,ex,ey);
  return(g);
}


void Gplot(TGraphErrors *g,int color=1,float mark=20.20,float size=1,float ykey=0,int frame=0,int exponent=0)
{
  //TGraphErrors *gtmp=(TGraphErrors*)g->Clone();
  TGraphErrors *gtmp=Gscale(g,pow(10,exponent));
  gtmp->SetName("gtmp");
  if(mark<0){for(int i=0;i<gtmp->GetN();++i)gtmp->SetPointError(i,0,0);}    
  if(mark<0)gtmp->Print();
  gtmp->SetLineColor(color);
  gtmp->SetMarkerColor(color);
  gtmp->SetMarkerSize(size);
  int mark1=int(mark);
  int mark2=int((mark-int(mark))*100+0.5);
  if(!mark2){
    gtmp->SetMarkerStyle(abs(mark1));
    if(frame)gtmp->Draw("ap");
    else gtmp->Draw("p");
    if(ykey>1e-6)keySymbol(0.2,ykey,gtmp->GetTitle(),color,abs(mark1),0.04,size);
  }
  if(!mark2)return;
  float x[1000],y[1000],ex[1000],ey[1000];
  double x1,y1;
  int N=gtmp->GetN();
  for(int k=0;k<N;++k){
    gtmp->GetPoint(k,x1,y1);
    x[k]=x1;
    y[k]=y1;
    ex[k]=gtmp->GetErrorX(k);
    ey[k]=gtmp->GetErrorY(k);
  }
  TGraphErrors *gtmp1=new TGraphErrors(N-1,x+1,y+1,ex+1,ey+1);
  TGraphErrors *gtmp2=new TGraphErrors(1,x,y,ex,ey);
  gtmp1->SetName("gtmp1");
  gtmp2->SetName("gtmp2");
  gtmp1->SetLineColor(color);
  gtmp2->SetLineColor(color);
  gtmp1->SetMarkerColor(color);
  gtmp2->SetMarkerColor(color);
  gtmp1->SetMarkerSize(size);
  gtmp2->SetMarkerSize(size);
  gtmp1->SetMarkerStyle(abs(mark1));
  gtmp2->SetMarkerStyle(abs(mark2));
  if(frame){gtmp1->Draw("ap");gtmp2->Draw("p");}
  else {gtmp1->Draw("p");gtmp2->Draw("p");}
  if(ykey>1e-6)keySymbol(0.2,ykey,gtmp->GetTitle(),color,abs(mark2),0.04,size);
  return;
}


void Gplot2(TGraphErrors *gy,TGraphErrors *gx,int color=1,int mark=20,float size=1,float ykey=0,int frame=0)
{
  float x[1000],y[1000],ex[1000],ey[1000];
  double x1,y1,x2,y2;
  int N=gx->GetN();
  for(int k=0;k<N;++k){
    gx->GetPoint(k,x1,y1);
    gy->GetPoint(k,x2,y2);
    x[k]=y1;
    y[k]=y2;
    ex[k]=gx->GetErrorY(k);
    ey[k]=gy->GetErrorY(k);
  }
  TGraphErrors *g=new TGraphErrors(N,x,y,ex,ey);
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(mark);
  g->SetMarkerSize(size);
  if(frame)g->Draw("ap");
  else g->Draw("p");
  if(ykey>1e-6)keySymbol(0.2,ykey,g->GetTitle(),color,mark,0.04,size);
  return;
}

TGraphErrors *Ggvsg(TGraphErrors *gy,TGraphErrors *gx, int i1=0,int i2=0)
{
  float x[1000],y[1000],ex[1000],ey[1000];
  double x1,y1,x2,y2;
  int N=gx->GetN();
  for(int k=i1;k<N;++k){
    gx->GetPoint(k,x1,y1);
    gy->GetPoint(k,x2,y2);
    x[k-i1]=y1;
    y[k-i1]=y2;
    ex[k-i1]=gx->GetErrorY(k);
    ey[k-i1]=gy->GetErrorY(k);
  }
  int n=N-i1;
  if(i2&&i2<N)n=i2-i1+1;
  TGraphErrors *g=new TGraphErrors(n,x,y,ex,ey);
  return(g);
}


void Gline(TGraphErrors *g,int color=1,int style=1,float width=1,float ykey=0,float size=1,int exponent=0)
{
  TGraphErrors *gtmp=Gscale(g,pow(10,exponent));
  //TGraphErrors *gtmp=(TGraphErrors*)g->Clone();
  gtmp->SetName("gtmp");
  for(int k=0;k<gtmp->GetN();++k)gtmp->SetPointError(k,0,0);
  gtmp->SetLineColor(color);
  gtmp->SetLineStyle(style);
  gtmp->SetLineWidth(width);
  gtmp->Draw("l");
  if(ykey>1e-6)keyLine(0.2,ykey,gtmp->GetTitle(),color,style,0.04,width,size);
  return;
}


TGraphErrors *Gscale(TGraphErrors *g1,float scale,float escale=0)
{
  float x[1000],y[1000],ex[1000],ey[1000];
  double x1,y1;
  int N=g1->GetN();
  for(int k=0;k<N;++k){
    g1->GetPoint(k,x1,y1);
    x[k]=x1;
    ex[k]=g1->GetErrorX(k);
    y[k]=scale*y1;
    ey[k]=sqrt(pow(g1->GetErrorY(k)/y1,2)+pow(escale/scale,2));
    ey[k]*=y[k];
  }
  TGraphErrors *g=new TGraphErrors(N,x,y,ex,ey);
  g->SetTitle(g1->GetTitle());
  g->SetMarkerColor(g1->GetMarkerColor());
  g->SetMarkerStyle(g1->GetMarkerStyle());
  g->SetLineColor(g1->GetLineColor());
  g->SetLineStyle(g1->GetLineStyle());
  g->SetLineWidth(g1->GetLineWidth());
  return(g);
}


TGraphErrors *Gabs(TGraphErrors *g1)
{
  float x[1000],y[1000],ex[1000],ey[1000];
  double x1,y1;
  int N=g1->GetN();
  for(int k=0;k<N;++k){
    g1->GetPoint(k,x1,y1);
    x[k]=x1;
    y[k]=fabs(y1);
    ex[k]=g1->GetErrorX(k);
    ey[k]=g1->GetErrorY(k);
  }
  TGraphErrors *g=new TGraphErrors(N,x,y,ex,ey);
  g->SetTitle(g1->GetTitle());
  return(g);
}


TGraphErrors *Gave(TGraphErrors *g1,TGraphErrors *g2)
{
  float x[1000],y[1000],ex[1000],ey[1000];
  double x1,y1,e1,x2,y2,e2;
  int N=g1->GetN();
  for(int k=0;k<N;++k){
    g1->GetPoint(k,x1,y1);
    g2->GetPoint(k,x2,y2);
    e1=g1->GetErrorY(k);
    e2=g2->GetErrorY(k);
    x[k]=x1;
    y[k]=(y1+y2)/2;
    ex[k]=g1->GetErrorX(k);
    ey[k]=sqrt(e1*e1+e2*e2)/2;
  }
  TGraphErrors *g=new TGraphErrors(N,x,y,ex,ey);
  g->SetTitle(g1->GetTitle());
  return(g);
}


TGraphErrors *Gave(int n,TGraphErrors **g)
{
  TGraphErrors *gave=g[0];
  for(int i=1;i<n;++i)gave=Gadd(gave,g[i]);
  gave=Gscale(gave,1./n);
  return(gave);
}


TGraphErrors *Gave(TGraphErrors *g0=0,TGraphErrors *g1=0,TGraphErrors *g2=0,TGraphErrors *g3=0,TGraphErrors *g4=0,TGraphErrors *g5=0,TGraphErrors *g6=0,TGraphErrors *g7=0,TGraphErrors *g8=0,TGraphErrors *g9=0)
{
  int n=0;
  TGraphErrors *g[10];
  if(g0){g[n]=g0;n++;}
  if(g1){g[n]=g1;n++;}
  if(g2){g[n]=g2;n++;}
  if(g3){g[n]=g3;n++;}
  if(g4){g[n]=g4;n++;}
  if(g5){g[n]=g5;n++;}
  if(g6){g[n]=g6;n++;}
  if(g7){g[n]=g7;n++;}
  if(g8){g[n]=g8;n++;}
  if(g9){g[n]=g9;n++;}
  TGraphErrors *gave=Gave(n,g);
  return(gave);
}


TGraph *GquadSum(TGraph *g1,TGraph *g2)
{
  float x[1000],y[1000];
  double x1,y1,x2,y2;
  int N=g1->GetN();
  for(int k=0;k<N;++k){
    g1->GetPoint(k,x1,y1);
    g2->GetPoint(k,x2,y2);
    x[k]=x1;
    y[k]=sqrt(y1*y1+y2*y2);
  }
  TGraph *g=new TGraph(N,x,y);
  g->SetTitle(g1->GetTitle());
  return(g);
}


TGraph *GquadSum(int n,TGraph **g1)
{
  float x[1000],y[1000];
  double x1,y1,x2,y2;
  int N=g1[0]->GetN();
  for(int k=0;k<N;++k){
    float sum=0;
    for(int i=0;i<n;++i){
      g1[i]->GetPoint(k,x1,y1);
      sum+=y1*y1;
      if(!i)x[k]=x1;
    }
    y[k]=sqrt(sum);
  }
  TGraph *g=new TGraph(N,x,y);
  g->SetTitle(g1[0]->GetTitle());
  return(g);
}


TGraphErrors *Gadd(TGraphErrors *g1,TGraphErrors *g2,int quad=0,int pm1=1,int pm2=1)
{
  float x[1000],y[1000],ex[1000],ey[1000];
  double x1,y1,e1,x2,y2,e2;
  int N=g1->GetN();
  for(int k=0;k<N;++k){
    g1->GetPoint(k,x1,y1);
    g2->GetPoint(k,x2,y2);
    e1=g1->GetErrorY(k);
    e2=g2->GetErrorY(k);
    x[k]=x1;
    ex[k]=g1->GetErrorX(k);
    if(quad){
      y[k]=sqrt(pm1*y1*y1+pm2*y2*y2);
      ey[k]=sqrt(pow(pm1*y1*e1,2)+pow(pm2*y2*e2,2))/y[k];
    }else{
      y[k]=pm1*y1+pm2*y2;
      ey[k]=sqrt(pow(pm1*e1,2)+pow(pm2*e2,2));
    }
  }
  TGraphErrors *g=new TGraphErrors(N,x,y,ex,ey);
  g->SetTitle(g1->GetTitle());
  return(g);
}


TGraphErrors *Gadd(TGraphErrors *g1,float f2,int quad=0,int pm1=1,int pm2=1)
{
  int N=g1->GetN();
  float xx[1000],ff[1000];for(int i=0;i<N;++i)ff[i]=f2;
  TGraphErrors *g2=new TGraphErrors(N,xx,ff,0,0);
  TGraphErrors *g=Gadd(g1,g2,quad,pm1,pm2);
  return(g);
}


TGraphErrors *Gdiff(TGraphErrors *g1,TGraphErrors *g2,int quad=0)
{
  TGraphErrors *g=Gadd(g1,g2,quad,1,-1);
  return(g);
}
/*
TGraphErrors *Gdiff(TGraphErrors *g1,TGraphErrors *g2,int quad=0)
{
  float x[1000],y[1000],ex[1000],ey[1000];
  double x1,y1,e1,x2,y2,e2;
  int N=g1->GetN();
  for(int k=0;k<N;++k){
    g1->GetPoint(k,x1,y1);
    g2->GetPoint(k,x2,y2);
    e1=g1->GetErrorY(k);
    e2=g2->GetErrorY(k);
    x[k]=x1;
    ex[k]=g1->GetErrorX(k);
    if(quad){
      y[k]=sqrt(y1*y1-y2*y2);
      ey[k]=sqrt(pow(y1*e1,2)+pow(y2*e2,2))/y[k];
    }else{
      y[k]=y1-y2;
      ey[k]=sqrt(e1*e1+e2*e2);
    }
  }
  TGraphErrors *g=new TGraphErrors(N,x,y,ex,ey);
  g->SetTitle(g1->GetTitle());
  return(g);
}
*/


TGraphErrors *Gprod(TGraphErrors *g1,TGraphErrors *g2)
{
  float x[1000],y[1000],ex[1000],ey[1000];
  double x1,y1,e1,x2,y2,e2;
  int N=g1->GetN();
  for(int k=0;k<N;++k){
    g1->GetPoint(k,x1,y1);
    g2->GetPoint(k,x2,y2);
    e1=g1->GetErrorY(k);
    e2=g2->GetErrorY(k);
    //printf("%f\t%f\t%f\t%f\t%f\t%f\n",x1,y1,e1,x2,y2,e2);
    x[k]=x1;
    ex[k]=g1->GetErrorX(k);
    if(y2!=0){
      y[k]=y1*y2;
      ey[k]=sqrt(pow(e1/y1,2)+pow(e2/y2,2))*fabs(y[k]);
    }
    else{y[k]=99999;ey[k]=0;}
  }
  TGraphErrors *g=new TGraphErrors(N,x,y,ex,ey);
  g->SetTitle(g1->GetTitle());
  return(g);
}


TGraphErrors *Gmult(TGraphErrors *g1,TGraphErrors *g2)
{
  float x[1000],y[1000],ex[1000],ey[1000];
  double x1,y1,e1,x2,y2,e2;
  int N=g1->GetN();
  for(int k=0;k<N;++k){
    g1->GetPoint(k,x1,y1);
    g2->GetPoint(k,x2,y2);
    e1=g1->GetErrorY(k);
    e2=g2->GetErrorY(k);
    //printf("%f\t%f\t%f\t%f\t%f\t%f\n",x1,y1,e1,x2,y2,e2);
    x[k]=x1;
    ex[k]=g1->GetErrorX(k);
    y[k]=y1*y2;
    ey[k]=sqrt(pow(e1*y2,2)+pow(e2*y1,2));
  }
  TGraphErrors *g=new TGraphErrors(N,x,y,ex,ey);
  g->SetTitle(g1->GetTitle());
  return(g);
}


TGraphErrors *Gdiv(TGraphErrors *g1,TGraphErrors *g2)
{
  float x[1000],y[1000],ex[1000],ey[1000];
  double x1,y1,e1,x2,y2,e2;
  int N=g1->GetN();
  for(int k=0;k<N;++k){
    g1->GetPoint(k,x1,y1);
    g2->GetPoint(k,x2,y2);
    e1=g1->GetErrorY(k);
    e2=g2->GetErrorY(k);
    //printf("%f\t%f\t%f\t%f\t%f\t%f\n",x1,y1,e1,x2,y2,e2);
    x[k]=x1;
    ex[k]=g1->GetErrorX(k);
    if(y2!=0){
      y[k]=y1/y2;
      ey[k]=sqrt(pow(e1/y1,2)+pow(e2/y2,2))*fabs(y[k]);
    }
    else{y[k]=99999;ey[k]=0;}
  }
  TGraphErrors *g=new TGraphErrors(N,x,y,ex,ey);
  g->SetTitle(g1->GetTitle());
  return(g);
}


TGraphErrors *Gdiv1(TGraphErrors *g1,TGraphErrors *g2)// g1/g2-1
{
  float x[1000],y[1000],ex[1000],ey[1000];
  double x1,y1,e1,x2,y2,e2;
  int N=g1->GetN();
  for(int k=0;k<N;++k){
    g1->GetPoint(k,x1,y1);
    g2->GetPoint(k,x2,y2);
    e1=g1->GetErrorY(k);
    e2=g2->GetErrorY(k);
    x[k]=x1;
    ex[k]=g1->GetErrorX(k);
    if(y2!=0){
      y[k]=y1/y2;
      ey[k]=sqrt(pow(e1/y1,2)+pow(e2/y2,2))*fabs(y[k]);
      y[k]-=1;
    }
    else{y[k]=99999;ey[k]=0;}
  }
  TGraphErrors *g=new TGraphErrors(N,x,y,ex,ey);
  g->SetTitle(g1->GetTitle());
  return(g);
}


TGraphErrors *Gsq(TGraphErrors *g1)
{
  float x[1000],y[1000],ex[1000],ey[1000];
  double x1,y1,e1;
  int N=g1->GetN();
  for(int k=0;k<N;++k){
    g1->GetPoint(k,x1,y1);
    e1=g1->GetErrorY(k);
    x[k]=x1;
    y[k]=y1*y1;
    ex[k]=g1->GetErrorX(k);
    ey[k]=2*y1*e1;
  }
  TGraphErrors *g=new TGraphErrors(N,x,y,ex,ey);
  g->SetTitle(g1->GetTitle());
  return(g);
}


TGraphErrors *Gsqrt(TGraphErrors *g1)
{
  float x[1000],y[1000],ex[1000],ey[1000];
  double x1,y1,e1;
  int N=g1->GetN();
  for(int k=0;k<N;++k){
    g1->GetPoint(k,x1,y1);
    e1=g1->GetErrorY(k);
    x[k]=x1;
    y[k]=sqrt(y1);
    ex[k]=g1->GetErrorX(k);
    if(y[k]==0)ey[k]=0;
    else ey[k]=e1/y[k]/2;
  }
  TGraphErrors *g=new TGraphErrors(N,x,y,ex,ey);
  g->SetTitle(g1->GetTitle());
  return(g);
}


double GcellX(TGraphErrors *g,int k)
{
  double x1,y1,e1;
  g->GetPoint(k,x1,y1);
  return(x1);
}


double GcellY(TGraphErrors *g,int k)
{
  double x1,y1,e1;
  g->GetPoint(k,x1,y1);
  return(y1);
}


void Gcell(TGraphErrors *g,int k,float *y,float *e)
{
  double x1,y1,e1;
  g->GetPoint(k,x1,y1);
  y[0]=y1;
  e[0]=g->GetErrorY(k);
  //printf("%f %f\t",y[0],e[0]);
  return;
}


void Gcell(TGraphErrors *g,int k,float *x,float *y,float *e)
{
  double x1,y1,e1;
  g->GetPoint(k,x1,y1);
  x[0]=x1;
  y[0]=y1;
  e[0]=g->GetErrorY(k);
  //printf("%f %f\t",y[0],e[0]);
  return;
}


void Glimits(float *xmin,float *xmax,float *ymin,float *ymax,TGraphErrors *g0=0,TGraphErrors *g1=0,TGraphErrors *g2=0,TGraphErrors *g3=0,TGraphErrors *g4=0,TGraphErrors *g5=0,TGraphErrors *g6=0,TGraphErrors *g7=0,TGraphErrors *g8=0,TGraphErrors *g9=0)
{
  Glimits(0,0.2,0.2,0.15,0.15,xmin,xmax,ymin,ymax,g0,g1,g2,g3,g4,g5,g6,g7,g8,g9);
}


void GlimitsGraph(float *xmin,float *xmax,float *ymin,float *ymax,TGraph *g0=0,TGraph *g1=0,TGraph *g2=0,TGraph *g3=0,TGraph *g4=0,TGraph *g5=0,TGraph *g6=0,TGraph *g7=0,TGraph *g8=0,TGraph *g9=0)
{
  GlimitsGraph(0,0.2,0.2,0.15,0.15,xmin,xmax,ymin,ymax,g0,g1,g2,g3,g4,g5,g6,g7,g8,g9);
}


void Glimits(float dxmin,float dxmax,float dymin,float dymax,float *xmin,float *xmax,float *ymin,float *ymax,TGraphErrors *g0=0,TGraphErrors *g1=0,TGraphErrors *g2=0,TGraphErrors *g3=0,TGraphErrors *g4=0,TGraphErrors *g5=0,TGraphErrors *g6=0,TGraphErrors *g7=0,TGraphErrors *g8=0,TGraphErrors *g9=0)
{
  Glimits(0,dxmin,dxmax,dymin,dymax,xmin,xmax,ymin,ymax,g0,g1,g2,g3,g4,g5,g6,g7,g8,g9);
}

void GlimitsGraph(float dxmin,float dxmax,float dymin,float dymax,float *xmin,float *xmax,float *ymin,float *ymax,TGraph *g0=0,TGraph *g1=0,TGraph *g2=0,TGraph *g3=0,TGraph *g4=0,TGraph *g5=0,TGraph *g6=0,TGraph *g7=0,TGraph *g8=0,TGraph *g9=0)
{
  GlimitsGraph(0,dxmin,dxmax,dymin,dymax,xmin,xmax,ymin,ymax,g0,g1,g2,g3,g4,g5,g6,g7,g8,g9);
}

void Glimits(int remember,float *xmin,float *xmax,float *ymin,float *ymax,TGraphErrors *g0=0,TGraphErrors *g1=0,TGraphErrors *g2=0,TGraphErrors *g3=0,TGraphErrors *g4=0,TGraphErrors *g5=0,TGraphErrors *g6=0,TGraphErrors *g7=0,TGraphErrors *g8=0,TGraphErrors *g9=0)
{
  Glimits(remember,0.2,0.2,0.15,0.15,xmin,xmax,ymin,ymax,g0,g1,g2,g3,g4,g5,g6,g7,g8,g9);
}


void GlimitsGraph(int remember,float *xmin,float *xmax,float *ymin,float *ymax,TGraph *g0=0,TGraph *g1=0,TGraph *g2=0,TGraph *g3=0,TGraph *g4=0,TGraph *g5=0,TGraph *g6=0,TGraph *g7=0,TGraph *g8=0,TGraph *g9=0)
{
  GlimitsGraph(remember,0.2,0.2,0.15,0.15,xmin,xmax,ymin,ymax,g0,g1,g2,g3,g4,g5,g6,g7,g8,g9);
}


void Glimits(int remember,float dxmin,float dxmax,float dymin,float dymax,float *xmin,float *xmax,float *ymin,float *ymax,TGraphErrors *g0=0,TGraphErrors *g1=0,TGraphErrors *g2=0,TGraphErrors *g3=0,TGraphErrors *g4=0,TGraphErrors *g5=0,TGraphErrors *g6=0,TGraphErrors *g7=0,TGraphErrors *g8=0,TGraphErrors *g9=0)
{
  const int Ng=10;
  TGraphErrors *g[Ng];
  g[0]=g0;g[1]=g1;g[2]=g2;g[3]=g3;g[4]=g4;g[5]=g5;g[6]=g6;g[7]=g7;g[8]=g8;g[9]=g9;
  double x,y,ex,ey;
  if(!remember){
    xmin[0]=999999;xmax[0]=-999999;
    ymin[0]=999999;ymax[0]=-999999;
  }
  for(int i=0;i<Ng;++i){
    if(!g[i])continue;
    for(int k=0;k<g[i]->GetN();++k){
      g[i]->GetPoint(k,x,y);
      ex=g[i]->GetErrorX(k);
      ey=g[i]->GetErrorY(k);
      //printf("[%d]\t%f %f\t%f %f\n",k,y,ey,y-ey,y+ey);
      if(x-ex<xmin[0])xmin[0]=x-ex;
      if(x+ex>xmax[0])xmax[0]=x+ex;
      if(y-ey<ymin[0])ymin[0]=y-ey;
      if(y+ey>ymax[0])ymax[0]=y+ey;
    }
  }
  float dx=xmax[0]-xmin[0];xmin[0]-=dxmin*dx;xmax[0]+=dxmax*dx;
  float dy=ymax[0]-ymin[0];ymin[0]-=dymin*dy;ymax[0]+=dymax*dy;
  return;
}


void GlimitsGraph(int remember,float dxmin,float dxmax,float dymin,float dymax,float *xmin,float *xmax,float *ymin,float *ymax,TGraph *g0=0,TGraph *g1=0,TGraph *g2=0,TGraph *g3=0,TGraph *g4=0,TGraph *g5=0,TGraph *g6=0,TGraph *g7=0,TGraph *g8=0,TGraph *g9=0)
{
  const int Ng=10;
  TGraph *g[Ng];
  g[0]=g0;g[1]=g1;g[2]=g2;g[3]=g3;g[4]=g4;g[5]=g5;g[6]=g6;g[7]=g7;g[8]=g8;g[9]=g9;
  double x,y;
  if(!remember){
    xmin[0]=999999;xmax[0]=-999999;
    ymin[0]=999999;ymax[0]=-999999;
  }
  for(int i=0;i<Ng;++i){
    if(!g[i])continue;
    for(int k=0;k<g[i]->GetN();++k){
      g[i]->GetPoint(k,x,y);
      if(x<xmin[0])xmin[0]=x;
      if(x>xmax[0])xmax[0]=x;
      if(y<ymin[0])ymin[0]=y;
      if(y>ymax[0])ymax[0]=y;
    }
  }
  float dx=xmax[0]-xmin[0];xmin[0]-=dxmin*dx;xmax[0]+=dxmax*dx;
  float dy=ymax[0]-ymin[0];ymin[0]-=dymin*dy;ymax[0]+=dymax*dy;
  return;
}


void Glimits2(float *xmin,float *xmax,float *ymin,float *ymax,float xmin0,float xmax0,float ymin0,float ymax0,TGraphErrors *g,int i0=-1,int i1=-1,int i2=-1,int i3=-1,int i4=-1,int i5=-1,int i6=-1,int i7=-1)
{
  const int Ni=8;
  int idx[Ni];
  idx[0]=i0;idx[1]=i1;idx[2]=i2;idx[3]=i3;idx[4]=i4;idx[5]=i5;idx[6]=i6;idx[7]=i7;
  double x,y,ex,ey;
  xmin[0]=xmin0;xmax[0]=xmax0;
  ymin[0]=ymin0;ymax[0]=ymax0;
  for(int i=0;i<Ni;++i){
    if(idx[i]<0)continue;
    g->GetPoint(idx[i],x,y);
    ex=g->GetErrorX(idx[i]);
    ey=g->GetErrorY(idx[i]);
    if(x-ex<xmin[0])xmin[0]=x-ex;
    if(x+ex>xmax[0])xmax[0]=x+ex;
    if(y-ey<ymin[0])ymin[0]=y-ey;
    if(y+ey>ymax[0])ymax[0]=y+ey;
  }
  float dx=0.20*(xmax[0]-xmin[0]);xmin[0]-=dx;xmax[0]+=dx;
  float dy=0.15*(ymax[0]-ymin[0]);ymin[0]-=dy;ymax[0]+=dy;
  return;
}


void Glimits3(float *ymin,float *ymax,float ymin0,float ymax0,float y0=1e6,float y1=1e6,float y2=1e6,float y3=1e6,float y4=1e6,float y5=1e6,float y6=1e6,float y7=1e6)
{
  const int Ny=8;
  float y[Ny]={y0,y1,y2,y3,y4,y5,y6,y7};
  ymin[0]=ymin0;ymax[0]=ymax0;
  for(int i=0;i<Ny;++i){
    if(y[i]>999999)continue;
    if(y[i]<ymin[0])ymin[0]=y[i];
    if(y[i]>ymax[0])ymax[0]=y[i];
  }
  float dy=0.10*(ymax[0]-ymin[0]);ymin[0]-=dy;ymax[0]+=dy;
  return;
}


TGraphErrors *Gcomb(TGraphErrors *gAA,TGraphErrors *gdAu,int i1=0,int i2=0)
{
  float x[1000],y[1000],ex[1000],ey[1000];
  double x1,y1,e1,x2,y2,e2;
  if(!i1&&!i2){TGraphErrors *g=gAA->Clone();return(g);}
  if(!i1&&i2){i1=i2;i2=0;}
  gdAu->GetPoint(i1,x1,y1);
  e1=gdAu->GetErrorY(i1);
  if(i2){
    gdAu->GetPoint(i2,x2,y2);
    e2=gdAu->GetErrorY(i2);
    y1=(y1+y2)/2;
    e1=sqrt(e1*e1+e2*e2)/2;
  }
  int N=1+gAA->GetN();
  x[0]=x1;
  y[0]=y1;
  ex[0]=gdAu->GetErrorX(i1);
  ey[0]=e1;
  for(int k=0;k<N;++k){
    gAA->GetPoint(k,x1,y1);
    x[k+1]=x1;
    y[k+1]=y1;
    ex[k+1]=gAA->GetErrorX(k);
    ey[k+1]=gAA->GetErrorY(k);
  }
  TGraphErrors *g=new TGraphErrors(N,x,y,ex,ey);
  g->SetTitle(gAA->GetTitle());
  return(g);
}

TGraphErrors *Gspread(int n,TGraphErrors **g)
{
  TGraphErrors *gave=Gave(n,g);
  TGraphErrors *gspread;
  for(int i=0;i<n;++i){
    TGraphErrors *gd=Gdiff(g[i],gave);
    gd=Gsq(gd);
    if(!i)gspread=gd;
    else gspread=Gadd(gspread,gd);
  }
  gspread=Gsqrt(gspread);
  gspread=Gscale(gspread,1./n);
  return(gspread);
}

TGraphErrors *Gspread(TGraphErrors *g0=0,TGraphErrors *g1=0,TGraphErrors *g2=0,TGraphErrors *g3=0,TGraphErrors *g4=0,TGraphErrors *g5=0,TGraphErrors *g6=0,TGraphErrors *g7=0,TGraphErrors *g8=0,TGraphErrors *g9=0)
{
  int n=0;
  TGraphErrors *g[10];
  if(g0){g[n]=g0;n++;}
  if(g1){g[n]=g1;n++;}
  if(g2){g[n]=g2;n++;}
  if(g3){g[n]=g3;n++;}
  if(g4){g[n]=g4;n++;}
  if(g5){g[n]=g5;n++;}
  if(g6){g[n]=g6;n++;}
  if(g7){g[n]=g7;n++;}
  if(g8){g[n]=g8;n++;}
  if(g9){g[n]=g9;n++;}
  TGraphErrors *gspread=Gspread(n,g);
  return(gspread);
}

TGraphErrors *GdiffHalf(TGraphErrors *g1,TGraphErrors *g2)
{
  TGraphErrors *g=Gdiff(g1,g2);
  g=Gabs(g);
  g=Gscale(g,0.5);
  return(g);
}

TGraph *Gdmax(TGraphErrors *g0,int n,TGraphErrors **g,int pm=0)
{
  const int N=g0->GetN();
  if(N>1000){printf("Gdmax error!\n");return(0);}
  float x[1000]={0},y[1000]={0};
  TGraph *gdmax=new TGraph(N,x,y);
  for(int i=0;i<n;++i){
    TGraphErrors *gd=Gdiff(g[i],g0);
    if(pm==1)gdmax=GmaxPlus(gdmax,gd);
    else if(pm==-1)gdmax=GmaxMinus(gdmax,gd);
    else gdmax=Gabsmax(gdmax,gd);
  }
  return(gdmax);
}

TGraphErrors *Gdmax(TGraphErrors *gg=0,TGraphErrors *g0=0,TGraphErrors *g1=0,TGraphErrors *g2=0,TGraphErrors *g3=0,TGraphErrors *g4=0,TGraphErrors *g5=0,TGraphErrors *g6=0,TGraphErrors *g7=0,TGraphErrors *g8=0,TGraphErrors *g9=0)
{
  int n=0;
  TGraphErrors *g[10];
  if(g0){g[n]=g0;n++;}
  if(g1){g[n]=g1;n++;}
  if(g2){g[n]=g2;n++;}
  if(g3){g[n]=g3;n++;}
  if(g4){g[n]=g4;n++;}
  if(g5){g[n]=g5;n++;}
  if(g6){g[n]=g6;n++;}
  if(g7){g[n]=g7;n++;}
  if(g8){g[n]=g8;n++;}
  if(g9){g[n]=g9;n++;}
  TGraphErrors *gdmax=Gdmax(gg,n,g);
  return(gdmax);
}

TGraphErrors *GdmaxPlusMinus(TGraphErrors *gg=0,int plusminus=0,TGraphErrors *g0=0,TGraphErrors *g1=0,TGraphErrors *g2=0,TGraphErrors *g3=0,TGraphErrors *g4=0,TGraphErrors *g5=0,TGraphErrors *g6=0,TGraphErrors *g7=0,TGraphErrors *g8=0,TGraphErrors *g9=0)
{
  int n=0;
  TGraphErrors *g[10];
  if(g0){g[n]=g0;n++;}
  if(g1){g[n]=g1;n++;}
  if(g2){g[n]=g2;n++;}
  if(g3){g[n]=g3;n++;}
  if(g4){g[n]=g4;n++;}
  if(g5){g[n]=g5;n++;}
  if(g6){g[n]=g6;n++;}
  if(g7){g[n]=g7;n++;}
  if(g8){g[n]=g8;n++;}
  if(g9){g[n]=g9;n++;}
  TGraphErrors *gdmax=Gdmax(gg,n,g,plusminus);
  return(gdmax);
}

/*
void Gwrite(char *file,TGraphErrors *g)
{
  FILE *fp=fopen(file,"w");
  double x,y,ex,ey;
  int N=g->GetN();
  fprintf(fp,"%d\n",N);
  for(int k=0;k<N;++k){
    g->GetPoint(k,x,y);
    ex=g->GetErrorX(k);
    ey=g->GetErrorY(k);
    fprintf(fp,"%d\t%f\t%f\t%f\t%f\n",k,x,ex,y,ey);
  }
  fclose(fp);
  return;
}
*/

void Gwrite(char *file,TGraphErrors *g0,TGraphErrors *g1=0,TGraphErrors *g2=0,TGraphErrors *g3=0,TGraphErrors *g4=0,TGraphErrors *g5=0,TGraphErrors *g6=0,TGraphErrors *g7=0,TGraphErrors *g8=0,TGraphErrors *g9=0)
{
  int n=0;
  TGraphErrors *g[10];
  if(g0){g[n]=g0;n++;}
  if(g1){g[n]=g1;n++;}
  if(g2){g[n]=g2;n++;}
  if(g3){g[n]=g3;n++;}
  if(g4){g[n]=g4;n++;}
  if(g5){g[n]=g5;n++;}
  if(g6){g[n]=g6;n++;}
  if(g7){g[n]=g7;n++;}
  if(g8){g[n]=g8;n++;}
  if(g9){g[n]=g9;n++;}
  Gwrite(file,n,g);
}

void Gwrite(char *file,int n=1,TGraphErrors **g)
{
  printf("writing %s\n",file);
  FILE *fp=fopen(file,"w");
  double x,y,ex,ey;
  for(int i=0;i<n;++i){
    int N=g[i]->GetN();
    fprintf(fp,"// %d %s\n",N,g[i]->GetTitle());
    for(int k=0;k<N;++k){
      g[i]->GetPoint(k,x,y);
      ex=g[i]->GetErrorX(k);
      ey=g[i]->GetErrorY(k);
      fprintf(fp,"%d\t%f\t%f\t%f\t%f\n",k,x,ex,y,ey);
    }
  }
  fclose(fp);
  return;
}

void Gwrite(char *file,TGraph *g0,TGraph *g1=0,TGraph *g2=0,TGraph *g3=0,TGraph *g4=0,TGraph *g5=0,TGraph *g6=0,TGraph *g7=0,TGraph *g8=0,TGraph *g9=0)
{
  int n=0;
  TGraph *g[10];
  if(g0){g[n]=g0;n++;}
  if(g1){g[n]=g1;n++;}
  if(g2){g[n]=g2;n++;}
  if(g3){g[n]=g3;n++;}
  if(g4){g[n]=g4;n++;}
  if(g5){g[n]=g5;n++;}
  if(g6){g[n]=g6;n++;}
  if(g7){g[n]=g7;n++;}
  if(g8){g[n]=g8;n++;}
  if(g9){g[n]=g9;n++;}
  Gwrite(file,n,g);
}

void Gwrite(char *file,int n=1,TGraph **g)
{
  printf("writing %s\n",file);
  FILE *fp=fopen(file,"w");
  double x,y;
  for(int i=0;i<n;++i){
    int N=g[i]->GetN();
    fprintf(fp,"// %d %s\n",N,g[i]->GetTitle());
    for(int k=0;k<N;++k){
      g[i]->GetPoint(k,x,y);
      fprintf(fp,"%d\t%f\t%f\n",k,x,y);
    }
  }
  fclose(fp);
  return;
}


void Gwrite2(char *file,int n=1,TGraph **g1,TGraph **g2)
{
  printf("writing %s\n",file);
  FILE *fp=fopen(file,"w");
  double x1,y1,x2,y2;
  for(int i=0;i<n;++i){
    int N=g1[i]->GetN();
    if(g2[i]->GetN()!=N){printf("Gwrite2 error!\n");return;}
    fprintf(fp,"// %d %s %s\n",N,g1[i]->GetTitle(),g2[i]->GetTitle());
    for(int k=0;k<N;++k){
      g1[i]->GetPoint(k,x1,y1);
      g2[i]->GetPoint(k,x2,y2);
      fprintf(fp,"%d\t%f\t%f\t%f\t%f\n",k,x1,y1,x2,y2);
    }
  }
  fclose(fp);
  return;
}


void graphSyst(TGraphErrors *g,TGraphErrors *g1,int fcolor=5,int bcolor=0,float dx=0,float dy=0,int ierrX=1){
  TGraphErrors *g2=Gscale(g1,-1);
  graphSyst(g,g1,g2,fcolor,bcolor,dx,dy,1,ierrX);
}


void graphSyst(TGraphErrors *g,float syst=0,int fcolor=5,int bcolor=0,float dx=0,float dy=0,int ierrX=1){
  TGraphErrors *g1=Gscale(g,syst);
  TGraphErrors *g2=Gscale(g,-syst);
  graphSyst(g,g1,g2,fcolor,bcolor,dx,dy,1,ierrX);
}


void graphSyst(TGraphErrors *g,TGraphErrors *g1,TGraphErrors *g2,int fcolor=5,int bcolor=0,float dx=0,float dy=0,int iadd=0,int ierrX=1){
  //fcolor!=0 draw box and fill with this color
  //bcolor!=0 draw box with this color
  //dx is used to draw box
  //dy!=0 draw caps
  //g1 or g2->errorX!=0 draw histogram-style boundary
  TLine *line=new TLine();
  float grafx[5],grafy[5];
  double x,y,y1,y2,ddx,ddy,y1old,y2old;
  int n=g->GetN();
  for(int i=0;i<n;++i){
    g1->GetPoint(i,x,y1);
    g2->GetPoint(i,x,y2);
    g->GetPoint(i,x,y);
    if(iadd){y1+=y;y2+=y;}
    ddx=g->GetErrorX(i);
    if(!ierrX)g->SetPointError(i,0,g->GetErrorY(i));
    if(ddx==0)ddx=dx;
    if(fcolor||bcolor){
      grafx[0]=x-ddx;grafy[0]=y1;
      grafx[1]=x+ddx;grafy[1]=y1;
      grafx[2]=x+ddx;grafy[2]=y2;
      grafx[3]=x-ddx;grafy[3]=y2;
      grafx[4]=x-ddx;grafy[4]=y1;
      TGraph *graf=new TGraph(5,grafx,grafy);
      if(fcolor){graf->SetFillColor(fcolor);graf->Draw("f");}
      if(bcolor){graf->SetLineWidth(0.5*(g1->GetLineWidth()+g2->GetLineWidth()));
	graf->SetLineColor(bcolor);graf->Draw("l");}
    }
    //if(g1->GetLineColor()){
    if(dy!=0){
      line->SetLineWidth(g1->GetLineWidth());
      line->SetLineColor(g1->GetLineColor());
      line->DrawLine(x-ddx,y1,x+ddx,y1);
      if(y1<y)ddy=dy;else ddy=-dy;
      line->DrawLine(x-ddx,y1,x-ddx,y1+ddy);
      line->DrawLine(x+ddx,y1,x+ddx,y1+ddy);
      //}
      //if(g2->GetLineColor()){
      line->SetLineWidth(g2->GetLineWidth());
      line->SetLineColor(g2->GetLineColor());
      line->DrawLine(x-ddx,y2,x+ddx,y2);
      if(y2<y)ddy=dy;else ddy=-dy;
      line->DrawLine(x-ddx,y2,x-ddx,y2+ddy);
      line->DrawLine(x+ddx,y2,x+ddx,y2+ddy);
    }
    //if(!i){y1old=y;y2old=y;}
    ddx=g1->GetErrorX(i);
    if(ddx!=0){
      line->SetLineWidth(g1->GetLineWidth());
      line->SetLineColor(g1->GetLineColor());
      if(i)line->DrawLine(x-ddx,y1,x-ddx,y1old);
      line->DrawLine(x-ddx,y1,x+ddx,y1);y1old=y1;
    }
    ddx=g2->GetErrorX(i);
    if(ddx!=0){
      line->SetLineWidth(g2->GetLineWidth());
      line->SetLineColor(g2->GetLineColor());
      if(i)line->DrawLine(x-ddx,y2,x-ddx,y2old);
      line->DrawLine(x-ddx,y2,x+ddx,y2);y2old=y2;
    }
  }
  if(g->GetMarkerStyle())g->Draw("p");
}


/*
void graphSyst(TGraphErrors *g,TGraphErrors *g1,TGraphErrors *g2,int fcolor=5,int bcolor=0,float dx=0,float dy=0){
  //fcolor!=0 fill with this color
  //bcolor!=0 draw box with this color
  //g1 or g2->color!=0 draw cap with its color and width
  //g1 or g2->error!=0 draw histogram-style boundary
  TLine *line=new TLine();
  float grafx[5],grafy[5];
  double x,y,y1,y2,ddx,ddy,y1old,y2old;
  int n=g->GetN();
  for(int i=0;i<n;++i){
    g1->GetPoint(i,x,y1);
    g2->GetPoint(i,x,y2);
    g->GetPoint(i,x,y);
    y1+=y;y2+=y;
    if(fcolor||bcolor){
      grafx[0]=x-dx;grafy[0]=y1;
      grafx[1]=x+dx;grafy[1]=y1;
      grafx[2]=x+dx;grafy[2]=y2;
      grafx[3]=x-dx;grafy[3]=y2;
      grafx[4]=x-dx;grafy[4]=y1;
      TGraph *graf=new TGraph(5,grafx,grafy);
      if(fcolor){graf->SetFillColor(fcolor);graf->Draw("f");}
      if(bcolor){graf->SetLineWidth(0.5*(g1->GetLineWidth()+g2->GetLineWidth()));
	graf->SetLineColor(bcolor);graf->Draw("l");}
    }
    if(g1->GetLineColor()){
      line->SetLineWidth(g1->GetLineWidth());
      line->SetLineColor(g1->GetLineColor());
      line->DrawLine(x-dx,y1,x+dx,y1);
      if(y1<y)ddy=dy;else ddy=-dy;
      line->DrawLine(x-dx,y1,x-dx,y1+ddy);
      line->DrawLine(x+dx,y1,x+dx,y1+ddy);
    }
    if(g2->GetLineColor()){
      line->SetLineWidth(g2->GetLineWidth());
      line->SetLineColor(g2->GetLineColor());
      line->DrawLine(x-dx,y2,x+dx,y2);
      if(y2<y)ddy=dy;else ddy=-dy;
      line->DrawLine(x-dx,y2,x-dx,y2+ddy);
      line->DrawLine(x+dx,y2,x+dx,y2+ddy);
    }
    if(!i){y1old=y;y2old=y;}
    ddx=g1->GetErrorX(i);
    if(ddx!=0){
      line->SetLineWidth(g1->GetLineWidth());
      line->SetLineColor(g1->GetLineColor());
      //line->DrawLine(x-ddx,y1,x-ddx,y1old);
      line->DrawLine(x-ddx,y1,x+ddx,y1);y1old=y1;
    }
    ddx=g2->GetErrorX(i);
    if(ddx!=0){
      line->SetLineWidth(g2->GetLineWidth());
      line->SetLineColor(g2->GetLineColor());
      //line->DrawLine(x-ddx,y2,x-ddx,y2old);
      line->DrawLine(x-ddx,y2,x+ddx,y2);y2old=y2;
    }
  }
  if(g->GetMarkerStyle())g->Draw("p");
}
*/

void GplotSyst(float y0,TGraph *gsyst,int fillcol=0,int lsty=1,int lwid=1,int top=0)
{
  float x[1000],y[1000];
  double x1,y1;
  int N=gsyst->GetN();
  for(int k=0;k<N;++k){
    gsyst->GetPoint(k,x1,y1);
    x[k]=x1;
    if(top)y[k]=y0-y1;
    else y[k]=y0+y1;
  }
  x[N]=x[N-1],y[N]=y0;
  x[N+1]=x[0],y[N+1]=y0;
  TGraph *g=new TGraph(N+2,x,y);
  if(fillcol<0){
    g->SetLineColor(-fillcol);
    g->SetLineStyle(lsty);
    g->SetLineWidth(lwid);
    g->Draw("l");
  }else{
    g->SetFillColor(fillcol);
    g->Draw("f");
  }
}

void GplotSyst(float y0,TGraph *gsyst1,TGraph *gsyst2,int fillcol=0,int lsty=1,int lwid=1,int top=0)
{
  TGraph *gsyst=Gmax(gsyst1,gsyst2);
  GplotSyst(y0,gsyst,fillcol,lsty,lwid,top);
}

void GplotSyst(TGraphErrors *g0,TGraph *gsyst,int fillcol=0,int lsty=1,int lwid=1,int top=0)
{
  float x[1000],y[1000];
  double x0,y0,dy,xdum;
  int N=g0->GetN();
  for(int k=0;k<N;++k){
    g0->GetPoint(k,x0,y0);
    gsyst->GetPoint(k,xdum,dy);
    x[k]=x0;
    y[k]=y0-dy;
    x[2*N-1-k]=x0;
    y[2*N-1-k]=y0+dy;
  }
  TGraph *g=new TGraph(2*N,x,y);
  if(fillcol<0){
    g->SetLineColor(-fillcol);
    g->SetLineStyle(lsty);
    g->SetLineWidth(lwid);
    g->Draw("l");
  }else{
    g->SetFillColor(fillcol);
    g->Draw("f");
  }
}


TGraph *Gmax(TGraph *g1,TGraph *g2)
{
  float x[1000],y[1000];
  double x1,y1,x2,y2;
  int N=g1->GetN();
  for(int k=0;k<N;++k){
    g1->GetPoint(k,x1,y1);
    g2->GetPoint(k,x2,y2);
    x[k]=x1;
    if(y1>y2)y[k]=y1;
    else y[k]=y2;
  }
  TGraph *g=new TGraph(N,x,y);
  g->SetTitle(g1->GetTitle());
  return(g);
}


TGraph *Gabsmax(TGraph *g1,TGraphErrors *g2)
{
  float x[1000],y[1000],ex[1000],ey[1000];
  double x1,y1,e1,x2,y2,e2;
  int N=g1->GetN();
  for(int k=0;k<N;++k){
    g1->GetPoint(k,x1,y1);
    g2->GetPoint(k,x2,y2);
    y1=fabs(y1);
    y2=fabs(y2);
    x[k]=x2;
    if(y1>y2)y[k]=y1;
    else y[k]=y2;
  }
  TGraph *g=new TGraph(N,x,y);
  g->SetTitle(g1->GetTitle());
  return(g);
}

TGraph *GmaxPlus(TGraph *g1,TGraph *g2)
{
  float x[1000],y[1000];
  double x1,y1,x2,y2;
  int N=g1->GetN();
  for(int k=0;k<N;++k){
    g1->GetPoint(k,x1,y1);
    g2->GetPoint(k,x2,y2);
    x[k]=x2;
    if(y1<0)y1=0;
    if(y2<0)y2=0;
    if(y1>y2)y[k]=y1;
    else y[k]=y2;
  }
  TGraph *g=new TGraph(N,x,y);
  g->SetTitle(g1->GetTitle());
  return(g);
}


TGraph *GmaxMinus(TGraph *g1,TGraph *g2)
{
  float x[1000],y[1000];
  double x1,y1,x2,y2;
  int N=g1->GetN();
  for(int k=0;k<N;++k){
    g1->GetPoint(k,x1,y1);
    g2->GetPoint(k,x2,y2);
    x[k]=x2;
    if(y1>0)y1=0;
    if(y2>0)y2=0;
    if(y1>y2)y[k]=y2;
    else y[k]=y1;
  }
  TGraph *g=new TGraph(N,x,y);
  g->SetTitle(g1->GetTitle());
  return(g);
}


