char *fileName(char *filename=""){
  TDatime *t=new TDatime();
  char *file=malloc(256*sizeof(char));
  sprintf(file,"%s[%d@%d]",filename,t->GetDate(),t->GetTime());
  printf("%s\n",file);
  return(file);
}

void stamp(char *string=""){
  TLatex *ltx=new TLatex();
  ltx->SetNDC(1);
  ltx->SetTextColor(17);
  ltx->SetTextSize(0.025);
  ltx->SetTextAlign(13);
  //ltx->DrawLatex(0.01,0.99,Form("%s/%s",gSystem->pwd(),gInterpreter->GetCurrentMacroName()));
  ltx->DrawLatex(0.01,0.99,Form("%s %s",gInterpreter->GetTopLevelMacroName(),string));
}

int Quit(char* prompt="Hit q or n to quit...")
{
  char character;
  fflush(stdin);
  printf("%s",prompt);
  fflush(stdout);
  character = getc(stdin);
  printf("\n");
  return ((character=='q')||
	  (character=='Q')||
	  (character=='n')||
	  (character=='N'));
}


/////////////////////////////////////////////////////////

char *StringReplace(char *str, char *orig, char *rep)
{
  char buffer[4096];
  char *p;
  if(!(p = strstr(str, orig)))  // Is 'orig' even in 'str'?
    return str;
  strncpy(buffer, str, p-str); // Copy characters from 'str' start to 'orig' st$
  buffer[p-str] = '\0';
  sprintf(buffer+(p-str), "%s%s", rep, p+strlen(orig));
  return buffer;
}

void StringCuts(char *o, const char *s, const char *s1)
{
  //o = malloc(1000);
  strcpy(o,s);
  strcat(o, " {"); strcat(o,s1); strcat(o,"}");
}

void StringCuts(char *o, const char *s, const char *s1, const char *s2)
{
  strcpy(o,s);
  strcat(o, " {");
  strcat(o, " ("); strcat(o,s1); strcat(o,")");
  strcat(o,"&&("); strcat(o,s2); strcat(o,")");
  strcat(o,  "}");
}

void StringCuts(char *o, const char *s, const char *s1, 
		const char *s2, const char *s3)
{
  strcpy(o,s);
  strcat(o, " {");
  strcat(o,  "("); strcat(o,s1); strcat(o,")");
  strcat(o,"&&("); strcat(o,s2); strcat(o,")");
  strcat(o,"&&("); strcat(o,s3); strcat(o,")");
  strcat(o,  "}");
}

void StringCuts(char *o, const char *s, const char *s1, const char *s2, 
		const char *s3, const char *s4)
{
  strcpy(o,s);
  strcat(o, " {");
  strcat(o,  "("); strcat(o,s1); strcat(o,")");
  strcat(o,"&&("); strcat(o,s2); strcat(o,")");
  strcat(o,"&&("); strcat(o,s3); strcat(o,")");
  strcat(o,"&&("); strcat(o,s4); strcat(o,")");
  strcat(o,  "}");
}


/////////////////////////////////////////////////////////

void Style(int i) 
{
  gStyle->Reset();

  //printf("\nStyle %d ...\n\n",i);
  switch (i) {
  case 1:  //talk style
    //gStyle->SetCanvasColor(10);              // white
    //gStyle->SetPadColor(10);                 // white
    gStyle->SetCanvasColor(0);              // white
    gStyle->SetPadColor(0);                 // white
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    //gStyle->SetFillColor(0);                 // clear (no-fill)
    gStyle->SetTitleColor(0);                // clear (no-fill)
    gStyle->SetStatColor(0);                 // clear (no-fill)
    gStyle->SetHistFillColor(0);             // clear (no-fill)
    //gStyle->SetFrameFillColor(10);           // white
    gStyle->SetFrameFillColor(0);           // white
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit();
    gStyle->SetLabelSize(0.06,"X");
    gStyle->SetLabelSize(0.06,"Y");
    gStyle->SetLineWidth(2);
    gStyle->SetTitleBorderSize(1);
    break;
  case 2:  //paper figure style
    //gStyle->SetCanvasColor(10);              // white
    //gStyle->SetPadColor(10);                 // white
    gStyle->SetCanvasColor(0);              // white
    gStyle->SetPadColor(0);                 // white
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetFillColor(0);                 // clear (no-fill)
    gStyle->SetTitleColor(0);                // clear (no-fill)
    gStyle->SetStatColor(0);                 // clear (no-fill)
    gStyle->SetHistFillColor(0);             // clear (no-fill)
    //gStyle->SetFrameFillColor(10);           // white
    gStyle->SetFrameFillColor(0);           // white
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit();
    gStyle->SetLabelSize(0.05,"X");
    gStyle->SetLabelSize(0.05,"Y");
    gStyle->SetLineWidth(1);
    gStyle->SetTitleBorderSize(1);
    break;
  default:
    gStyle->Reset();
    break;
  }
}

/////////////////////////////////////////////////////////
void XYpads(Float_t leftOffset, Float_t bottomOffset, const Int_t nx, const Int_t ny, 
	    Float_t left, Float_t right, Float_t top, Float_t bottom, Float_t xgap, Float_t ygap,
	    TPad** XYpad, Float_t* xyCent,
	    Float_t* XpadSize, Float_t* XpadLeft, Float_t* XpadRite,
	    Float_t* YpadSize, Float_t* YpadTop, Float_t* YpadBot)
//all variables are in parent frame, except XpadLeft, XpadRight, YpadTop, and YpadBot are in local frame of pad.
{
  float xsize=(1-leftOffset-left-right-(nx-1)*xgap)/nx;
  float ysize=(1-top-bottom-bottomOffset-(ny-1)*ygap)/ny;
  float rite=xgap/2;
  float x1=leftOffset,x2=leftOffset;
  for(int i=0;i<nx;++i){
    if(i)left=xgap/2;
    x2+=left;
    //XpadCent[i]=x2+xsize/2;//this is in parent frame
    x2+=xsize;
    if(i==nx-1)rite=right;
    x2+=rite;
    if(x2>1)x2=1;
    XpadSize[i]=x2-x1;
    XpadLeft[i]=left/XpadSize[i];//this is in pad frame
    XpadRite[i]=1-rite/XpadSize[i];//this is in pad frame
    //XpadCent[i]=(1-rite/XpadSize[i]+left/XpadSize[i])/2;//this is in pad frame
    float tp=top,bot=ygap/2;
    float y1=1,y2=1;
    for(int j=0;j<ny;++j){
      if(j)tp=ygap/2;
      y1-=tp;
      //YpadCent[i]=y1-ysize/2;//this is in parent frame
      y1-=ysize;
      if(j==ny-1)bot=bottom;
      y1-=bot;
      YpadSize[j]=y2-y1;
      YpadTop[i]=1-tp/YpadSize[i];//this is in pad frame
      YpadBot[i]=bot/YpadSize[i];//this is in pad frame
      //YpadCent[i]=(1-tp/YpadSize[i]+bot/YpadSize[i])/2;//this is in pad frame
      //XYpad[i][j]=new TPad(Form("pad%d_%d",i,j),"",x1,y1,x2,y2);
      //XYpad[i][j]->SetLeftMargin(left/XpadSize[i]);XYpad[i][j]->SetRightMargin(rite/XpadSize[i]);
      //XYpad[i][j]->SetTopMargin(tp/YpadSize[j]);XYpad[i][j]->SetBottomMargin(bot/YpadSize[j]);
      //XYpad[i][j]->Draw();
      if(y1<0)y1=0;
      XYpad[nx*i+j]=new TPad(Form("pad%d_%d",i,j),"",x1,y1,x2,y2);
      XYpad[nx*i+j]->SetLeftMargin(left/XpadSize[i]);XYpad[nx*i+j]->SetRightMargin(rite/XpadSize[i]);
      XYpad[nx*i+j]->SetTopMargin(tp/YpadSize[j]);XYpad[nx*i+j]->SetBottomMargin(bot/YpadSize[j]);
      XYpad[nx*i+j]->Draw();
      y2=y1;
    }
    x1=x2;
  }
  xyCent[0]=(leftOffset+left+1-right)/2;
  xyCent[1]=(bottomOffset+bottom+1-top)/2;
  return;
}

void Xpads(Float_t bottomOffset, const Int_t n, 
	   Float_t left, Float_t right, Float_t top, Float_t bottom, Float_t gap,
	   TPad** Xpad, Float_t* xyCent, Float_t* XpadSize, Float_t* XpadLeft, Float_t* XpadRite)
//all variables are in parent frame, except XpadLeft and XpadRite are in local frame of pad.
{
  float xsize=(1-left-right-(n-1)*gap)/n;
  float rite=gap/2;
  float x1=0,x2=0;
  for(int i=0;i<n;++i){
    if(i)left=gap/2;
    x2+=left;
    //XpadCent[i]=x2+xsize/2;//this is in parent frame
    x2+=xsize;
    if(i==n-1)rite=right;
    x2+=rite;
    XpadSize[i]=x2-x1;
    XpadLeft[i]=left/XpadSize;//this is in pad frame
    XpadRite[i]=1-rite/XpadSize[i];//this is in pad frame
    //XpadCent[i]=(1-rite/XpadSize[i]+left/XpadSize[i])/2;//this is in pad frame
    if(x2>1)x2=1;
    Xpad[i]=new TPad(Form("pad%d",i),"",x1,bottomOffset,x2,1);
    Xpad[i]->SetLeftMargin(left/XpadSize[i]);Xpad[i]->SetRightMargin(rite/XpadSize[i]);
    Xpad[i]->SetTopMargin(top/(1-bottomOffset));Xpad[i]->SetBottomMargin(bottom/(1-bottomOffset));
    Xpad[i]->Draw();
    x1=x2;
  }
  xyCent[0]=(left+1-right)/2;
  xyCent[1]=(bottomOffset+bottom+1-top)/2;
  return;
}

void Ypads(Float_t leftOffset, const Int_t n, 
	   Float_t left, Float_t right, Float_t top, Float_t bottom, Float_t gap,
	   TPad** Ypad, Float_t* xyCent, Float_t* YpadSize, Float_t* YpadTop, Float_t* YpadBot)
//all variables are in parent frame, except XpadTop and XpadBot are in local frame of pad.
{
  float ysize=(1-top-bottom-(n-1)*gap)/n;
  float bot=gap/2;
  float y1=1,y2=1;
  for(int i=0;i<n;++i){
    if(i)top=gap/2;
    y1-=top;
    //YpadCent[i]=y1-ysize/2;//this is in parent frame
    y1-=ysize;
    if(i==n-1)bot=bottom;
    y1-=bot;
    YpadSize[i]=y2-y1;
    YpadTop[i]=1-top/YpadSize[i];//this is in pad frame
    YpadBot[i]=bot/YpadSize[i];//this is in pad frame
    //YpadCent[i]=(1-top/YpadSize[i]+bot/YpadSize[i])/2;//this is in pad frame
    if(y1<0)y1=0;
    Ypad[i]=new TPad(Form("pad%d",i),"",leftOffset,y1,1,y2);
    Ypad[i]->SetLeftMargin(left/(1-leftOffset));Ypad[i]->SetRightMargin(right/(1-leftOffset));
    Ypad[i]->SetTopMargin(top/YpadSize[i]);Ypad[i]->SetBottomMargin(bot/YpadSize[i]);
    Ypad[i]->Draw();
    y2=y1;
  }
  xyCent[0]=(leftOffset+left+1-right)/2;
  xyCent[1]=(bottom+1-top)/2;
  return;
}
/*
void equalPads(const Int_t nx, const Int_t ny, 
	       Float_t left, Float_t right, Float_t top, Float_t bottom, Float_t xGap, Float_t yGap, 
	       TPad* pad, Float_t* padxSize, Float_t* padySize, Float_t* padxCent, Float_t* padyCent)
{
  float xsize=(1-left-right-(nx-1)*xGap)/nx;
  float ysize=(1-top-bottom-(ny-1)*yGap)/ny;
  float x1=left;
  for(int ix=0;ix<nx;++ix){
    float y1=1-top;
    for(int iy=0;iy<ny;++iy){
      pad[ix][iy]=new TPad(Form("pad_%d_%d",ix,iy),"",0,1-size1,1,1);

      padSize[ix]
    padSizeY[0]=(1+2*top-bottom)/3;


  size2=(1-top-bottom)/3,size3=(1-top+2*bottom)/3;
  TCanvas *cfig=new TCanvas("cfig","",10,10,500,750);
  pad1=new TPad("pad1","",0,1-size1,1,1);
  pad2=new TPad("pad2","",0,size3,1,1-size1);
  pad3=new TPad("pad3","",0,0,1,size3);
  pad1->SetLeftMargin(left);pad1->SetRightMargin(right);
  pad1->SetBottomMargin(0);pad1->SetTopMargin(top/size1);
  pad1->Draw();
  pad2->SetLeftMargin(left);pad2->SetRightMargin(right);
  pad2->SetBottomMargin(0);pad2->SetTopMargin(0);
  pad2->Draw();
  pad3->SetLeftMargin(left);pad3->SetRightMargin(right);
  pad3->SetBottomMargin(bottom/size3);pad3->SetTopMargin(0);
  pad3->Draw();
*/
/////////////////////////////////////////////////////////
void keyLine(Float_t x, Float_t y, const Char_t* tt, 
             Int_t color=1, Int_t style=1, Float_t tSize=0.04, 
	     Int_t lWid=1, Float_t lSize=1)
{
  gPad->Update();
  Float_t x1 = gPad->GetFrame()->GetX1();
  Float_t x2 = gPad->GetFrame()->GetX2();
  Float_t dx = x2-x1;
  x = x1 + dx*x;
  Float_t y1 = gPad->GetFrame()->GetY1();
  Float_t y2 = gPad->GetFrame()->GetY2();
  Float_t dy = y2-y1;
  y = y1 + dy*y;
  Float_t xx=x+0.05*lSize*dx;
  if(gPad->GetLogx()){x=pow(10,x);xx=pow(10,xx);}
  if(gPad->GetLogy())y=pow(10,y);
  TLine *l = new TLine(x,y,xx,y);
  l->SetLineColor(color);
  l->SetLineStyle(style);
  l->SetLineWidth(lWid);
  l->Draw();
  TLatex *t = new TLatex();
  //t->SetNDC(1);
  t->SetTextAlign(12);
  t->SetTextSize(tSize);
  t->SetTextColor(color);
  //printf("%f %f\n",x+0.07*lSize*dx,y);
  //t->DrawLatex(x+0.07*lSize*dx,y,tt);
  t->DrawLatex(xx+0.02*dx,y,tt);
  //t->DrawLatex(x+0.05*lSize+0.02,y,tt);
}

void keySymb(Float_t x, Float_t y, const Char_t* tt, 
               Int_t color=1, Int_t marker=1, Float_t tSize=0.04, Float_t mSize=1)
{
  gPad->Update();
  Float_t x1 = gPad->GetFrame()->GetX1();
  Float_t x2 = gPad->GetFrame()->GetX2();
  Float_t dx = x2-x1;
  x = x1 + dx*x;
  Float_t y1 = gPad->GetFrame()->GetY1();
  Float_t y2 = gPad->GetFrame()->GetY2();
  Float_t dy = y2-y1;
  y = y1 + dy*y;
  TMarker *l = new TMarker(x+0.025*dx,y,marker);
  l->SetMarkerColor(color);
  l->SetMarkerSize(mSize);
  l->Draw();
  //TText *t = new TText(x+0.06*dx,y,tt);
  TLatex *t = new TLatex();
  t->SetTextAlign(12);
  t->SetTextSize(tSize);
  t->SetTextColor(color);
  //t->Draw();
  t->DrawLatex(x+0.06*dx,y,tt);
}

TMarker *keySymbol(Float_t x, Float_t y, const Char_t* tt, 
		   Int_t color=1, Int_t marker=1, Float_t tSize=0.04, Float_t mSize=1)
{
  TMarker *m = new TMarker(x,y,marker);
  m->SetNDC(1);
  m->SetMarkerColor(color);
  m->SetMarkerSize(mSize);
  m->Draw();
  TLatex *t = new TLatex();
  t->SetNDC(1);
  t->SetTextAlign(12);
  t->SetTextColor(color);
  t->SetTextSize(tSize);
  t->DrawLatex(x+0.005*mSize+0.02,y,tt);
  return(m);
}

void keyText(Float_t x, Float_t y, const Char_t* tt, 
	  Int_t color=1, Float_t tSize=0.04, Int_t tAlign=12)
{
  gPad->Update();
  Float_t x1 = gPad->GetFrame()->GetX1();
  Float_t x2 = gPad->GetFrame()->GetX2();
  Float_t dx = x2-x1;
  x = x1 + dx*x;
  Float_t y1 = gPad->GetFrame()->GetY1();
  Float_t y2 = gPad->GetFrame()->GetY2();
  Float_t dy = y2-y1;
  y = y1 + dy*y;
  TText *txt = new TText();
  txt->SetTextAlign(tAlign);
  txt->SetTextSize(tSize);
  txt->SetTextColor(color);
  txt->DrawText(x,y,tt);
}

void wLatex(Float_t x, Float_t y, const Char_t* tt, 
	    Int_t color=1, Int_t tAlign=12, Float_t tSize, Float_t tAngle)
{
  TLatex *t = new TLatex();
  t->SetNDC(1);
  t->SetTextColor(color);
  t->SetTextAlign(tAlign);
  t->SetTextSize(tSize);
  t->SetTextAngle(tAngle);
  t->DrawLatex(x,y,tt);
}

/////////////////////////////////////////////////////////

int utilityPrint=1; //print absolute error
//int utilityPrint=-1; //print relative error
int utilityPrint=0;

void graphSystCap(int n, Float_t *x, Float_t *y, Float_t *ex, Float_t *ey, 
		  Float_t *syPlus=0, Float_t *syMinus=0, 
		  Float_t dxaxis=1, Float_t dyaxis=1, 
		  Int_t symb=20, Float_t sSize=1, Int_t color=1, Int_t lwidth=1)
{
  if (utilityPrint) printf("\n__________graphSystCap__________symb=%d sSize=%f color=%d__________\n",symb,sSize,color);
  if (utilityPrint) {
    for (int i=0; i<n; ++i) {
      printf("x[%d]=%f; y[%d]=%f; ",i,x[i],i,y[i]);
      if (utilityPrint<0) {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]/x[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]/y[i]);
	if (syPlus) printf("sy1[%d]=%f; ",i,syPlus[i]/y[i]);
	if (syMinus) printf("sy2[%d]=%f; ",i,syMinus[i]/y[i]);
      }
      else {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]);
	if (syPlus) printf("sy1[%d]=%f; ",i,syPlus[i]);
	if (syMinus) printf("sy2[%d]=%f; ",i,syMinus[i]);
      }
      printf("\n");
    }
  }
  if(symb){
    TGraphErrors *graf = new TGraphErrors(n,x,y,ex,ey);
    graf->SetMarkerStyle(symb);
    graf->SetMarkerSize(sSize);
    graf->SetMarkerColor(color);
    graf->SetLineColor(color);
    graf->SetLineWidth(lwidth);
    graf->Draw("p");
  }
  if (!syPlus) return;
  float scale=0.01*sSize;
  if (scale<0.01) scale=0.01;
  //double xmin,xmax,ymin,ymax;
  //gPad->GetRangeAxis(xmin,ymin,xmax,ymax);
  //float xmin=gPad->GetUxmin();
  //float xmax=gPad->GetUxmax();
  //float ymin=gPad->GetUymin();
  //float ymax=gPad->GetUymax();
  //printf("%f %f %f %f\n",xmin,ymin,xmax,ymax);
  //float dx=(xmax-xmin)*scale;
  //float dy=(ymax-ymin)*scale;
  float dx=dxaxis*scale;
  float dy=dyaxis*scale;
  TLine *line=new TLine();
  line->SetLineColor(color);
  line->SetLineWidth(lwidth);
  for (int i=0; i<n; ++i) {
    float x1=x[i]-dx;
    float x2=x[i]+dx;
    float y1=y[i]+syPlus[i];
    if (syMinus) float y2=y[i]+syMinus[i];
    else         float y2=y[i]-syPlus[i];
    line->DrawLine(x1,y1,x2,y1);
    line->DrawLine(x1,y2,x2,y2);
    line->DrawLine(x1,y1,x1,y1-dy);
    line->DrawLine(x2,y1,x2,y1-dy);
    line->DrawLine(x1,y2,x1,y2+dy);
    line->DrawLine(x2,y2,x2,y2+dy);
  }
}

void tfBand(TF1 *f1,TF1 *f2,float x1=0,float x2=1,int bcol=5,int lwid=1,int lcol=1)
{
  const int n=1000;
  float xband[2*n], yband[2*n];
  float dx=(x2-x1)/n;
  float x=x1+0.5*dx;
  for (int i=0; i<n; ++i) {
    xband[i]=x;
    xband[2*n-1-i]=x;
    yband[i]=f1->Eval(x);
    yband[2*n-1-i]=f2->Eval(x);
    x+=dx;
  }
  TGraph *grafband = new TGraph(2*n,xband,yband);
  grafband->SetFillColor(bcol);
  grafband->Draw("f");
  if(!lwid)return;
  f1->SetLineColor(lcol);f1->SetLineWidth(lwid);
  f2->SetLineColor(lcol);f2->SetLineWidth(lwid);
  f1->Draw("same");
  f2->Draw("same");
}

void graphSystBand(int n, Float_t *x, Float_t *y, Float_t *ex, Float_t *ey, 
		   Float_t *syPlus=0, Float_t *syMinus=0, Int_t colorBand=5,
		   Int_t symb=20, Float_t sSize=1, Int_t color=1, Int_t lwidth=1)
{
  if (utilityPrint) printf("\n__________graphSysErr2band__________symb=%d sSize=%f color=%d colorBand__________\n",symb,sSize,color,colorBand);
  if (utilityPrint) {
    for (int i=0; i<n; ++i) {
      printf("x[%d]=%f; y[%d]=%f; ",i,x[i],i,y[i]);
      if (utilityPrint<0) {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]/x[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]/y[i]);
	if (syPlus ) printf("sy1[%d]=%f; ",i,syPlus [i]/y[i]);
	if (syMinus) printf("sy2[%d]=%f; ",i,syMinus[i]/y[i]);
      }
      else {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]);
	if (syPlus ) printf("sy1[%d]=%f; ",i,syPlus [i]);
	if (syMinus) printf("sy2[%d]=%f; ",i,syMinus[i]);
      }
      printf("\n");
    }
  }
  if (syPlus) {
    float xband[1000], yband[1000];
    for (int i=0; i<n; ++i) {
      xband[i]=x[i];
      xband[2*n-1-i]=x[i];
      float y1=y[i]+syPlus[i];
      if (syMinus) float y2=y[i]+syMinus[i];
      else         float y2=y[i]-syPlus[i];
      if (y1<y2) {float yy=y2; y2=y1; y1=yy;} 
      yband[i]=y1;
      yband[2*n-1-i]=y2;
      //printf("%d %d\n",i,2*n-1-i);
    }
    //for(int iii=0; iii<2*n;++iii)printf("%f %f\n",xband[iii],yband[iii]);
    TGraph *grafband = new TGraph(2*n,xband,yband);
    grafband->SetFillColor(colorBand);
    grafband->Draw("f");
  }
  if(!symb)return;
  TGraphErrors *graf = new TGraphErrors(n,x,y,ex,ey);
  graf->SetMarkerStyle(symb);
  graf->SetMarkerSize(sSize);
  graf->SetMarkerColor(color);
  graf->SetLineColor(color);
  graf->SetLineWidth(lwidth);
  graf->Draw("p");
}

void graphSystLine(int n, Float_t *x, Float_t *y, Float_t *ex, Float_t *ey, 
		   Float_t *syPlus=0, Float_t *syMinus=0, Int_t colorLine=5,
		   Int_t symb=20, Float_t sSize=1, Int_t color=1, Int_t lwidth=1)
{
  if (utilityPrint) printf("\n__________graphSystLine__________symb=%d sSize=%f color=%d colorLine=%d__________\n",symb,sSize,color,colorLine);
  if (utilityPrint) {
    for (int i=0; i<n; ++i) {
      printf("x[%d]=%f; y[%d]=%f; ",i,x[i],i,y[i]);
      if (utilityPrint<0) {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]/x[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]/y[i]);
	if (syPlus ) printf("sy1[%d]=%f; ",i,syPlus [i]/y[i]);
	if (syMinus) printf("sy2[%d]=%f; ",i,syMinus[i]/y[i]);
      }
      else {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]);
	if (syPlus ) printf("sy1[%d]=%f; ",i,syPlus [i]);
	if (syMinus) printf("sy2[%d]=%f; ",i,syMinus[i]);
      }
      printf("\n");
    }
  }
  if (syPlus) {
    float xline[1000], yline[1000];
    for (int i=0; i<n; ++i) {
      xline[i]=x[i];
      xline[2*n-1-i]=x[i];
      float y1=y[i]+syPlus[i];
      if (syMinus) float y2=y[i]+syMinus[i];
      else         float y2=y[i]-syPlus[i];
      if (y1<y2) {float yy=y2; y2=y1; y1=yy;} 
      yline[i]=y1;
      yline[2*n-1-i]=y2;
      //printf("%d %d\n",i,2*n-1-i);
    }
    xline[2*n]=xline[0];
    yline[2*n]=yline[0];
    //for(int iii=0; iii<2*n+1;++iii)printf("%f %f\n",xline[iii],yline[iii]);
    TGraph *grafline = new TGraph(2*n+1,xline,yline);
    grafline->SetLineColor(colorLine);
    grafline->Draw("l");
  }
  if(!symb)return;
  TGraphErrors *graf = new TGraphErrors(n,x,y,ex,ey);
  graf->SetMarkerStyle(symb);
  graf->SetMarkerSize(sSize);
  graf->SetMarkerColor(color);
  graf->SetLineColor(color);
  graf->SetLineWidth(lwidth);
  graf->Draw("p");
}

void graphSystBox(int n, Float_t *x, Float_t *y, Float_t *ex, Float_t *ey,
		  Float_t *syPlus=0, Float_t *syMinus=0, Int_t colorBox=5,
		  Float_t dxaxis=1, Int_t symb=20, Float_t sSize=1, Int_t color=1, Int_t lwidth=1)
{
  if (utilityPrint) printf("\n__________graphSystBox__________symb=%d sSize=%f color=%d colorBox=%d__________\n",symb,sSize,color,colorBox);
  if (utilityPrint) {
    for (int i=0; i<n; ++i) {
      printf("x[%d]=%f; y[%d]=%f; ",i,x[i],i,y[i]);
      if (utilityPrint<0) {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]/x[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]/y[i]);
	if (syPlus ) printf("sy1[%d]=%f; ",i,syPlus [i]/y[i]);
	if (syMinus) printf("sy2[%d]=%f; ",i,syMinus[i]/y[i]);
      }
      else {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]);
	if (syPlus ) printf("sy1[%d]=%f; ",i,syPlus [i]);
	if (syMinus) printf("sy2[%d]=%f; ",i,syMinus[i]);
      }
      printf("\n");
    }
  }
  if (syPlus) {
    float scale=0.015*sSize;
    if (scale<0.01) scale=0.01;
    float dx=dxaxis*scale;
    float xbox[5], ybox[5];
    for (int i=0; i<n; ++i) {
      float y1=y[i]+syPlus[i];
      if (syMinus) float y2=y[i]+syMinus[i];
      else         float y2=y[i]-syPlus[i];
      xbox[0]=x[i]-dx;ybox[0]=y1;
      xbox[1]=x[i]+dx;ybox[1]=y1;
      xbox[2]=x[i]+dx;ybox[2]=y2;
      xbox[3]=x[i]-dx;ybox[3]=y2;
      xbox[4]=x[i]-dx;ybox[4]=y1;
      TGraph *grafbox=new TGraph(5,xbox,ybox);
      if(colorBox>=0){
	grafbox->SetFillColor(colorBox);
	grafbox->Draw("f");
      }
      else{
	grafbox->SetLineWidth(-colorBox);
	grafbox->SetLineColor(color);
	grafbox->Draw("l");
      }
    }
  }
  if(!symb)return;
  TGraphErrors *graf = new TGraphErrors(n,x,y,ex,ey);
  graf->SetMarkerStyle(symb);
  graf->SetMarkerSize(sSize);
  graf->SetMarkerColor(color);
  graf->SetLineColor(color);
  graf->SetLineWidth(lwidth);
  graf->Draw("p");
}

void graphSystBox(int n, Float_t *x, Float_t *y, Float_t *ex, Float_t *ey,
		  Float_t *syPlus=0, Float_t *syMinus=0, Int_t colorBox=5,
		  Float_t* dx, Int_t symb=20, Float_t sSize=1, Int_t color=1, Int_t lwidth=1)
{
  if (utilityPrint) printf("\n__________graphSystBox__________symb=%d sSize=%f color=%d colorBox=%d__________\n",symb,sSize,color,colorBox);
  if (utilityPrint) {
    for (int i=0; i<n; ++i) {
      printf("x[%d]=%f; y[%d]=%f; ",i,x[i],i,y[i]);
      if (utilityPrint<0) {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]/x[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]/y[i]);
	if (syPlus ) printf("sy1[%d]=%f; ",i,syPlus [i]/y[i]);
	if (syMinus) printf("sy2[%d]=%f; ",i,syMinus[i]/y[i]);
      }
      else {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]);
	if (syPlus ) printf("sy1[%d]=%f; ",i,syPlus [i]);
	if (syMinus) printf("sy2[%d]=%f; ",i,syMinus[i]);
      }
      printf("\n");
    }
  }
  if (syPlus) {
    float xbox[5], ybox[5];
    for (int i=0; i<n; ++i) {
      float y1=y[i]+syPlus[i];
      if (syMinus) float y2=y[i]+syMinus[i];
      else         float y2=y[i]-syPlus[i];
      xbox[0]=x[i]-dx[i];ybox[0]=y1;
      xbox[1]=x[i]+dx[i];ybox[1]=y1;
      xbox[2]=x[i]+dx[i];ybox[2]=y2;
      xbox[3]=x[i]-dx[i];ybox[3]=y2;
      xbox[4]=x[i]-dx[i];ybox[4]=y1;
      TGraph *grafbox=new TGraph(5,xbox,ybox);
      if(colorBox>=0){
	grafbox->SetFillColor(colorBox);
	grafbox->Draw("f");
      }
      else{
	grafbox->SetLineWidth(-colorBox);
	grafbox->SetLineColor(color);
	grafbox->Draw("l");
      }
    }
  }
  if(!symb)return;
  TGraphErrors *graf = new TGraphErrors(n,x,y,ex,ey);
  graf->SetMarkerStyle(symb);
  graf->SetMarkerSize(sSize);
  graf->SetMarkerColor(color);
  graf->SetLineColor(color);
  graf->SetLineWidth(lwidth);
  graf->Draw("p");
}

void graphSystArrow(int n, Float_t *x, Float_t *y, Float_t *ex, Float_t *ey, 
		    Float_t *sySingleSided=0,
		    Int_t symb=20, Float_t sSize=1, Int_t color=1, Int_t lwidth=1)
{
  if (utilityPrint) printf("\n__________graphSystArrow__________symb=%d sSize=%f color=%d__________\n",symb,sSize,color);
  if (utilityPrint) {
    for (int i=0; i<n; ++i) {
      printf("x[%d]=%f; y[%d]=%f; ",i,x[i],i,y[i]);
      if (utilityPrint<0) {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]/x[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]/y[i]);
	if (sySingleSided) printf("sy[%d]=%f; ",i,sySingleSided[i]/y[i]);
      }
      else {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]);
	if (sySingleSided) printf("sy[%d]=%f; ",i,sySingleSided[i]);
      }
      printf("\n");
    }
  }
  if (sySingleSided) {
    TArrow *arrow=new TArrow();
    arrow->SetLineColor(color);arrow->SetFillColor(color);arrow->SetLineWidth(3);
    for (int i=0; i<n; ++i) {
      arrow->DrawArrow(x[i],y[i],x[i],y[i]+sySingleSided[i],0.01);
    }
  }
  if(!symb)return;
  TGraphErrors *graf = new TGraphErrors(n,x,y,ex,ey);
  graf->SetMarkerStyle(symb);
  graf->SetMarkerSize(sSize);
  graf->SetMarkerColor(color);
  graf->SetLineColor(color);
  graf->SetLineWidth(lwidth);
  graf->Draw("p");
}


/////////////////////////////////////////////////////////
void graphLine(int n, Float_t *x, Float_t *y, Int_t style=1, Int_t color=1, Int_t width=1, Int_t marker=1, Int_t size=1)
{
  graphLine(n,x,y,0,0,style,color,width,marker,size);
}

void graphLine(int n, Float_t *x, Float_t *y, Float_t *ex, Float_t *ey, Int_t style=1, Int_t color=1, Int_t width=1, Int_t marker=1, Int_t size=1)
{
  TGraphErrors *graf = new TGraphErrors(n,x,y,ex,ey);
  graf->SetLineColor(color);
  graf->SetLineWidth(width);
  graf->SetMarkerStyle(marker);
  graf->SetMarkerSize(size);
  graf->SetMarkerColor(color);
  graf->SetLineStyle(1);
  graf->Draw("p");
  if(!style)return;
  TGraph *graf2 = new TGraph(n,x,y);
  graf2->SetLineColor(color);
  graf2->SetLineWidth(width);
  graf2->SetMarkerStyle(marker);
  graf2->SetMarkerSize(size);
  graf2->SetMarkerColor(color);
  graf2->SetLineStyle(style);
  graf2->Draw("l");
}

void graphAsymLine(int n, Float_t *x, Float_t *y, Float_t *explus, Float_t *exminus, Float_t *eyplus, Float_t *eyminus, Int_t style=1, Int_t color=1, Int_t width=1, Int_t marker=1, Int_t size=1)
{
  TGraphAsymmErrors *graf = new TGraphAsymmErrors(n,x,y,exminus,explus,eyminus,eyplus);
  graf->SetLineColor(color);
  graf->SetLineWidth(width);
  graf->SetMarkerStyle(marker);
  graf->SetMarkerSize(size);
  graf->SetMarkerColor(color);
  graf->SetLineStyle(1);
  graf->Draw("p");
  if(!style)return;
  TGraph *graf2 = new TGraph(n,x,y);
  graf2->SetLineColor(color);
  graf2->SetLineWidth(width);
  graf2->SetMarkerStyle(marker);
  graf2->SetMarkerSize(size);
  graf2->SetMarkerColor(color);
  graf2->SetLineStyle(style);
  graf2->Draw("l");
}

void graphBand(int n, Float_t *x, Float_t *y1, Float_t *y2, 
	       Int_t colorBand=5, Int_t lcolor=0, Int_t lwidth=1, Int_t lstyle=1)
{
  float xband[1000], yband[1000];
  for (int i=0; i<n; ++i) {
    xband[i]=x[i];
    yband[i]=y1[i];
    xband[2*n-1-i]=x[i];
    yband[2*n-1-i]=y2[i];
  }
  TGraph *grafband = new TGraph(2*n,xband,yband);
  grafband->SetFillColor(colorBand);
  grafband->Draw("f");
  if(!lcolor)return;
  graphLine(n,x,y1,lstyle,lcolor,lwidth);
  graphLine(n,x,y2,lstyle,lcolor,lwidth);
}


/////////////////////////////////////////////////////////
void graphSysErr(int n, Float_t *x, Float_t *y, Float_t *ex, Float_t *ey,
		 Float_t *syPlus=0, Float_t *syMinus=0, 
		 Int_t symb=20, Float_t sSize=1, Int_t color=1)
{
  if (utilityPrint) printf("\n__________graphSysErr__________symb=%d sSize=%f color=%d__________\n",symb,sSize,color);
  if (utilityPrint) {
    for (int i=0; i<n; ++i) {
      printf("x[%d]=%f; y[%d]=%f; ",i,x[i],i,y[i]);
      if (utilityPrint<0) {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]/x[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]/y[i]);
	if (syPlus) printf("sy1[%d]=%f; ",i,syPlus[i]/y[i]);
	if (syMinus) printf("sy2[%d]=%f; ",i,syMinus[i]/y[i]);
      }
      else {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]);
	if (syPlus) printf("sy1[%d]=%f; ",i,syPlus[i]);
	if (syMinus) printf("sy2[%d]=%f; ",i,syMinus[i]);
      }
      printf("\n");
    }
  }
  TGraphErrors *graf = new TGraphErrors(n,x,y,ex,ey);
  graf->SetMarkerStyle(symb);
  graf->SetMarkerSize(sSize);
  graf->SetMarkerColor(color);
  graf->SetLineColor(color);
  graf->Draw("p");
  if (!syPlus) return;
  float scale=0.01*sSize;
  if (scale<0.01) scale=0.01;
  float xmin=gPad->GetUxmin();
  float xmax=gPad->GetUxmax();
  float dx=(xmax-xmin)*scale;
  float ymin=gPad->GetUymin();
  float ymax=gPad->GetUymax();
  float dy=(ymax-ymin)*scale;
  TLine *line=new TLine();
  line->SetLineColor(color);
  for (int i=0; i<n; ++i) {
    float x1=x[i]-dx;
    float x2=x[i]+dx;
    float y1=y[i]+syPlus[i];
    if (syMinus) float y2=y[i]+syMinus[i];
    else         float y2=y[i]-syPlus[i];
    line->DrawLine(x1,y1,x2,y1);
    line->DrawLine(x1,y2,x2,y2);
    line->DrawLine(x1,y1,x1,y1-dy);
    line->DrawLine(x2,y1,x2,y1-dy);
    line->DrawLine(x1,y2,x1,y2+dy);
    line->DrawLine(x2,y2,x2,y2+dy);
  }
}

void graphSysErr2(int n, Float_t *x, Float_t *y, Float_t *ex, Float_t *ey, 
		  Float_t *syPlus=0, Float_t *syMinus=0, 
		  Float_t dxaxis=1, Float_t dyaxis=1, 
		  Int_t symb=20, Float_t sSize=1, Int_t color=1)
{//same as graphSysErr except dxaxis,dyaxis NOT obtained automatically.
  if (utilityPrint) printf("\n__________graphSysErr2__________symb=%d sSize=%f color=%d__________\n",symb,sSize,color);
  if (utilityPrint) {
    for (int i=0; i<n; ++i) {
      printf("x[%d]=%f; y[%d]=%f; ",i,x[i],i,y[i]);
      if (utilityPrint<0) {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]/x[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]/y[i]);
	if (syPlus) printf("sy1[%d]=%f; ",i,syPlus[i]/y[i]);
	if (syMinus) printf("sy2[%d]=%f; ",i,syMinus[i]/y[i]);
      }
      else {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]);
	if (syPlus) printf("sy1[%d]=%f; ",i,syPlus[i]);
	if (syMinus) printf("sy2[%d]=%f; ",i,syMinus[i]);
      }
      printf("\n");
    }
  }
  TGraphErrors *graf = new TGraphErrors(n,x,y,ex,ey);
  graf->SetMarkerStyle(symb);
  graf->SetMarkerSize(sSize);
  graf->SetMarkerColor(color);
  graf->SetLineColor(color);
  graf->Draw("p");
  if (!syPlus) return;
  float scale=0.01*sSize;
  if (scale<0.01) scale=0.01;
  /*
  float xmin=gPad->GetUxmin();
  float xmax=gPad->GetUxmax();
  float dx=(xmax-xmin)*scale;
  float ymin=gPad->GetUymin();
  float ymax=gPad->GetUymax();
  float dy=(ymax-ymin)*scale;
  */
  float dx=dxaxis*scale;
  float dy=dyaxis*scale;
  TLine *line=new TLine();
  line->SetLineColor(color);
  for (int i=0; i<n; ++i) {
    float x1=x[i]-dx;
    float x2=x[i]+dx;
    float y1=y[i]+syPlus[i];
    if (syMinus) float y2=y[i]+syMinus[i];
    else         float y2=y[i]-syPlus[i];
    line->DrawLine(x1,y1,x2,y1);
    line->DrawLine(x1,y2,x2,y2);
    line->DrawLine(x1,y1,x1,y1-dy);
    line->DrawLine(x2,y1,x2,y1-dy);
    line->DrawLine(x1,y2,x1,y2+dy);
    line->DrawLine(x2,y2,x2,y2+dy);
  }
}

void graphSysErr3(int n, Float_t *x, Float_t *y, Float_t *ex, Float_t *ey, 
		  Float_t *syUpper=0, Float_t *syLower=0, 
		  Float_t dxaxis=1, Float_t dyaxis=1, 
		  Int_t symb=20, Float_t sSize=1, Int_t color=1)
{//same as graphSysErr2 except syUpper,syLower are bounds, NOT +/-.
  if (utilityPrint) printf("\n__________graphSysErr3__________symb=%d sSize=%f color=%d__________\n",symb,sSize,color);
  if (utilityPrint) {
    for (int i=0; i<n; ++i) {
      printf("x[%d]=%f; y[%d]=%f; ",i,x[i],i,y[i]);
      if (utilityPrint<0) {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]/x[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]/y[i]);
	if (syUpper) printf("sy1[%d]=%f; ",i,syUpper[i]/y[i]-1);
	if (syLower) printf("sy2[%d]=%f; ",i,syLower[i]/y[i]-1);
      }
      else {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]);
	if (syUpper) printf("sy1[%d]=%f; ",i,syUpper[i]-y[i]);
	if (syLower) printf("sy2[%d]=%f; ",i,syLower[i]-y[i]);
      }
      printf("\n");
    }
  }
  TGraphErrors *graf = new TGraphErrors(n,x,y,ex,ey);
  graf->SetMarkerStyle(symb);
  graf->SetMarkerSize(sSize);
  graf->SetMarkerColor(color);
  graf->SetLineColor(color);
  graf->Draw("p");
  if (!syUpper) return;
  float scale=0.01*sSize;
  if (scale<0.01) scale=0.01;
  /*
  float xmin=gPad->GetUxmin();
  float xmax=gPad->GetUxmax();
  float dx=(xmax-xmin)*scale;
  float ymin=gPad->GetUymin();
  float ymax=gPad->GetUymax();
  float dy=(ymax-ymin)*scale;
  */
  float dx=dxaxis*scale;
  float dy=dyaxis*scale;
  TLine *line=new TLine();
  line->SetLineColor(color);
  for (int i=0; i<n; ++i) {
    float x1=x[i]-dx;
    float x2=x[i]+dx;
    float y1=syUpper[i];
    if (syLower) float y2=syLower[i];
    else         float y2=y[i]-(syUpper[i]-y[i]);
    if (y1<y2) {float yy=y2; y2=y1; y1=yy;} 
    line->DrawLine(x1,y1,x2,y1);
    line->DrawLine(x1,y1,x1,y1-dy);
    line->DrawLine(x2,y1,x2,y1-dy);
    line->DrawLine(x1,y2,x2,y2);
    line->DrawLine(x1,y2,x1,y2+dy);
    line->DrawLine(x2,y2,x2,y2+dy);
  }
}

void graphSysErr3a(int n, Float_t *x, Float_t *y, Float_t *ex, Float_t *ey, 
		   Float_t *syUpper=0, Float_t *syLower=0, 
		   Float_t dxaxis=1, Float_t dyaxis=1, 
		   Int_t symb=20, Float_t sSize=1, Int_t *color)
{
  for (int i=0; i<n; ++i) {
    graphSysErr3(n, x, y, ex, ey, syUpper, syLower, dxaxis, dyaxis, 
		 symb, sSize, color[i]);
  }
}

void graphSysErr2logy(int n, Float_t *x, Float_t *y, Float_t *ex, Float_t *ey, 
		  Float_t *syPlus=0, Float_t *syMinus=0, 
		  Float_t dxaxis=1, Float_t dyaxis=1, 
		  Int_t symb=20, Float_t sSize=1, Int_t color=1)
{
  if (!syPlus) return;
  float syUpTmp[1000], syLoTmp[1000];
  for (int i=0; i<n; ++i) {
    if (y[i]+syPlus[i]>0) syUpTmp[i]=log10(y[i]+syPlus[i]);
    else syUpTmp[i]=-100;
    if (syMinus) {
      if (y[i]+syMinus[i]>0) syLoTmp[i]=log10(y[i]+syMinus[i]);
      else syLoTmp[i]=-100;
    }
    else {
      if (y[i]-syPlus[i]>0) syLoTmp[i]=log10(y[i]-syPlus[i]);
      else syLoTmp[i]=-100;
    }
    //printf("%d %f %f\n",i,syUpTmp[i],syLoTmp[i]);
  }
  graphSysErr3(n,x,y,ex,ey,syUpTmp,syLoTmp,dxaxis,3*log10(dyaxis),symb,sSize,color);
}

void graphSysErr3logy(int n, Float_t *x, Float_t *y, Float_t *ex, Float_t *ey, 
		  Float_t *syUpper=0, Float_t *syLower=0, 
		  Float_t dxaxis=1, Float_t dyaxis=1, 
		  Int_t symb=20, Float_t sSize=1, Int_t color=1)
{
  if (!syUpper) return;
  float syUpTmp[1000], syLoTmp[1000];
  for (int i=0; i<n; ++i) {
    if (syUpper[i]>0) syUpTmp[i]=log10(syUpper[i]);
    else syUpTmp[i]=-100;
    if (syLower) {
      if (syLower[i]>0) syLoTmp[i]=log10(syLower[i]);
      else syLoTmp[i]=-100;
    }
    else {
      if (2*y[i]-syUpper[i]>0) syLoTmp[i]=log10(2*y[i]-syUpper[i]);
      else syLoTmp[i]=-100;
    }
    //printf("%d %f %f\n",i,syUpTmp[i],syLoTmp[i]);
  }
  graphSysErr3(n,x,y,ex,ey,syUpTmp,syLoTmp,dxaxis,3*log10(dyaxis),symb,sSize,color);
}

void graphSysErr2band(int n, Float_t *x, Float_t *y, Float_t *ex, Float_t *ey, 
		      Float_t *syPlus=0, Float_t *syMinus=0, 
		      Float_t dxaxis=1, Int_t symb=20, Float_t sSize=1, 
		      Int_t color=1, Int_t colorband=5)
{
  if (utilityPrint) printf("\n__________graphSysErr2band__________symb=%d sSize=%f color=%d colorband=%d__________\n",symb,sSize,color,colorband);
  if (utilityPrint) {
    for (int i=0; i<n; ++i) {
      printf("x[%d]=%f; y[%d]=%f; ",i,x[i],i,y[i]);
      if (utilityPrint<0) {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]/x[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]/y[i]);
	if (syPlus ) printf("sy1[%d]=%f; ",i,syPlus [i]/y[i]);
	if (syMinus) printf("sy2[%d]=%f; ",i,syMinus[i]/y[i]);
      }
      else {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]);
	if (syPlus ) printf("sy1[%d]=%f; ",i,syPlus [i]);
	if (syMinus) printf("sy2[%d]=%f; ",i,syMinus[i]);
      }
      printf("\n");
    }
  }
  if (syPlus) {
    float scale=0.01*sSize;
    if (scale<0.01) scale=0.01;
    float dx=dxaxis*scale;
    float xband[1000], yband[1000];
    for (int i=0; i<n; ++i) {
      float xtmp=x[i];
      //if (i==0) xtmp-=dx;
      //if (i==n-1) xtmp+=dx;
      xband[i]=xtmp;
      xband[2*n-1-i]=xtmp;
      float y1=y[i]+syPlus[i];
      if (syMinus) float y2=y[i]+syMinus[i];
      else         float y2=y[i]-syPlus[i];
      if (y1<y2) {float yy=y2; y2=y1; y1=yy;} 
      yband[i]=y1;
      yband[2*n-1-i]=y2;
      //printf("%d %d\n",i,2*n-1-i);
    }
    //for(int iii=0; iii<2*n;++iii)printf("%f %f\n",xband[iii],yband[iii]);
    TGraph *grafband = new TGraph(2*n,xband,yband);
    grafband->SetFillColor(colorband);
    grafband->Draw("f");
  }
  if (ey) {
    TGraphErrors *graf = new TGraphErrors(n,x,y,ex,ey);
    graf->SetMarkerStyle(symb);
    graf->SetMarkerSize(sSize);
    graf->SetMarkerColor(color);
    graf->SetLineColor(color);
    graf->Draw("p");
  }
}

void graphSysErr3band(int n, Float_t *x, Float_t *y, Float_t *ex, Float_t *ey, 
		      Float_t *syUpper=0, Float_t *syLower=0, 
		      Float_t dxaxis=1, Int_t symb=20, Float_t sSize=1, 
		      Int_t color=1, Int_t colorband=5)
{
  if (utilityPrint) printf("\n__________graphSysErr3band__________symb=%d sSize=%f color=%d colorband=%d__________\n",symb,sSize,color,colorband);
  if (utilityPrint) {
    for (int i=0; i<n; ++i) {
      printf("x[%d]=%f; y[%d]=%f; ",i,x[i],i,y[i]);
      if (utilityPrint<0) {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]/x[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]/y[i]);
	if (syUpper) printf("sy1[%d]=%f; ",i,syUpper[i]/y[i]-1);
	if (syLower) printf("sy2[%d]=%f; ",i,syLower[i]/y[i]-1);
      }
      else {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]);
	if (syUpper) printf("sy1[%d]=%f; ",i,syUpper[i]-y[i]);
	if (syLower) printf("sy2[%d]=%f; ",i,syLower[i]-y[i]);
      }
      printf("\n");
    }
  }
  if (syUpper) {
    float scale=0.01*sSize;
    if (scale<0.01) scale=0.01;
    float dx=dxaxis*scale;
    float xband[1000], yband[1000];
    for (int i=0; i<n; ++i) {
      float xtmp=x[i];
      //if (i==0) xtmp-=dx;
      //if (i==n-1) xtmp+=dx;
      xband[i]=xtmp;
      xband[2*n-1-i]=xtmp;
      float y1=syUpper[i];
      if (syLower) float y2=syLower[i];
      else         float y2=y[i]-(syUpper[i]-y[i]);
      if (y1<y2) {float yy=y2; y2=y1; y1=yy;} 
      yband[i]=y1;
      yband[2*n-1-i]=y2;
      //printf("%d %d\n",i,2*n-1-i);
    }
    //for(int iii=0; iii<2*n;++iii)printf("%f %f\n",xband[iii],yband[iii]);
    TGraph *grafband = new TGraph(2*n,xband,yband);
    grafband->SetFillColor(colorband);
    grafband->Draw("f");
  }
  if (ey) {
    TGraphErrors *graf = new TGraphErrors(n,x,y,ex,ey);
    graf->SetMarkerStyle(symb);
    graf->SetMarkerSize(sSize);
    graf->SetMarkerColor(color);
    graf->SetLineColor(color);
    graf->Draw("p");
  }
}

void graphSysErr3(int n, Float_t *x, Float_t *y, 
		  Float_t *ex, Float_t *e1y, Float_t *e2y, 
		  Float_t *syUpper=0, Float_t *syLower=0, 
		  Float_t dxaxis=1, Float_t dyaxis=1, 
		  Int_t symb=20, Float_t sSize=1, Int_t color=1)
{ //same as graphSysErr2 except syUpper,syLower are bounds, NOT +/-.
  //same as the other graphSysErr3 except stat. errors are asymmetric.
  if (utilityPrint) printf("\n__________graphSysErr3__________symb=%d sSize=%f color=%d__________\n",symb,sSize,color);
  if (utilityPrint) {
    for (int i=0; i<n; ++i) {
      printf("x[%d]=%f; y[%d]=%f; ",i,x[i],i,y[i]);
      if (utilityPrint<0) {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]/x[i]);
	if (e1y) printf("ey1[%d]=%f; ",i,e1y[i]/y[i]);
	if (e2y) printf("ey2[%d]=%f; ",i,e2y[i]/y[i]);
	if (syUpper) printf("sy1[%d]=%f; ",i,syUpper[i]/y[i]-1);
	if (syLower) printf("sy2[%d]=%f; ",i,syLower[i]/y[i]-1);
      }
      else {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]);
	if (e1y) printf("ey1[%d]=%f; ",i,e1y[i]);
	if (e2y) printf("ey2[%d]=%f; ",i,e2y[i]);
	if (syUpper) printf("sy1[%d]=%f; ",i,syUpper[i]-y[i]);
	if (syLower) printf("sy2[%d]=%f; ",i,syLower[i]-y[i]);
      }
      printf("\n");
    }
  }
  if (e1y&&e2y) {
    float y_tmp[1000], ey_tmp[1000];
    for (int i=0; i<n; ++i) {
      y_tmp[i]=y[i]+(e1y[i]+e2y[i])/2;
      ey_tmp[i]=fabs(e1y[i]-e2y[i])/2;
    }
    TGraphErrors *graf = new TGraphErrors(n,x,y_tmp,ex,ey_tmp);
    //graf->SetMarkerStyle(1);
    //graf->SetMarkerSize(0.01);
    graf->SetMarkerColor(color);
    graf->SetLineColor(color);
    graf->Draw(" ");
    TGraphErrors *graf = new TGraphErrors(n,x,y,ex,0);
    if (symb>23) graf->SetMarkerStyle(symb-4);
    else graf->SetMarkerStyle(symb);
    graf->SetMarkerSize(sSize);
    graf->SetMarkerColor(0);
    graf->SetLineColor(color);
    graf->Draw("p");
    TGraphErrors *graf = new TGraphErrors(n,x,y,ex,0);
    graf->SetMarkerStyle(symb);
    graf->SetMarkerSize(sSize);
    graf->SetMarkerColor(color);
    graf->SetLineColor(color);
    graf->Draw("p");
  }
  else if (e1y) {
    TGraphErrors *graf = new TGraphErrors(n,x,y,ex,e1y);
    //graf->SetMarkerStyle(symb);
    //graf->SetMarkerSize(sSize);
    graf->SetMarkerColor(color);
    graf->SetLineColor(color);
    //graf->Draw("p");
    graf->Draw(" ");
    TGraphErrors *graf = new TGraphErrors(n,x,y,ex,0);
    if (symb>23) graf->SetMarkerStyle(symb-4);
    else graf->SetMarkerStyle(symb);
    graf->SetMarkerSize(sSize);
    graf->SetMarkerColor(0);
    graf->SetLineColor(color);
    graf->Draw("p");
    TGraphErrors *graf = new TGraphErrors(n,x,y,ex,0);
    graf->SetMarkerStyle(symb);
    graf->SetMarkerSize(sSize);
    graf->SetMarkerColor(color);
    graf->SetLineColor(color);
    graf->Draw("p");
  }
  if (!syUpper) return;
  float scale=0.01*sSize;
  if (scale<0.01) scale=0.01;
  /*
  float xmin=gPad->GetUxmin();
  float xmax=gPad->GetUxmax();
  float dx=(xmax-xmin)*scale;
  float ymin=gPad->GetUymin();
  float ymax=gPad->GetUymax();
  float dy=(ymax-ymin)*scale;
  */
  float dx=dxaxis*scale;
  float dy=dyaxis*scale;
  TLine *line=new TLine();
  line->SetLineColor(color);
  for (int i=0; i<n; ++i) {
    float x1=x[i]-dx;
    float x2=x[i]+dx;
    float y1=syUpper[i];
    if (syLower) float y2=syLower[i];
    else         float y2=y[i]-(syUpper[i]-y[i]);
    if (y1<y2) {float yy=y2; y2=y1; y1=yy;} 
    line->DrawLine(x1,y1,x2,y1);
    line->DrawLine(x1,y1,x1,y1-dy);
    line->DrawLine(x2,y1,x2,y1-dy);
    line->DrawLine(x1,y2,x2,y2);
    line->DrawLine(x1,y2,x1,y2+dy);
    line->DrawLine(x2,y2,x2,y2+dy);
  }
}

void graphSysErr3band(int n, Float_t *x, Float_t *y, 
		      Float_t *ex, Float_t *e1y, Float_t *e2y, 
		      Float_t *syUpper=0, Float_t *syLower=0, 
		      Float_t dxaxis=1, Int_t symb=20, Float_t sSize=1, 
		      Int_t color=1, Int_t colorband=5)
{ //same as the other graphSysErr3band except stat. errors are asymmetric.
  if (utilityPrint) printf("\n__________graphSysErr3band__________symb=%d sSize=%f color=%d colorband=%d__________\n",symb,sSize,color,colorband);
  if (utilityPrint) {
    for (int i=0; i<n; ++i) {
      printf("x[%d]=%f; y[%d]=%f; ",i,x[i],i,y[i]);
      if (utilityPrint<0) {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]/x[i]);
	if (e1y) printf("ey1[%d]=%f; ",i,e1y[i]/y[i]);
	if (e2y) printf("ey2[%d]=%f; ",i,e2y[i]/y[i]);
	if (syUpper) printf("sy1[%d]=%f; ",i,syUpper[i]/y[i]-1);
	if (syLower) printf("sy2[%d]=%f; ",i,syLower[i]/y[i]-1);
      }
      else {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]);
	if (e1y) printf("ey1[%d]=%f; ",i,e1y[i]);
	if (e2y) printf("ey2[%d]=%f; ",i,e2y[i]);
	if (syUpper) printf("sy1[%d]=%f; ",i,syUpper[i]-y[i]);
	if (syLower) printf("sy2[%d]=%f; ",i,syLower[i]-y[i]);
      }
      printf("\n");
    }
  }
  if (syUpper) {
    float scale=0.01*sSize;
    if (scale<0.01) scale=0.01;
    float dx=dxaxis*scale;
    float xband[1000], yband[1000];
    for (int i=0; i<n; ++i) {
      float xtmp=x[i];
      //if (i==0) xtmp-=dx;
      //if (i==n-1) xtmp+=dx;
      xband[i]=xtmp;
      xband[2*n-1-i]=xtmp;
      float y1=syUpper[i];
      if (syLower) float y2=syLower[i];
      else         float y2=y[i]-(syUpper[i]-y[i]);
      if (y1<y2) {float yy=y2; y2=y1; y1=yy;} 
      yband[i]=y1;
      yband[2*n-1-i]=y2;
      //printf("%d %d\n",i,2*n-1-i);
    }
    //for(int iii=0; iii<2*n;++iii)printf("%f %f\n",xband[iii],yband[iii]);
    TGraph *grafband = new TGraph(2*n,xband,yband);
    grafband->SetFillColor(colorband);
    grafband->Draw("f");
  }
  if (e1y&&e2y) {
    float y_tmp[1000], ey_tmp[1000];
    for (int i=0; i<n; ++i) {
      y_tmp[i]=y[i]+(e1y[i]+e2y[i])/2;
      ey_tmp[i]=fabs(e1y[i]-e2y[i])/2;
    }
    TGraphErrors *graf = new TGraphErrors(n,x,y_tmp,ex,ey_tmp);
    //graf->SetMarkerStyle(1);
    //graf->SetMarkerSize(0.01);
    graf->SetMarkerColor(color);
    graf->SetLineColor(color);
    graf->Draw(" ");
    TGraphErrors *graf = new TGraphErrors(n,x,y,ex,0);
    if (symb>23) graf->SetMarkerStyle(symb-4);
    else graf->SetMarkerStyle(symb);
    graf->SetMarkerSize(sSize);
    graf->SetMarkerColor(0);
    graf->SetLineColor(color);
    graf->Draw("p");
    TGraphErrors *graf = new TGraphErrors(n,x,y,ex,0);
    graf->SetMarkerStyle(symb);
    graf->SetMarkerSize(sSize);
    graf->SetMarkerColor(color);
    graf->SetLineColor(color);
    graf->Draw("p");
  }
  else (e1y) {
    TGraphErrors *graf = new TGraphErrors(n,x,y,ex,e1y);
    //graf->SetMarkerStyle(symb);
    //graf->SetMarkerSize(sSize);
    graf->SetMarkerColor(color);
    graf->SetLineColor(color);
    //graf->Draw("p");
    graf->Draw(" ");
    TGraphErrors *graf = new TGraphErrors(n,x,y,ex,0);
    if (symb>23) graf->SetMarkerStyle(symb-4);
    else graf->SetMarkerStyle(symb);
    graf->SetMarkerSize(sSize);
    graf->SetMarkerColor(0);
    graf->SetLineColor(color);
    graf->Draw("p");
    TGraphErrors *graf = new TGraphErrors(n,x,y,ex,0);
    graf->SetMarkerStyle(symb);
    graf->SetMarkerSize(sSize);
    graf->SetMarkerColor(color);
    graf->SetLineColor(color);
    graf->Draw("p");
    //printf("df %d\n\n",symb);
  }
}

/////////////////////////////////////////////////////////
void histSysErrBand(TH1D *h0, TH1D *h1, TH1D *h2, Int_t color=5, Int_t lcolor=0, Int_t lwidth=1, Float_t dedge=0.5, Int_t ich1=0, Int_t ich2=0, Int_t marker=0, Int_t mcolor=0)
{
  if (!ich1) ich1=1;
  if (!ich2) ich2=h1->GetNbinsX();
  for (int ich=ich1; ich<=ich2; ++ich) {
    float dedge1=0, dedge2=1;
    if (ich==ich1) dedge1=dedge;
    if (ich==ich2) dedge2=1-dedge;
    float grafxtmp[5], grafytmp[5];
    grafxtmp[0]=h1->GetBinLowEdge(ich)+dedge1*h1->GetBinWidth(ich);
    grafxtmp[1]=grafxtmp[0];
    grafxtmp[2]=h1->GetBinLowEdge(ich)+dedge2*h1->GetBinWidth(ich);
    grafxtmp[3]=grafxtmp[2];
    grafxtmp[4]=grafxtmp[0];
    grafytmp[0]=h1->GetBinContent(ich);
    grafytmp[1]=h2->GetBinContent(ich);
    grafytmp[2]=h2->GetBinContent(ich);
    grafytmp[3]=h1->GetBinContent(ich);
    grafytmp[4]=h1->GetBinContent(ich);
    TGraph *graftmp=new TGraph(5,grafxtmp,grafytmp);
    if(color){graftmp->SetFillColor(color); graftmp->Draw("f");}
    if (lwidth) {graftmp->SetLineWidth(lwidth);graftmp->SetLineColor(lcolor); graftmp->Draw("l");}
  }
  if(!marker||!h0)return;
  const int npts=ich2-ich1+1;
  float grafxpts[npts], grafypts[npts], grafepts[npts];
  for (int ich=ich1; ich<=ich2; ++ich) {
    grafxpts[ich-ich1]=h0->GetBinCenter(ich);
    grafypts[ich-ich1]=h0->GetBinContent(ich);
    grafepts[ich-ich1]=h0->GetBinError(ich);
  }
  TGraphErrors *grafpts=new TGraphErrors(npts,grafxpts,grafypts,0,grafepts);
  grafpts->SetMarkerColor(mcolor); grafpts->SetLineColor(mcolor);
  grafpts->SetMarkerStyle(marker); grafpts->Draw("p");
}
void histSysErrBand(TH1F *h0, TH1F *h1, TH1F *h2, Int_t color=5, Int_t lcolor=0, Int_t lwidth=1, Float_t dedge=0.5, Int_t ich1=0, Int_t ich2=0){histSysErrBand((TH1D*)h0,(TH1D*)h1,(TH1D*)h2,color,lcolor,lwidth,dedge,ich1,ich2);}

void histSysErrBand2(TH1D *h0, TH1D *h1, TH1D *h2, Int_t color=5, Int_t lcolor=0, Int_t lwidth=1, Int_t symb=20, Int_t scolor=2, Float_t dedge=0.5, Int_t ich1=0, Int_t ich2=0)
{
  h0->SetMarkerStyle(symb);
  //h0->SetMarkerSize(sSize);
  h0->SetMarkerColor(scolor);
  h0->SetLineColor(scolor);
  h0->Draw("e");
  histSysErrBand(h0, h1, h2, color, lcolor, lwidth, dedge, ich1, ich2);
  h0->Draw("esame");
}
void histSysErrBand(TH1F *h0, TH1F *h1, TH1F *h2, Int_t color=5, Int_t lcolor=0, Int_t lwidth=1, Int_t symb=20, Int_t scolor=2, Float_t dedge=0.5, Int_t ich1=0, Int_t ich2=0){histSysErrBand((TH1D*)h0, (TH1D*)h1, (TH1D*)h2, Int_t color=5, Int_t lcolor=0, Int_t lwidth=1, Int_t symb=20, Int_t scolor=2, Float_t dedge=0.5, Int_t ich1=0, Int_t ich2=0);}

void hist2SysErrBand(TH1D *h0, TH1D *h1u, TH1D *h1l, TH1D *h2u, TH1D *h2l, Int_t color=5, Int_t lcolor=0, Int_t lwidth=1, Float_t dedge=0.5, Int_t ich1=0, Int_t ich2=0)
{
  if (!ich1) ich1=1;
  if (!ich2) ich2=h0->GetNbinsX();
  for (int ich=ich1; ich<=ich2; ++ich) {
    float dedge1=0, dedge2=1;
    if (ich==ich1) dedge1=dedge;
    if (ich==ich2) dedge2=1-dedge;
    float grafxtmp[5], grafytmp[5];
    grafxtmp[0]=h0->GetBinLowEdge(ich)+dedge1*h0->GetBinWidth(ich);
    grafxtmp[1]=grafxtmp[0];
    grafxtmp[2]=h0->GetBinLowEdge(ich)+dedge2*h0->GetBinWidth(ich);
    grafxtmp[3]=grafxtmp[2];
    grafxtmp[4]=grafxtmp[0];
    float y0=h0->GetBinContent(ich);
    float y1u=h1u->GetBinContent(ich);
    float y1l=h1l->GetBinContent(ich);
    float y2u=h2u->GetBinContent(ich);
    float y2l=h2l->GetBinContent(ich);
    float d1u, d1l, d2u, d2l;
    if (y1u>=y0&&y1l<=y0) {d1u=y1u-y0; d1l=y1l-y0;}
    else if (y1u<=y0&&y1l>=y0) {d1u=y1l-y0; d1l=y1u-y0;}
    else {printf("%d y0=%f y1u=%f y1l=%f\n",ich,y0,y1u,y1l); d1u=0; d1l=0;}
    if (y2u>=y0&&y2l<=y0) {d2u=y2u-y0; d2l=y2l-y0;}
    else if (y2u<=y0&&y2l>=y0) {d2u=y2l-y0; d2l=y2u-y0;}
    else {printf("%d y0=%f y2u=%f y2l=%f\n",ich,y0,y2u,y2l); d2u=0; d2l=0;}
    float du= sqrt(d1u*d1u+d2u*d2u);
    float dl=-sqrt(d1l*d1l+d2l*d2l);
    //float du=d1u+d2u;
    //float dl=d1l+d2l;
    grafytmp[0]=y0+du;
    grafytmp[1]=y0+dl;
    grafytmp[2]=y0+dl;
    grafytmp[3]=y0+du;
    grafytmp[4]=y0+du;
    TGraph *graftmp=new TGraph(5,grafxtmp,grafytmp);
    graftmp->SetFillColor(color); graftmp->Draw("f");
    if (lcolor) {graftmp->SetLineWidth(lwidth); graftmp->Draw("l");}
  }
}
void hist2SysErrBand(TH1F *h0, TH1F *h1u, TH1F *h1l, TH1F *h2u, TH1F *h2l, Int_t color=5, Int_t lcolor=0, Int_t lwidth=1, Float_t dedge=0.5, Int_t ich1=0, Int_t ich2=0){hist2SysErrBand((TH1D*)h0, (TH1D*)h1u, (TH1D*)h1l, (TH1D*)h2u, (TH1D*)h2l, Int_t color=5, Int_t lcolor=0, Int_t lwidth=1, Float_t dedge=0.5, Int_t ich1=0, Int_t ich2=0);}
/////////////////////////////////////////////////////////

void printVector(int n, float *v1)
{
  printf("\n");
  for (int i=0; i<n; ++i) {
    printf("%d",i);
    printf(" %f",v1[i]);
    printf("\n");
  }
  printf("\n");
}

void printVector(int n, int *v1)
{
  printf("\n");
  for (int i=0; i<n; ++i) {
    printf("%d",i);
    printf(" %d",v1[i]);
    printf("\n");
  }
  printf("\n");
}

void printVector(int n, float *v1, float *v2)
{
  printf("\n");
  for (int i=0; i<n; ++i) {
    printf("%d",i);
    printf(" %f",v1[i]);
    printf(" %f",v2[i]);
    printf("\n");
  }
  printf("\n");
}

void printVector(int n, float *v1, float *v2, float *v3)
{
  printf("\n");
  for (int i=0; i<n; ++i) {
    printf("%d",i);
    printf(" %f",v1[i]);
    printf(" %f",v2[i]);
    printf(" %f",v3[i]);
    printf("\n");
  }
  printf("\n");
}

void printVector(int n, float *v1, float *v2, float *v3, float *v4)
{
  printf("\n");
  for (int i=0; i<n; ++i) {
    printf("%d",i);
    printf(" %f",v1[i]);
    printf(" %f",v2[i]);
    printf(" %f",v3[i]);
    printf(" %f",v4[i]);
    printf("\n");
  }
  printf("\n");
}

void printVector(int n, float *v1, float *v2, float *v3, float *v4, float *v5)
{
  printf("\n");
  for (int i=0; i<n; ++i) {
    printf("%d",i);
    printf(" %f",v1[i]);
    printf(" %f",v2[i]);
    printf(" %f",v3[i]);
    printf(" %f",v4[i]);
    printf(" %f",v5[i]);
    printf("\n");
  }
  printf("\n");
}

void printRatio(int n, float *v1, float *v2)
{
  printf("\n");
  for (int i=0; i<n; ++i) {
    printf("%d",i);
    printf(" %f",v1[i]);
    printf(" %f",v2[i]);
    printf(" %f",v2[i]/v1[i]);
    printf("\n");
  }
  printf("\n");
}

void graphLimits(float dx=0.2,float dy=0.2,float *xmin,float *xmax,float *ymin,float *ymax,TGraphErrors *g0=0,TGraphErrors *g1=0,TGraphErrors *g2=0,TGraphErrors *g3=0,TGraphErrors *g4=0,TGraphErrors *g5=0,TGraphErrors *g6=0,TGraphErrors *g7=0)
{
  graphLimits(dx,dx,dy,dy,xmin,xmax,ymin,ymax,g0,g1,g2,g3,g4,g5,g6,g7);
}

void graphLimits(float dx1=0.2,float dx2=0.2,float dy1=0.2,float dy2=0.2,float *xmin,float *xmax,float *ymin,float *ymax,TGraphErrors *g0=0,TGraphErrors *g1=0,TGraphErrors *g2=0,TGraphErrors *g3=0,TGraphErrors *g4=0,TGraphErrors *g5=0,TGraphErrors *g6=0,TGraphErrors *g7=0)
{
  const int Ng=8;
  TGraphErrors *g[Ng];
  g[0]=g0;g[1]=g1;g[2]=g2;g[3]=g3;g[4]=g4;g[5]=g5;g[6]=g6;g[7]=g7;
  double x,y,ex,ey;
  xmin[0]=999999;xmax[0]=-999999;
  ymin[0]=999999;ymax[0]=-999999;
  for(int i=0;i<Ng;++i){
    if(!g[i])continue;
    for(int k=0;k<g[i]->GetN();++k){
      g[i]->GetPoint(k,x,y);
      ex=g[i]->GetErrorX(k);
      ey=g[i]->GetErrorY(k);
      if(x-ex<xmin[0])xmin[0]=x-ex;
      if(x+ex>xmax[0])xmax[0]=x+ex;
      if(y-ey<ymin[0])ymin[0]=y-ey;
      if(y+ey>ymax[0])ymax[0]=y+ey;
    }
  }
  float dx1=dx1*(xmax[0]-xmin[0]),dx2=dx2*(xmax[0]-xmin[0]);xmin[0]-=dx1;xmax[0]+=dx2;
  float dy1=dy1*(ymax[0]-ymin[0]),dy2=dy2*(ymax[0]-ymin[0]);ymin[0]-=dy1;ymax[0]+=dy2;
  return;
}

void arrayLimits(float dy=0.2,float *ymin,float *ymax,int n=1,float *a0,float *e0p=0,float *e0m=0,float *a1=0,float *e1p=0,float *e1m=0,float *a2=0,float *e2p=0,float *e2m=0,float *a3=0,float *e3p=0,float *e3m=0,float *a4=0,float *e4p=0,float *e4m=0)
{
  arrayLimits(dy,dy,ymin,ymax,n,a0,e0p,e0m,a1,e1p,e1m,a2,e2p,e2m,a3,e3p,e3m,a4,e4p,e4m);
}

void arrayLimits(float dymin=0.2,float dymax=0.2,float *ymin,float *ymax,int n=1,float *a0,float *e0p=0,float *e0m=0,float *a1=0,float *e1p=0,float *e1m=0,float *a2=0,float *e2p=0,float *e2m=0,float *a3=0,float *e3p=0,float *e3m=0,float *a4=0,float *e4p=0,float *e4m=0)
{
  const int Na=5;
  float *a[Na],*ep[Na],*em[Na];
  a[0]=a0;ep[0]=e0p;em[0]=e0m;
  a[1]=a1;ep[1]=e1p;em[1]=e1m;
  a[2]=a2;ep[2]=e2p;em[2]=e2m;
  a[3]=a3;ep[3]=e3p;em[3]=e3m;
  a[4]=a4;ep[4]=e4p;em[4]=e4m;
  ymin[0]=999999;ymax[0]=-999999;
  for(int i=0;i<Na;++i){
    if(!a[i])continue;
    for(int k=0;k<n;++k){
      float e_p=0,e_m=0;
      if(ep[i])e_p=ep[i][k];
      if(em[i])e_m=fabs(em[i][k]);
      if(a[i][k]+e_p>ymax[0])ymax[0]=a[i][k]+e_p;
      if(a[i][k]-e_m<ymin[0])ymin[0]=a[i][k]-e_m;
    }
  }
  float dymin=dymin*(ymax[0]-ymin[0]),dymax=dymax*(ymax[0]-ymin[0]);
  ymin[0]-=dymin;ymax[0]+=dymax;
  return;
}
