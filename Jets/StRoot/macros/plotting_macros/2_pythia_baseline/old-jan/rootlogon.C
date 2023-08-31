{
//TStyle *gStyle= new TStyle("STAR_jan","Jan's style");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetOptDate(0);
  
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.09);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);
  
  gStyle->SetTitleSize(0.052,"Y");
  gStyle->SetTitleOffset(0.95,"Y");
  gStyle->SetTitleSize(0.06,"X");
  gStyle->SetTitleOffset(0.95,"X");
  gStyle->SetLabelSize(0.03,"X");
  gStyle->SetLabelSize(0.03,"Y");
  
  gStyle->SetLineWidth(0);
  gStyle->SetPadTickY(1);

 
//for png, gif graphics 
	const Int_t hatchesLineWidth_PNG=1;
	const Int_t lineLineWidthLine_PNG=2;
	const Int_t graphLineWidth_PNG=2;
	const Int_t fillStyle_PNG=3001;
	const Int_t gridStyle_PNG=3;
	const Color_t gridColor_PNG=kGray+1;
	const Int_t gridLineWidth_PNG=0;

//for pdf, ps graphics
	const Int_t hatchesLineWidth_PDF=0;
	const Int_t lineLineWidthLine_PDF=2;
	const Int_t graphLineWidth_PDF=1;
	const Int_t fillStyle_PDF=3002;//3244;
	const Int_t gridStyle_PDF=3;
	const Color_t gridColor_PDF=kGray;
	const Int_t gridLineWidth_PDF=0;
	
	
	const Int_t hatchesLineWidthG[]={hatchesLineWidth_PNG,hatchesLineWidth_PDF};
	const Int_t lineLineWidthG[]={lineLineWidthLine_PNG,lineLineWidthLine_PDF};
	const Int_t graphLineWidthG[]={graphLineWidth_PNG,graphLineWidth_PDF};
	const Int_t fillStyleG[]={fillStyle_PNG,fillStyle_PDF};
	const Int_t gridStyleG[]={gridStyle_PNG,gridStyle_PDF};
	const Color_t gridColorG[]={gridColor_PNG,gridColor_PDF};
	const Int_t gridLineWidthG[]={gridLineWidth_PNG,gridLineWidth_PDF};
}
