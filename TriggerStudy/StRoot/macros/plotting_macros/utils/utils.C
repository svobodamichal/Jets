#include "common.h"
#include "fit_functions.C" //declarations of fitting functions

typedef unsigned short ushort;
typedef TGraphAsymmErrors Graph;
enum ErrorType { down=-1, none, up };


//================================================================================
//calculate index of given R or pTlead in the C_nR[] or C_pTls[] arrays
//================================================================================
short calculate_index(float x, ushort idxtype /*0 - R, 1 - pTleading*/)
{
	float eps=0.001;
	int idx=0;
	
	if(idxtype==0) //Rindex
	for(int r=0;r<C_nR; r++)
	{
		if(x-eps<C_nR[r] && x+eps>C_nR[r]) idx=r;
	}
	else if(idxtype==1) //pTlead index
	for(int p=0;p<C_npTlead; p++)
	{
		if(x-eps<C_pTls[p] && x+eps>C_pTls[p]) idx=p;
	}
	return idx;
}

//================================================================================
//calculate number of bins of given histogram between values pTmin and pTmax and save the position of the first (nonzero) bin above pTmin
//================================================================================
ushort calc_nbins(TH1* histo, float pTmin, float pTmax, short &firstbin, bool nonzero=true /*start with nonzero bin*/)
{
	short nbins=0;
	firstbin=0;
	for(int bn=1; bn<=histo->GetNbinsX(); bn++)
	{
		float center=histo->GetBinCenter(bn);
		float val=histo->GetBinContent(bn);
		//cout<<"x:"<<center<<" y:"<<val<<endl;
		if(center<pTmin || center>pTmax) continue;
		if(firstbin==0)
		{
			//cout<<"firstbin="<<bn;
			if(nonzero && val<1E-20)continue;
			firstbin=bn;
			//cout<<"firstbin="<<bn;
		}
		nbins++;
	}
	return nbins;
}

//================================================================================
//calculate error of a ratio through a MC propagation
//================================================================================
void calc_sigma_ratio(Double_t mu1, Double_t mu2, Double_t sigma1, Double_t sigma2, Double_t &sigma_up, Double_t &sigma_down, Double_t &center, Int_t N=100, bool useMPV=1, TH1D* hratio=NULL, TF1* gauss1=NULL,TF1* gauss2=NULL,float left=0, float right=0)
{

   if(mu1>0 && mu2>0){
   //cout<<"calculating ratio error"<<endl;
	//hratio=new TH1D("hratio","",100,left,right);
   
   for(int i=0; i<N; i++)
   {
      Double_t rnd1=gRandom->Gaus(mu1,sigma1);
      Double_t rnd2=gRandom->Gaus(mu2,sigma2);
      Double_t xi=0;
      xi=(TMath::Abs(rnd2)>0) ? rnd1/rnd2 : 0;
      hratio->Fill(xi);
   }
   double mu=hratio->GetMean(); //mean value
	double mpv=hratio->GetBinCenter(hratio->GetMaximumBin()); //most probable value
   
	center=(useMPV) ? mpv : mu; //use most probable value instead of mean - MPV is more proper choice and will lead to more accurate fits
   //fit

	gauss1->SetParameter(0,N);
	gauss1->FixParameter(1,center);
	gauss1->SetParameter(2,TMath::Abs(right-center)/2);
	gauss1->SetLineColor(kRed);

	gauss2->SetParameter(0,N);
	gauss2->FixParameter(1,center);
	gauss2->SetParameter(2,TMath::Abs(center-left)/2);
	gauss2->SetLineColor(kMagenta);

	hratio->Fit(gauss1,"","",center,right);
	hratio->Fit(gauss2,"","",left,center);  

	sigma_up=TMath::Abs(gauss1->GetParameter(2));
	sigma_down=TMath::Abs(gauss2->GetParameter(2));

	//too large sigmas (>sigmaMAX) come from bad fit => we will replace them
	const float sigmaMAX=1.0;  //cut value
	const float sigmaMIN=0.0001;  //cut value
	if(sigma_up>sigmaMAX || sigma_up<sigmaMIN) sigma_up=(sigma1/mu1+sigma2/mu2)*mu1/mu2;
	if(sigma_down>sigmaMAX || sigma_down<sigmaMIN) sigma_down=(sigma1/mu1+sigma2/mu2)*mu1/mu2;

	//cout<<"sigma_up: "<<sigma_up<<"sigma_down: "<<sigma_down<<endl;
	return;
   }
   
   else
   {
      sigma_down=0;
      sigma_up=0;
      return;
   }
}

//================================================================================
//REBIN histogram "hIN" to the same binning as "hTEMPLATE" and give it the name "name".
//If the input histogram is normalized by the bin width, set input_scaled_to_width=1, otherwise use input_scaled_to_width=0
//================================================================================
TH1D* rebinhisto(TH1D* hIN, TH1D* hTEMPLATE, TString name="", bool input_scaled_to_width=1)
{
	TH1D* houtput=(TH1D*) hTEMPLATE->Clone(name);
	houtput->Reset();
	int nbins_old=hIN->GetNbinsX();
	int nbins_new=houtput->GetNbinsX();
	
	for(int i=1; i<=nbins_old; i++)
	{
		double low_edge_old=hIN->GetBinLowEdge(i);
		double width_old=hIN->GetBinWidth(i);
		double upper_edge_old=low_edge_old+width_old;
		double bin_content_old=hIN->GetBinContent(i);
		double bin_error_old=hIN->GetBinError(i);
		if(input_scaled_to_width) //make sure we have counts, not counts/bin_width
		{
			bin_content_old=bin_content_old*width_old; 
			bin_error_old=bin_error_old*width_old; 
		}
		double bin_error_old2=bin_error_old*bin_error_old;
		
		for(int j=1; j<=nbins_new; j++)
		{
			bool overlap=false;
			double low_edge_new=houtput->GetBinLowEdge(j);
			double width_new=houtput->GetBinWidth(j);
			double upper_edge_new=low_edge_new+width_new;
			double bin_content_new=houtput->GetBinContent(j);
			double bin_error_new=houtput->GetBinError(j);
			double bin_error_new2=bin_error_new*bin_error_new;
			
			/*calculate new bin content and bin error of the new bins
			the situation may look like this:
			                a)             b)            c.1)           c.2)
			old:          |   |          |   |          |   |          |   |
			new:           ||           |      |      |   |               |   |
			*/
			//a) new bin lies completely within the old bin
			if(low_edge_new>=low_edge_old && upper_edge_new<=upper_edge_old) 
			{
				bin_content_new+=bin_content_old*width_new/width_old;
				bin_error_new2+=bin_error_old2*width_new/width_old; // because sigma_i^2/sigma_j^2 = n_i/n_j and sigma_tot^2=sigma_1^2+sigma_2^2+...
				overlap=true;
			}
			//b) old bin lies completely within the new bin
			else if(low_edge_new<=low_edge_old && upper_edge_new>=upper_edge_old) 
			{
				bin_content_new+=bin_content_old;
				bin_error_new2+=bin_error_old2; 
				overlap=true;
			}
			//c.1) the new bin and old bin have a partial overlap - new one starts below the old one's lower edge
			else if(low_edge_new<=low_edge_old && upper_edge_new>low_edge_old) 
			{
				double overlap_width=upper_edge_new-low_edge_old;
				bin_content_new+=bin_content_old*overlap_width/width_old;
				bin_error_new2+=bin_error_old2*overlap_width/width_old; 
				overlap=true;
			}
			//c.2) the new bin and old bin have a partial overlap - new one starts above the old one's lower edge
			else if(low_edge_new<upper_edge_old && upper_edge_new>=upper_edge_old) 
			{
				double overlap_width=upper_edge_old-low_edge_new;
				bin_content_new+=bin_content_old*overlap_width/width_old;
				bin_error_new2+=bin_error_old2*overlap_width/width_old; 
				overlap=true;
			}
			
			if(!overlap) continue; //this bin has nothing in common with the old bin
			
			//Fill the new bin content...
			houtput->SetBinContent(j,bin_content_new);
			//...and error
			if(bin_error_new2>0) bin_error_new=TMath::Sqrt(bin_error_new2);
			else bin_error_new=0;
			houtput->SetBinError(j,bin_error_new);
				
		}//loop over new bins
	}//loop over old bins
	
	if(input_scaled_to_width) //if desired normalize the histogram by the bin width
		houtput->Sumw2();
		houtput->Scale(1.0,"width");
	return houtput;
}

/*{  
	TH1D* hOLD=(TH1D*) hIN->Clone("hold");
	hOLD->Sumw2();
	const int nbns=401;
	for(int i=0; i<4; i++)
	{
		if(hOLD->GetNbinsX()<nbns)continue;
		hOLD->Rebin(2);
		if(input_scaled_to_width)hOLD->Scale(0.5);
	}
	Double_t yield[nbns];
	Double_t error[nbns];
	Double_t binwidth[nbns];
	TH1D* hNEW=(TH1D*) hTEMPLATE->Clone(name.Data());
	hNEW->Reset("MICE");
	hNEW->Sumw2();
	for(int bn_new=1;bn_new<hNEW->GetNbinsX()+1;bn_new++){
		int bincount=0;
		for(int bn_old=1; bn_old<hOLD->GetNbinsX()+1;bn_old++)
		{
			if(hNEW->FindBin(hOLD->GetBinCenter(bn_old))!=bn_new)continue;
			//cout<<"bin new "<<bn_new<<" pT: "<<hNEW->GetBinCenter(bn_new)<<"bin old "<<bn_old<<" pT: "<<hOLD->GetBinCenter(bn_old)<<endl;
			yield[bincount]=hOLD->GetBinContent(bn_old);
			error[bincount]=hOLD->GetBinError(bn_old);
			binwidth[bincount]=hOLD->GetBinWidth(bn_old);
			if(input_scaled_to_width)
			{
				yield[bincount]=yield[bincount]*binwidth[bincount];
				error[bincount]=error[bincount]*binwidth[bincount];
			}
			bincount++;
		}//loop over bins of original histogram
		if(bincount==0)continue;

		double yieldt=0;
		double errt2=0;
		float bwidth=0;
		for(int i=0;i<bincount;i++)
		{
			yieldt+=yield[i];
			errt2+=error[i]*error[i];
			bwidth+=binwidth[i];
		}

		if(input_scaled_to_width)
		{
			hNEW->SetBinContent(bn_new,yieldt/bwidth);
			hNEW->SetBinError(bn_new,TMath::Sqrt(errt2)/bwidth);
		}
		else
		{
			hNEW->SetBinContent(bn_new,yieldt);
			hNEW->SetBinError(bn_new,TMath::Sqrt(errt2));
		}
	}//loop over bins of new histogram 

	delete hOLD;
	return hNEW;
} */

//================================================================================
//BIN CENTER CORRECTION (fit) - move center of the bin to the position of the mean pT (= move it left for falling spectra)
//the new center is calculated based on a fit to the spectrum
//works with arrays wich can be then used for a TGraph
//================================================================================
void bin_shift_LR(TF1* FitF, int nbins, const double* xarr, const double* yarr, const double* xerr_low=0, const double* xerr_high=0, const double* yerr_low=0, const double* yerr_high=0,const int niter=5)
{
  cout<<"BIN CENTER CORRECTION"<<endl;
   for(int itr=1;itr<niter;itr++)
   {
   TGraphAsymmErrors* TGE=new TGraphAsymmErrors(nbins,xarr,yarr,xerr_low,xerr_high,yerr_low,yerr_high);
   TGE->Fit(FitF);
   
   
	//loop over points in graph
	for(int i=0; i<nbins; i++) 
	{
	Double_t xleft=xarr[i]-xerr_low[i];
	Double_t xright=xarr[i]+xerr_high[i];
   
   Double_t int_left=FitF->Integral(xleft,xarr[i]);
   Double_t int_right=FitF->Integral(xarr[i],xright);
   double int_ratio=int_right/int_left;
      if(itr==niter-1) cout<<"bin"<<i<<" integral ratio:"<<int_ratio<<endl;
   xarr[i]=xleft+(xarr[i]-xleft)*int_ratio;
   xerr_low[i]=xarr[i]-xleft;
   xerr_high[i]=xright-xarr[i];
   }//loop over points in graph

   }//loop over iterations
	return;
}
//================================================================================
//BIN CENTER CORRECTION (interpolation) - move center of the bin to the position of the mean pT (= move it left for falling spectra)
//the new center is calculated based on a linear interpolationbetween the points of the spectrum
//works with arrays wich can be then used for a TGraph
//================================================================================
void bin_shift_LR(int nbins, const double* xarr, const double* yarr, const double* xerr_low=0, const double* xerr_high=0, const double* yerr_low=0, const double* yerr_high=0,const int niter=2)
{
	cout<<"BIN CENTER CORRECTION"<<endl;

	//iteration loop
	for(int iter=0; iter<niter; iter++){
	
	//bin loop
	for(int bin=1;bin<nbins-1;bin++)
	{
		double x1,y1,x2,y2,x3,y3,dxL,dyL,dxR,dyR,w,xL,xR,xLnew,xRnew,yPartL,yPartR, k;
		//x,y values in bins i-1, i, i+1
		x1=xarr[bin-1];
		y1=yarr[bin-1];
		x2=xarr[bin];
		y2=yarr[bin];
		x3=xarr[bin+1];
		y3=yarr[bin+1];
		
		//left and right x-error bar in ith bin (=half of the bin width in the first iteration)
		xL=xerr_low[bin];
		xR=xerr_high[bin];
		//bin width
		w=xL+xR;
		
		//diferences between bin centers and bin contents
		dyL=TMath::Abs(y1-y2);
		dyR=TMath::Abs(y2-y3);
		dxL=TMath::Abs(x2-x1);
		dxR=TMath::Abs(x3-x2);
		
		//part of dy which "belongs" to the ith bin
		yPartL=(dxL>0) ? ((dyL/dxL)*xL) : 0;
		yPartR=(dxR>0) ? ((dyR/dxR)*xR) : 0;
		//their ratio
		k=(yPartR>0) ? (yPartL/yPartR) : 0;
		
		//new x-errorbars and new bin center
		xRnew=(k*w)/(1.0+k);
		xLnew=w-xRnew;
		
		//update the arrays
		xarr[bin]=x2-xL+xLnew;
		xerr_low[bin]=xLnew;
		xerr_high[bin]=xRnew;
	}//bin loop
	}//iteration loop
	return;
}

//================================================================================
//BIN CENTER CORRECTION (up and down) - move center of the bin up or down so a smooth function going through the new centers has the same integral as the original bins
//works with a histogram
//================================================================================
TH1D* bin_shift_UD_hist(TF1* FitF, TH1D* hIN,TString outname="", const int niter=5, const float alpha=0.5)
{
   if(outname=="")outname=Form("%s_%s",hIN->GetTitle(),"recentered");
   TH1D* hNEW=(TH1D*) hIN->Clone(outname);
   for(int itr=1;itr<=niter;itr++)
   {  
      TH1D* htemp=(TH1D*)hNEW->Clone("htemp");
      hNEW->Reset("MICE");
      htemp->Fit(FitF);
      int nbins=htemp->GetNbinsX();
      for(int i=1; i<=nbins; i++) 
      {
         double yIN=hIN->GetBinContent(i);
         double yNEW=htemp->GetBinContent(i);
         if(!(yNEW>0))continue;
         double error=htemp->GetBinError(i);
         double width=htemp->GetBinWidth(i);
         double xleft=htemp->GetBinLowEdge(i);
         double xright=xleft+width;
         double xcenter=htemp->GetBinCenter(i);
         
         double intH=yIN*width;
         //double intHNEW=yIN*width;
         double intF=FitF->Integral(xleft,xright);
         
         //scale y so bin center matches the fit function (just in case the fit was unable to go through all bin centers)
         double yf=FitF->Eval(xcenter);
         double prescale=1+(yf/yNEW-1)*alpha;
         yNEW=yNEW*prescale;
         
         //scale y up or down
         double scale=1+(intH/intF-1)*alpha;
         yNEW=yNEW*scale;
         error=error*scale;
         hNEW->SetBinContent(i,yNEW);
         hNEW->SetBinError(i,error);
         
        // if(xcenter>30 && xcenter<40)
		//cout<<"iter: "<<itr<<" original y: "<<yIN<<" new y: "<<yNEW<<" I f: "<<intF<<" I h: "<<intH<<" ratio H/F: "<<intH/intF<<" scale: "<<scale<<endl;
      }//bin loop
      delete htemp;
   }//iteration loop
   return hNEW;
}


//================================================================================
//Calculate average and spread between several histograms (passed as an array of pointers) for each bin.
//================================================================================
void calculate_avrg_spread(unsigned int nhistos, TH1D** histoArray, TH1D* haverage, TH1D* hsigmap, TH1D* hsigman, unsigned short divideByN=0, bool calculate_avg=true, short first_histo=0)
{
	//haverage, hsigmap, hsigman have to be empty histograms (except the case when calculate_avg=false) with the same binning as histograms in histoArray
	//calculate_avg=true: calculate the average from the input array
	//calculate_avg=false: use the content of the "TH1D *haverage" as the average value
	//first_histo: first histogram of the array from which we want to start with the calculation - use values >0 if you want to skip some of the histos in the array 
	//   (e.g.: we want to skip the first one if it is the average which we want to use as our external average by setting calculate_avg=false)
	//divideByN: divide the spread (=sqrt(Sum_i delta^2_i)) with x where 0:x=1, 1:x=sqrt(N), 2:x=sqrt(N*(N-1))
	
	if(nhistos==0) return;
	
	//calculate average (otherwise use the external one)
	if(calculate_avg){
		haverage->Reset("MICE");
		for(int bin=1; bin<haverage->GetNbinsX();bin++)
		{
			float avg=0;
			for(int i=first_histo; i<nhistos+first_histo; i++)
			{
				float val=histoArray[i]->GetBinContent(bin);
				avg+=val;
				//cout<<"bin "<<bin<<" sum:"<<avg<<endl;
			}
			avg=(float)avg/nhistos;
			haverage->SetBinContent(bin,avg);
		
			//set statistical error (from the first histogram of the histoArray
			float serror=(histoArray[0]->GetBinContent(bin)>0) ? histoArray[0]->GetBinError(bin)/histoArray[0]->GetBinContent(bin) : 0;
			serror=serror*haverage->GetBinContent(bin);
			haverage->SetBinError(bin,serror);
		}
	}

	for(int bin=1; bin<haverage->GetNbinsX();bin++)
	{
		float sigma_plus=0;
		float sigma_minus=0;
		float nplus=0;
		float nminus=0;
		float val1=haverage->GetBinContent(bin);
		if(!(val1>0))continue;
		for(int i=first_histo; i<nhistos+first_histo; i++)
		{
			float val2=histoArray[i]->GetBinContent(bin);
			if(!(val2>0))continue;
			float dif=val2-val1;
			if(dif>0) 
			{
				sigma_plus+=(dif*dif);
				nplus++;
			}
			else
			{
				sigma_minus+=(dif*dif);
				nminus++;
			}
		}//histogram list loop
		sigma_plus=TMath::Sqrt(sigma_plus);
		sigma_minus=TMath::Sqrt(sigma_minus);
		
		//normalization of the error
		float normp, normn;
		switch (divideByN)
		{
			case 1 : 
				normp=(nplus>0) ? (TMath::Sqrt(nplus)) : 1;
				normn=(nminus>0) ? (TMath::Sqrt(nminus)) : 1;
				break;
			case 2 :
				normp=((nplus-1)>0) ? (TMath::Sqrt(nplus*(nplus-1))) : 1;
				normn=((nminus-1)>0) ? (TMath::Sqrt(nminus*(nminus-1))) : 1;
				break;
			default :
				normp=1;
				normn=1;
		}
		
		//cout<<"sigma before normaliztaion:"<<sigma_plus<<endl;
		sigma_plus=sigma_plus/normp;
		sigma_minus=sigma_minus/normn;
		//cout<<"sigma after normaliztaion:"<<sigma_plus<<endl;
		
		hsigmap->SetBinContent(bin,sigma_plus);
		hsigman->SetBinContent(bin,sigma_minus);
		
		
		
	}//bin loop
	return;
}

//================================================================================
//Calculate average and spread between several histograms (passed as an array of pointers) for each bin.
//NOTE: The array of histograms is now a 2D array (histo[i][j]) instead of 1D as in case of "calculate_avrg_spread()"
//================================================================================
void calculate_avrg_spread2D(unsigned int nhistos1, unsigned int nhistos2, TH1D* histoArray[52][52], TH1D* haverage, TH1D* hsigmap, TH1D* hsigman, unsigned short divideByN=0, bool calculate_avg=true, short first_histo=0)
{
	//haverage, hsigmap, hsigman have to be empty histograms (except the case when calculate_avg=false) with the same binning as histograms in histoArray
	//calculate_avg=true: calculate the average from the input array
	//calculate_avg=false: use the content of the "TH1D *haverage" as the average value
	//first_histo: first histogram of the array from which we want to start with the calculation - use values >0 if you want to skip some of the histos in the array 
	//   (e.g.: we want to skip the first one if it is the average which we want to use as our external average by setting calculate_avg=false)
	//divideByN: divide the spread (=sqrt(Sum_i delta^2_i)) with x where 0:x=1, 1:x=sqrt(N), 2:x=sqrt(N*(N-1))
	int nhistos=nhistos1*nhistos2;
	if(nhistos==0) return;
	
	//calculate average (otherwise use the external one)
	if(calculate_avg){
		haverage->Reset("MICE");
		for(int bin=1; bin<haverage->GetNbinsX();bin++)
		{
			float avg=0;
			for(int i=0; i<nhistos1; i++){
			for(int j=0; j<nhistos2; j++)
			{
				if(i==0 && j<first_histo)continue;
				float val=histoArray[i][j]->GetBinContent(bin);
				avg+=val;
				//cout<<"bin "<<bin<<" sum:"<<avg<<endl;
			}}
			avg=(float)avg/nhistos;
			haverage->SetBinContent(bin,avg);
		
			//set statistical error (from the first histogram of the histoArray
			float serror=(histoArray[0][0]->GetBinContent(bin)>0) ? histoArray[0][0]->GetBinError(bin)/histoArray[0][0]->GetBinContent(bin) : 0;
			serror=serror*haverage->GetBinContent(bin);
			haverage->SetBinError(bin,serror);
		}
	}

	for(int bin=1; bin<haverage->GetNbinsX();bin++)
	{
		float sigma_plus=0;
		float sigma_minus=0;
		float nplus=0;
		float nminus=0;
		float val1=haverage->GetBinContent(bin);
		if(!(val1>0))continue;
		
		for(int i=0; i<nhistos1; i++){
		for(int j=0; j<nhistos2; j++)
		{
			float val2=histoArray[i][j]->GetBinContent(bin);
			if(!(val2>0))continue;
			float dif=val2-val1;
			if(dif>0) 
			{
				sigma_plus+=(dif*dif);
				nplus++;
			}
			else
			{
				sigma_minus+=(dif*dif);
				nminus++;
			}
		}}//histogram list loop
		sigma_plus=TMath::Sqrt(sigma_plus);
		sigma_minus=TMath::Sqrt(sigma_minus);
		
		//normalization of the error
		float normp, normn;
		switch (divideByN)
		{
			case 1 : 
				normp=(nplus>0) ? (TMath::Sqrt(nplus)) : 1;
				normn=(nminus>0) ? (TMath::Sqrt(nminus)) : 1;
				break;
			case 2 :
				normp=((nplus-1)>0) ? (TMath::Sqrt(nplus*(nplus-1))) : 1;
				normn=((nminus-1)>0) ? (TMath::Sqrt(nminus*(nminus-1))) : 1;
				break;
			default :
				normp=1;
				normn=1;
		}
		
		//cout<<"sigma before normaliztaion:"<<sigma_plus<<endl;
		sigma_plus=sigma_plus/normp;
		sigma_minus=sigma_minus/normn;
		//cout<<"sigma after normaliztaion:"<<sigma_plus<<endl;
		
		hsigmap->SetBinContent(bin,sigma_plus);
		hsigman->SetBinContent(bin,sigma_minus);
		
		
		
	}//bin loop
	return;
}
//================================================================================
//interpolation function (x)
//given the y-coordinate "yvalue", look at the TGraph "graph" and find corresponding x-coordinate (by interpolating two closest points around the "yval")
//if "errt" is set to "up" or "down", corresponding error (readed from Graph "graph_err" - so we can use several kinds of errors) is taken into account
//CAUTION: We assume a monotonous shape of the graph!
//================================================================================

double interpolate_x(Graph* graph, double yval, ErrorType errt=none, Graph* graph_err=NULL, bool falling=1/*falling or rising spectrum*/,bool relative_err=0/*do we want to use the relative size of the error or its absolute value?*/, bool debug=0)
{
		if(debug) cout<<"interpolating for y="<<yval<<endl;
		//first, find the corresponding 2 points of the 2nd spectrum between which we want to interpolate
		double xss1,xss2,yss1,yss2;
		int pss1=0;
		int pss2=0;
		int npoints=graph->GetN();
		if(debug) cout<<"n points:"<<npoints<<endl;
		for(int point=0; point<npoints; point++)
		{
			//if(debug) cout<<"point #"<<point<<endl;
			double xtmp,xtmp2;
			double ytmp,ytmp2;
			graph->GetPoint(point,xtmp,ytmp);
			if(relative_err && errt!=none)//we want to use the relative size of the error, not its absolute value)
			{
				graph_err->GetPoint(point,xtmp2,ytmp2);
			}
			if(errt==up)
			{
				float yerr=graph_err->GetErrorYhigh(point);
				if(relative_err) //we want to use the relative size of the error, not its absolute value)
					yerr=(ytmp2>0) ? (yerr/ytmp2)*ytmp : 0;
				ytmp=ytmp+yerr;
			}
			else if(errt==down)
			{
				float yerr=graph_err->GetErrorYlow(point);
				if(relative_err) //we want to use the relative size of the error, not its absolute value)
					yerr=(ytmp2>0) ? (yerr/ytmp2)*ytmp : 0;
				ytmp=ytmp-yerr;
			}
			if(debug) cout<<"ytmp:"<<ytmp<<" xtmp:"<<xtmp<<endl;
			if((ytmp<yval && falling) || (ytmp>yval && !falling))  //we are already bellow/above yval => we found the 2nd point (and the first point is hence the previous one which is saved in *ss1 variables)
			{
				if(point==0) return xtmp; //yval is higher/lower than the first point of the graph!
				xss2=xtmp;
				yss2=ytmp;
				pss2=point;

				if(debug)
				{
					cout<<"edge points (x,y):"<<endl;
					cout<<"("<<xss1<<","<<yss1<<") - ("<<xss2<<","<<yss2<<")"<<endl;
				}
				break;
			}
			else
			{
				if(point==npoints-1) return xtmp;//yval is lower/higher than the last point of the graph!
				xss1=xtmp;
				yss1=ytmp;
				pss1=point;
			}
		}//loop over points in 2nd spectrum
	

		double deltay=TMath::Abs(yss1-yss2);
		double deltax=TMath::Abs(xss1-xss2);
		double xval_lin;
		double fraction;

		//this should work both for falling and rising spectrum
		fraction=(deltay>0) ? (TMath::Abs(yval-yss1)/deltay) : 0;
		xval_lin=xss1+(fraction*(xss2-xss1));
		

		return xval_lin;
}
//================================================================================
//interpolation function (y)
//given the x-coordinate "xvalue", look at the TGraph "graph" and find corresponding y-coordinate (by interpolating two closest points around the "xval")
//if "errt" is set to "up" or "down", corresponding error (readed from Graph "graph_err" - so we can use several kinds of errors) is taken into account
//================================================================================

double interpolate_y(Graph* graph, double xval, ErrorType errt=none, Graph* graph_err=NULL)
{
		//first, find the corresponding 2 points of the 2nd spectrum between which we want to interpolate
		double xss1,xss2,yss1,yss2;
		int pss1=0;
		int pss2=0;
		int npoints=graph->GetN();
		for(int point=0; point<npoints; point++)
		{
			double xtmp;
			double ytmp;
			graph->GetPoint(point,xtmp,ytmp);
			if(errt==up)
			{
				float xerr=graph_err->GetErrorXhigh(point);
				xtmp=xtmp+xerr;
			}
			else if(errt==down)
			{
				float xerr=graph_err->GetErrorXlow(point);
				xtmp=xtmp-xerr;
			}
			
			if(xtmp>xval)  //we are already above xval => we found the 2nd point (and the first point is hence the previous one which is saved in *ss1 variables)
			{
				if(point==0) return ytmp; //xval is lower than the first point of the graph!
				xss2=xtmp;
				yss2=ytmp;
				pss2=point;

				//cout<<"edge points (x,y):"<<endl;
				//cout<<"("<<xss1<<","<<yss1<<") - ("<<xss2<<","<<yss2<<")"<<endl;
				break;
			}
			else
			{
				if(point==npoints-1) return ytmp;//x val is higher than the last point of the graph!
				xss1=xtmp;
				yss1=ytmp;
				pss1=point;
			}
		}//loop over points in 2nd spectrum

		double deltay=TMath::Abs(yss1-yss2);
		double deltax=TMath::Abs(xss1-xss2);
		double yval_lin;
		double fraction;

		//this should work both for falling and rising spectrum
		fraction=(deltax>0) ? ((xval-xss1)/deltax) : 0;
		yval_lin=yss1+(fraction*(yss2-yss1));
		
		return yval_lin;
}
//================================================================================
//create a TGraphAsymmErrors from a histogram
//================================================================================

Graph* histo2graph(TH1D* histo)
{
	int ntmp=histo->GetNbinsX();
	const int nbins=ntmp;
	double x[nbins];
	double y[nbins];
	double xl[nbins];
	double xr[nbins];
	double yl[nbins];
	double yr[nbins];
	
	for(int bin=0; bin<nbins; bin++)
	{
		x[bin]=histo->GetBinCenter(bin+1);
		y[bin]=histo->GetBinContent(bin+1);
		xl[bin]=histo->GetBinWidth(bin+1)/2.0;
		xr[bin]=histo->GetBinWidth(bin+1)/2.0;
		yl[bin]=histo->GetBinError(bin+1);
		yr[bin]=histo->GetBinError(bin+1);
	}
	Graph* graph=new TGraphAsymmErrors(nbins, x, y, xl, xr, yl, yr);
	return graph;
}

//================================================================================
//move center points of the graph "grtomove" to the new positions given by the points from "grtemplate" and rescale the y-errors accordingly
//================================================================================
void move_graph_points(Graph* grtomove, Graph* grtemplate)
{
	for(int pt=0; pt<grtemplate->GetN(); pt++)
	{
		double xtemp,xtemp_down,xtemp_up,ytemp;
		double xmove,xmove_down,xmove_up,ymove,ymove_up,ymove_down;

		grtemplate->GetPoint(pt,xtemp,ytemp);
		grtomove->GetPoint(pt,xmove,ymove);
		double scaler=(TMath::Abs(ymove>0)) ? ytemp/ymove : 0;
		
	   xtemp_down=grtemplate->GetErrorXlow(pt);
	   xtemp_up=grtemplate->GetErrorXhigh(pt);

	   xmove_down=grtomove->GetErrorXlow(pt);
	   xmove_up=grtomove->GetErrorXhigh(pt);
	   ymove_down=grtomove->GetErrorYlow(pt);
	   ymove_up=grtomove->GetErrorYhigh(pt);
		
		//cout<<pt<<" scaler:"<<scaler<<" yold:"<<ymove_down<<" ynew:"<<ymove_down*scaler<<endl;
		
		grtomove->SetPoint(pt, xtemp, ytemp);
		grtomove->SetPointError(pt, xtemp_down, xtemp_up, ymove_down*scaler, ymove_up*scaler);

	}
	return;
}

//================================================================================
//scale every point of the graph and corresponding y errors by a constant
//================================================================================
void scale_graph(Graph* graph, double scaler, bool divide=0)
{
	for(int pt=0; pt<graph->GetN(); pt++)
	{
		double x,x_down,x_up,y,y_up,y_down;

		graph->GetPoint(pt,x,y);
	   x_down=graph->GetErrorXlow(pt);
	   x_up=graph->GetErrorXhigh(pt);
	   y_down=graph->GetErrorYlow(pt);
	   y_up=graph->GetErrorYhigh(pt);

       if(!divide)
       {
		graph->SetPoint(pt, x, y*scaler);
		graph->SetPointError(pt, x_down, x_up, y_down*scaler, y_up*scaler);
       }
       else
       {
           
            graph->SetPoint(pt, x, y/scaler);
            graph->SetPointError(pt, x_down, x_up, y_down/scaler, y_up/scaler);
       }

	}
	return;
}

//================================================================================
//scale every point of the graph and corresponding y errors by the constants from the "scaler" array (one graph point = one constant)
//================================================================================
void scale_graph(Graph* graph, double* scaler, bool divide=0)
{
	for(int pt=0; pt<graph->GetN(); pt++)
	{
		double x,x_down,x_up,y,y_up,y_down;

		graph->GetPoint(pt,x,y);
	   x_down=graph->GetErrorXlow(pt);
	   x_up=graph->GetErrorXhigh(pt);
	   y_down=graph->GetErrorYlow(pt);
	   y_up=graph->GetErrorYhigh(pt);

       if(!divide)
       {
		graph->SetPoint(pt, x, y*scaler[pt]);
		graph->SetPointError(pt, x_down, x_up, y_down*scaler[pt], y_up*scaler[pt]);
       }
       else
       {
            graph->SetPoint(pt, x, y/scaler[pt]);
            graph->SetPointError(pt, x_down, x_up, y_down/scaler[pt], y_up/scaler[pt]);
       }

	}
	return;
}

//================================================================================
//replace y errors of the graph "grold" by the errors from "grtemplate" and rescale them according to the ratio of y values yold/ytemp
//================================================================================
void copy_graph_errors(Graph* grold, Graph* grtemplate)
{
	for(int pt=0; pt<grtemplate->GetN(); pt++)
	{
		double xtemp,ytemp, ytemp_up, ytemp_down;
		double xold,yold,xold_up,xold_down;

		grtemplate->GetPoint(pt,xtemp,ytemp);
		grold->GetPoint(pt,xold,yold);
		double scaler=(TMath::Abs(ytemp>0)) ? yold/ytemp : 0;
		
	   xold_down=grold->GetErrorXlow(pt);
		xold_up=grold->GetErrorXhigh(pt);
	   ytemp_down=grtemplate->GetErrorYlow(pt);
	   ytemp_up=grtemplate->GetErrorYhigh(pt);
		
		//cout<<pt<<" scaler:"<<scaler<<" yold:"<<ymove_down<<" ynew:"<<ymove_down*scaler<<endl;
		
		grold->SetPointError(pt, xold_down, xold_up, ytemp_down*scaler, ytemp_up*scaler);

	}
	return;
}
//================================================================================
//integral of a graph
//================================================================================
double graph_integral(Graph* graph, ErrorType err=none, int skipNbins=0)
{
	int ntmp=graph->GetN();
	const int nbins=ntmp;
	double integral=0;
	
	for(int bn=skipNbins; bn<nbins; bn++)
	{
		double x,y,xl,xr,yl,yr;
		graph->GetPoint(bn,x,y);
		xl=graph->GetErrorXlow(bn);
		xr=graph->GetErrorXhigh(bn);
		yl=graph->GetErrorYlow(bn);
		yr=graph->GetErrorYhigh(bn);
		
		double width=xl+xr;
		
		if(err==up) y=y+yr;
		else if(err==down) y=y-yl;
		
		integral+=width*y;
	}
	return integral;
}


//================================================================================
//ratio of two graphs with different bin centers (e.g. due to bin shift correction)
//================================================================================
Graph* ratio_of_graphs(Graph* gr1,Graph* gr2)
{
	int ntmp=gr1->GetN();
	const int npoints=ntmp;
	double x[npoints];
	double y[npoints];
	double xl[npoints];
	double xr[npoints];
	double yl[npoints];
	double yr[npoints];
	
	for(int pt=0; pt<npoints;pt++)
	{
		double x1,y1,x2,y2,xl1,xr1;
		gr1->GetPoint(pt,x1,y1);
		xl1=gr1->GetErrorXlow(pt);
		xr1=gr1->GetErrorXhigh(pt);
		
		//find interpolated value of y2 corresponding to x1
		y2=interpolate_y(gr2,x1);
		double ratio=(y2>0) ? (y1/y2) : 0;
		
		x[pt]=x1;
		y[pt]=ratio;
		xl[pt]=xl1;
		xr[pt]=xr1;
		
		//do not calculate error bars, that would be tricky (we will use error bars calculated from ratios without bin shift correction)
		yl[pt]=0;
		yr[pt]=0;
		
	}//bin loop
	
	Graph* graph=new TGraphAsymmErrors(npoints, x, y, xl, xr, yl, yr);
	return graph;
	
}
//===============================================================
//combine errors from 2 graphs
//===============================================================
TGraphAsymmErrors* CombineGraphErrors(TGraphAsymmErrors *gr1,TGraphAsymmErrors *gr2,bool setErrX=1)
{
	int npoints=gr1->GetN();
	const int npts=50;
	double gx[npts];
	double gx_l[npts];
	double gx_h[npts];
	double gy[npts];
	double gy_l[npts];
	double gy_h[npts];
	
	for(int i=0;i<npoints;++i){
		double x1,y1;
		gr1->GetPoint(i,x1,y1);
		double x1_down;
		double x1_up;
		
		if(setErrX==0) //set x errors to 0
		{
			x1_down=0; 
			x1_up=0; 
		}
		else
		{
			x1_down=gr1->GetErrorXlow(i);
			x1_up=gr1->GetErrorXhigh(i);
		}
		double y1_down=gr1->GetErrorYlow(i);
		double y1_up=gr1->GetErrorYhigh(i);

		double y2_down=gr2->GetErrorYlow(i);
		double y2_up=gr2->GetErrorYhigh(i);
	
		gx[i]=x1;
		gx_l[i]=x1_down;
		gx_h[i]=x1_up;
		gy[i]=y1;
		gy_l[i]=TMath::Sqrt(y1_down*y1_down+y2_down*y2_down);
		gy_h[i]=TMath::Sqrt(y1_up*y1_up+y2_up*y2_up);
		
	}
	TGraphAsymmErrors* graph=new TGraphAsymmErrors(npts, gx, gy, gx_l, gx_h, gy_l, gy_h);
	return graph;
		
}


//================================================================================
//functions for quick settings of histogram/tgraph (graphical) properties
//================================================================================
void tgraph_set_atributes(TGraphAsymmErrors* tgraph, float xmin, float xmax, float ymin, float ymax, Color_t line_color=kBlue, Color_t marker_color=kBlue, 
								  Color_t fill_color=kBlue, short line_width=1, short marker_size=1, short marker_style=24, int fill_style=0, int line_style=1)
{
	tgraph->GetXaxis()->SetLimits(xmin,xmax);
	tgraph->GetHistogram()->SetMaximum(ymax);
	tgraph->GetHistogram()->SetMinimum(ymin);
	tgraph->SetLineColor(line_color);
	tgraph->SetMarkerColor(marker_color);
	tgraph->SetFillColor(fill_color);
	tgraph->SetLineWidth(line_width);
	tgraph->SetMarkerSize(marker_size);
	tgraph->SetMarkerStyle(marker_style);
	tgraph->SetFillStyle(fill_style);
	tgraph->SetLineStyle(line_style);
	
	return;
	
}
void tgraph_set_atributes(TGraphErrors* tgraph, float xmin, float xmax, float ymin, float ymax, Color_t line_color=kBlue, Color_t marker_color=kBlue, 
								  Color_t fill_color=kBlue, short line_width=1, short marker_size=1, short marker_style=24, int fill_style=0, int line_style=1)
{
	tgraph->GetXaxis()->SetLimits(xmin,xmax);
	tgraph->GetHistogram()->SetMaximum(ymax);
	tgraph->GetHistogram()->SetMinimum(ymin);
	tgraph->SetLineColor(line_color);
	tgraph->SetMarkerColor(marker_color);
	tgraph->SetFillColor(fill_color);
	tgraph->SetLineWidth(line_width);
	tgraph->SetMarkerSize(marker_size);
	tgraph->SetMarkerStyle(marker_style);
	tgraph->SetFillStyle(fill_style);
	tgraph->SetLineStyle(line_style);

	return;
	
}

void histo_set_atributes(TH1* histo, TString title, TString x_title, TString y_title, float xmin, float xmax, float ymin, float ymax, Color_t line_color=kBlue, Color_t marker_color=kBlue, short line_width=1, short marker_size=1, short marker_style=24)
{
	histo->SetTitle(title);
	histo->GetXaxis()->SetRangeUser(xmin, xmax);
	histo->GetYaxis()->SetRangeUser(ymin, ymax);
	histo->GetXaxis()->SetTitle(x_title);
	histo->GetYaxis()->SetTitle(y_title);
	histo->SetLineColor(line_color);
	histo->SetMarkerColor(marker_color);
	histo->SetLineWidth(line_width);
	histo->SetMarkerSize(marker_size);
	histo->SetMarkerStyle(marker_style);
	
	return;
}

void histo_set_atributes(TH1* histo, Color_t line_color=kBlue, Color_t marker_color=kBlue, short line_width=1, short marker_size=1, short marker_style=24)
{
	histo->SetLineColor(line_color);
	histo->SetMarkerColor(marker_color);
	histo->SetLineWidth(line_width);
	histo->SetMarkerSize(marker_size);
	histo->SetMarkerStyle(marker_style);
	
	return;
}


//================================================================================
//Check if given unfolding qualities (backfolded measured ratio, curvature, ratio of successive iterations) pass the QA cuts
//================================================================================

bool qa_check(float test_ch,float test_bf,float test_curv, C_systems system, C_unfoldings unf, short bining, C_tests test, short r, short pTl)
{
	float epsilon=1E-6;
	float test_bf_cut=backfold_cut[test][system][bining][unf][r][pTl];//drop solutions with chi2 of ratio of backfoled/measured > this cut
	float test_ch_cut=change_cut[test][system][bining][unf][r][pTl];//drop solutions with chi2 of ratio of successive iterations > this cut (Bayes only)
	float test_curv_cut=curv_cut[test][system][bining][unf][r][pTl]; //drop solutions with large oscialations
	
	if(test_bf_cut<0+epsilon) test_bf_cut= backfold_cut_default[unf];
	if(test_ch_cut<0+epsilon) test_ch_cut= change_cut_default[unf];
	if(test_curv_cut<0+epsilon) test_curv_cut= curv_cut_default[unf];
	
	bool show=1;
	//cout<<"bf test:"<<test_bf<<" cut:"<<test_bf_cut<<endl;
	//cout<<"change test:"<<test_ch<<" cut:"<<test_ch_cut<<endl;
	//cout<<"curvature test:"<<test_curv<<" cut:"<<test_curv_cut<<endl;
	if (test_bf>test_bf_cut || test_ch>test_ch_cut|| test_curv>test_curv_cut) show=0;
	return show;
}


