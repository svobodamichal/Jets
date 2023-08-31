//Tests for comparing two histograms (with the same binning)

double chi2_test(const TH1D* count1, const TH1D* count2, float stop=50, float start=-50)
{
   // calculate the chi^2. count1 and count2 are the counts
   float chi2= 0.0;
	float nbins=0;
   int n= count1->GetNbinsX();
	//cout<<"chi2 test"<<endl;
   for (int i = 1 ; i <= n ; i++) {
      float cent=count1->GetBinCenter(i);
      if(cent<start || cent>stop)continue;
      float psum  = (count1->GetBinContent(i) + count2->GetBinContent(i));
      float pdiff = (count1->GetBinContent(i) - count2->GetBinContent(i));
		//cout<<"bin:"<<i<<" diff: "<<pdiff<<endl;
		//cout<<" h1: "<<count1->GetBinContent(i)<<" h2: "<<count2->GetBinContent(i)<<endl;
      if (psum > 1.0)
         chi2 = chi2 + (pdiff*pdiff)/psum;
      else
         chi2 = chi2 + (pdiff*pdiff);
		if(psum>0) nbins++;
   }
   return chi2/nbins;
}

double Canberra_test(const TH1D* count1, const TH1D* count2, float stop=50, float start=-50)
{
   float sum=0;
	float nbins=0;
   int n= count1->GetNbinsX();
   for (Int_t i = 1 ; i <= n ; i++) {
      float cent=count1->GetBinCenter(i);
		if(cent<start || cent>stop)continue;		
		//cout<<"Camberra cent:"<<cent<<endl;
      float adiff = TMath::Abs((count1->GetBinContent(i) - count2->GetBinContent(i)));
		//cout<<"adiff:"<<adiff<<endl;
      float rdiff =0;
      float min=TMath::Min(count1->GetBinContent(i),count2->GetBinContent(i));
		//cout<<"min:"<<min<<endl;
      if (min>0) 
		{
			rdiff =adiff/min;
			nbins++;
		}
		//cout<<"rdiff:"<<rdiff<<endl;
      sum+=rdiff;
		//cout<<"sum:"<<sum<<endl;
   }
return sum/nbins;
}


double KS_test(const TH1D* count1, const TH1D* count2, float stop=50, float start=-50)
{
   int nbins= count1->GetNbinsX();
   float n1=count1->Integral();
   float n2=count2->Integral();

   float sum1=0;
   float sum2=0;
   float max=0;
   for (Int_t i = 1 ; i <= nbins ; i++)
   {
      float cent=count1->GetBinCenter(i);
		float y1=count1->GetBinContent(i);
		float y2=count2->GetBinContent(i);

      if(cent<start || cent>stop)continue;
		//cout<<"KS cent:"<<cent<<endl;
      sum1+=(float) y1/n1;
      sum2+=(float) y2/n2;
      float diff=TMath::Abs(sum1-sum2);
      max=TMath::Max(max,diff);
   }
   return max;
}

double Curv_test(const TH1D* count1, float stop=50, float start=-50)
{

   // calculate the curvature. count1 are the counts
   float curv= 0.0;
	float nbins=0;
   int n= count1->GetNbinsX();
   for (int i = 2 ; i < n ; i++) {
		float centm1=count1->GetBinCenter(i-1);
		float centp1=count1->GetBinCenter(i+1);
		if(centm1<start || centp1>stop)continue;
      float wim1=count1->GetBinContent(i-1)/count1->GetBinWidth(i-1);
		float wi=count1->GetBinContent(i)/count1->GetBinWidth(i);
		float wip1=count1->GetBinContent(i+1)/count1->GetBinWidth(i+1);
      if(!wi>0)continue;
		float diff=((wip1-wi)-(wi-wim1))/wi;
      float difsqr=TMath::Power(diff,2);
      curv+=difsqr;
		if(difsqr>0)nbins++;
   }
   return curv/nbins;
}
