enum binNorm { COUNTS, CPU };  //counts, counts per unit bin width

TH1D* rebin_histogram(TH1D* hinput, TH1D* htemplate, TString outname, binNorm input_bin_norm=COUNTS , binNorm output_bin_norm=COUNTS);
TH2D* rebin_histogram2D(TH2D* hinput, TH1D* htemplate_x, TH1D* htemplate_y,TString outname);
void convert_histo_from_Counts_to_CPU(TH1D* h1);
void convert_histo_from_CPU_to_Counts(TH1D* h1);

//**********************************
//useful functions
//**********************************

//function for rebinning 1D histogram hinput to the bin structure of the htemplate
TH1D* rebin_histogram(TH1D* hinput, TH1D* htemplate, TString outname, binNorm input_bin_norm,binNorm output_bin_norm)
{
	TH1D* houtput=(TH1D*) htemplate->Clone(outname);
	houtput->Reset();
	int nbins_old=hinput->GetNbinsX();
	int nbins_new=houtput->GetNbinsX();
	
	for(int i=1; i<=nbins_old; i++)
	{
		double low_edge_old=hinput->GetBinLowEdge(i);
		double width_old=hinput->GetBinWidth(i);
		double upper_edge_old=low_edge_old+width_old;
		double bin_content_old=hinput->GetBinContent(i);
		double bin_error_old=hinput->GetBinError(i);
		if(input_bin_norm==CPU) //make sure we have counts, not counts/bin_width
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
	
	if(output_bin_norm==CPU) //if desired normalize the histogram by the bin width
		convert_histo_from_Counts_to_CPU(houtput);
	return houtput;
}

//divide the content of the histogram bins and errors with the bin widths (this is same as h1->Scale(1,"width");)
void convert_histo_from_Counts_to_CPU(TH1D* h1)
{
	h1->Sumw2();
	h1->Scale(1,"width");
	/*
	int nbins=h1->GetNbinsX();
	for(int i=1; i<=nbins; i++)
	{
		double width=h1->GetBinWidth(i);
		double bin_content=h1->GetBinContent(i);
		double bin_error=h1->GetBinError(i);
		h1->SetBinContent(i,bin_content/width);
		h1->SetBinError(i,bin_error/width);
	}
	*/
}

//convert histograms with bin contents normalized by the bin width to histograms with bins containing raw counts
//multiply the content of the histogram bins and errors with the bin widths - this will convert bin contents of a type dN/dE to N
void convert_histo_from_CPU_to_Counts(TH1D* h1)
{
	int nbins=h1->GetNbinsX();
	for(int i=1; i<=nbins; i++)
	{
		double width=h1->GetBinWidth(i);
		double bin_content=h1->GetBinContent(i);
		double bin_error=h1->GetBinError(i);
		h1->SetBinContent(i,bin_content*width);
		h1->SetBinError(i,bin_error*width);
	}

}

TH2D* rebin_histogram2D(TH2D* hinput, TH1D* htemplate_x, TH1D* htemplate_y,TString outname)
{
	const int nbinsXnew=htemplate_x->GetNbinsX();
	const int nbinsXold=hinput->GetNbinsX();
	//const int nbinsYnew=htemplate_y->GetNbinsY();
	const int nbinsYnew=htemplate_y->GetNbinsX();
	const int nbinsYold=hinput->GetNbinsY();
	
	const double *binsXnew=htemplate_x->GetXaxis()->GetXbins()->GetArray();
	const double *binsXold=hinput->GetXaxis()->GetXbins()->GetArray(); //doesn't work, not used. If needed see binsYold
	/*double binsXold[nbinsXold];
	for (int binx = 1; binx <=nbinsXold;binx++){
		cout << hinput->GetXaxis()->GetBinLowEdge(binx) << endl;
		binsXold[binx-1]=hinput->GetXaxis()->GetBinLowEdge(binx);
	}
	*/
	const double *binsYnew=htemplate_y->GetXaxis()->GetXbins()->GetArray();
	//const double *binsYold=hinput->GetYaxis()->GetXbins()->GetArray(); //this method doesn't work for bins with automatically generated bin edges 
	//static const int NBINSYOLD = nbinsYold;
	//double binsYold_tmp[NBINSYOLD]; //works with 801 hard number
	int size = nbinsYold; 
	double *binsYold_tmp;
	binsYold_tmp = new double[size+1];
	for (int biny = 0; biny <=nbinsYold;biny++){
		//cout << "here" << endl;
		binsYold_tmp[biny] = hinput->GetYaxis()->GetBinLowEdge(biny + 1);
		//cout << binsYold_tmp[biny] << endl;
	}
	double *binsYold = binsYold_tmp;
	
	TH2D* hnewXoldY = new TH2D("hnewXoldY","new x-binning, old y-binning",nbinsXnew,binsXnew,nbinsYold,binsYold);
	for(int biny=1; biny<=nbinsYold;biny++)
	{
		TH1D* hprojectionX=(TH1D*) hinput->ProjectionX(Form("x_%i",biny),biny,biny);
		TH1D* hrebinX=(TH1D*) rebin_histogram(hprojectionX, htemplate_x, Form("hrebinX_%i",biny));
		for(int binx=1;binx<=nbinsXnew;binx++)
		{
			double cont=hrebinX->GetBinContent(binx);
			double err=hrebinX->GetBinError(binx);
			hnewXoldY->SetBinContent(binx,biny,cont);
			hnewXoldY->SetBinError(binx,biny,err);
		}
		delete hprojectionX;
		delete hrebinX;
	}
	delete binsYold_tmp;
	
	TH2D* hnewXnewY = new TH2D("hnewXnewY","new x-binning, new y-binning",nbinsXnew,binsXnew,nbinsYnew,binsYnew);
	//for(int binx=1; binx<=nbinsYold;binx++)
	for(int binx=1; binx<=nbinsXnew;binx++)	
	{
		TH1D* hprojectionY=(TH1D*) hnewXoldY->ProjectionY(Form("y_%i",binx),binx,binx);
		TH1D* hrebinY=(TH1D*) rebin_histogram(hprojectionY, htemplate_y, Form("hrebinY_%i",binx));
		for(int biny=1;biny<=nbinsYnew;biny++)
		{
			double cont=hrebinY->GetBinContent(biny);
			double err=hrebinY->GetBinError(biny);
			hnewXnewY->SetBinContent(binx,biny,cont);
			hnewXnewY->SetBinError(binx,biny,err);
            if (err > 0 && cont/err < sqrt(10)) {hnewXnewY->SetBinContent(binx, biny,0);hnewXnewY->SetBinError(binx,biny,0);};

        }
		delete hprojectionY;
		delete hrebinY;
	}
	
	//return hnewXoldY;
	return hnewXnewY;
}

TH2D* rebin_histogram2D_1Axis(TH2D* hinput, TH1D* htemplate, TString axis="X" ,TString outname="hrebinned",binNorm input_bin_norm=COUNTS /*what is the bin content of the input histogram? COUNTS (raw counts) or CPU (counts per axis unit)*/, 
binNorm output_bin_norm=COUNTS /*desired structure of the output histogram*/)
{
	const int nbinsnew=htemplate->GetNbinsX();
	const int nbinsXold=hinput->GetNbinsX();
	const int nbinsYold=hinput->GetNbinsY();
	
	const double *binsnew=htemplate->GetXaxis()->GetXbins()->GetArray();
	const double *binsXold=hinput->GetXaxis()->GetXbins()->GetArray();
	const double *binsYold=hinput->GetYaxis()->GetXbins()->GetArray();
		
	TH2D* hrebinned;
	if(axis=="X") hrebinned= new TH2D(outname,outname,nbinsnew,binsnew,nbinsYold,binsYold);
	else hrebinned= new TH2D(outname,outname,nbinsXold,binsXold,nbinsnew,binsnew);
	
	if(axis=="X")
	{
	for(int biny=1; biny<=nbinsYold;biny++)
	{
		TH1D* hprojectionX=(TH1D*) hinput->ProjectionX(Form("x_%i",biny),biny,biny);
		TH1D* hrebinX=(TH1D*) rebin_histogram(hprojectionX, htemplate, Form("hrebinX_%i",biny),input_bin_norm,output_bin_norm);
		for(int binx=1;binx<=nbinsnew;binx++)
		{
			double cont=hrebinX->GetBinContent(binx);
			double err=hrebinX->GetBinError(binx);
			hrebinned->SetBinContent(binx,biny,cont);
			hrebinned->SetBinError(binx,biny,err);
		}
		delete hprojectionX;
		delete hrebinX;
	}
		
	}
	else
	{
		for(int binx=1; binx<=nbinsXold;binx++)
	{
		TH1D* hprojectionY=(TH1D*) hinput->ProjectionY(Form("x_%i",binx),binx,binx);
		TH1D* hrebinY=(TH1D*) rebin_histogram(hprojectionY, htemplate, Form("hrebinY_%i",binx),input_bin_norm,output_bin_norm);
		for(int biny=1;biny<=nbinsnew;biny++)
		{
			double cont=hrebinY->GetBinContent(biny);
			double err=hrebinY->GetBinError(biny);
			hrebinned->SetBinContent(binx,biny,cont);
			hrebinned->SetBinError(binx,biny,err);
		}
		delete hprojectionY;
		delete hrebinY;
	}
	}
	
	return hrebinned;
	
}
//EOF
