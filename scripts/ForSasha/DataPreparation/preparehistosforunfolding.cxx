#include "util.h"

// Will take input histograms from .root file and add together histos for central and peripheral collisions. Output: .root file.

using namespace std;

void preparehistosforunfolding(string prod = "trig")
{
    // Open the input ROOT file and create the output ROOT file
    TFile *f1 = TFile::Open("output.root");
    TFile *fout1 = new TFile(Form("histos_HT2_%s.root", prod.c_str()), "recreate");

    int j = 0;
    int k = 0;
    int l = 0;

    TString cf[7][7];

    TList* list1 = (TList *)f1->Get("stPicoHFJetMaker");

    double R;
    array<int, 7> pTarr = {0, 3, 4, 5, 6, 7, 9};
    array<double, 3> Rarr = {0.2, 0.3, 0.4};
    array<int, 7> centbins = {1, 2, 3, 4, 5, 6, 7};
    array<string, 7> centrality = {"0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-80%"};

    TH1D* hfull[7][7];
    TH1D* hresfull[7][7];

    for(double &R : Rarr)
    {
        k = 0;
        for(int &centbin : centbins)
        {
            j = 0;
            for(int &pTlead : pTarr)
            {
                cf[j][k] = Form("hfpT_pTl%i_R0%.0f_centbin%i", pTlead, R*10, centbin);

                hfull[j][k] = (TH1D*)list1->FindObject(cf[j][k]);
                if (hfull[j][k]->GetEntries() == 0) {
                    cout << "empty histogram: " << hfull[j][k]->GetName() << ", skipping! " << endl;
                    continue;
                }

                hresfull[j][k] = (TH1D*)hfull[j][k]->Clone(hfull[j][k]->GetName());
                hresfull[j][k]->SetTitle("");
                hresfull[j][k]->SetXTitle("p_{T,jet}^{reco} [GeV/#it{c}]");
                hresfull[j][k]->Rebin(10);

                j++;
            }
            k++;
        }

        // No merging of centrality bins 0 and 1, and 7 and 8, since they are already merged

        for (int i = 0; i < 7; i++) {
            TString namecent = Form("hfpT_pTl%d_R0%.0f_peripheral", pTarr[i], R*10);
            hresfull[i][6]->Write(namecent);

            namecent = Form("hfpT_pTl%d_R0%.0f_central", pTarr[i], R*10);
            hresfull[i][0]->Write(namecent);
        }

        l++;
    }

    fout1->Close();
}
