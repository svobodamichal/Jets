#include "TF1.h"
#include "TFile.h"
#include "TRandom.h"
#include "TPythia6.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "ThrmRecPythia.h"
#include "ThrmFourVector.h"
#include "Riostream.h"


void fragment(UInt_t PID, Double_t pT, TClonesArray *pararr)
{

// SETTING RANDOM SEED
TDatime dt;
UInt_t curtime = dt.Get(); //current time
UInt_t procid = gSystem->GetPid(); //process id
UInt_t seed = curtime - procid; //uniqe seed
gRandom->SetSeed(seed);
		
pythia = new TPythia6();
pythia->SetMRPY(1, seed); //set random seed
 
Double_t phi_parton = 2*TMath::Pi()*gRandom->Uniform(-1., +1.);
Double_t eta_parton = gRandom->Uniform(-1., +1.);
Double_t theta = 2.0*TMath::ATan(TMath::Exp(-1*eta_parton));
Double_t E = pT/TMath::Sin(theta);
      
pythia->Py1ent(1, PID, E, theta, phi_parton); 
pythia->Pyexec();
Int_t fin = pythia->ImportParticles(simarr, "Final");
return;
}

void analyse()
{
//initial setup
Double_t parton_pT=1; //desired parton pT in GeV
Double_t parton_PID=1; //PID of fragmenting parton: 1=d, 2=u, 3=s, 4=c, 21=g	

TClonesArray *simarr = new TClonesArray("TParticle", 1000); //array for saving the created particles

//run PYTHIA fragmentation
fragment(parton_PID,parton_pT,simarr)

//loop over created particles
Int_t nparticles = simarr->GetEntries();
for(Int_t ipart = 0; ipart < nparticles; ipart++)
	{
	  TParticle *particle = (TParticle*)simarr->At(ipart);
	  /*
	  //if you need the TLorentzVector
	  TLorentzVector partlv; //jet constituent lorentz vector
	  particle->Momentum(partlv); //create TLorentzVector from TParticle
	  partlv.SetPtEtaPhiM(partlv.Pt(), partlv.Eta(), partlv.Phi(), 0); //set M=0, probably not necessary
	  */
	  
	  Double_t pT=particle->Pt();
	  UInt_t pdg=particle->GetPdgCode();
	  cout<<"particle PDG:"<<pdg<<" pT: "<<pT<<endl;
	}//end of loop

	return;
}