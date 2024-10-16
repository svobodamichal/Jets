static const Double_t pTmin          = 0.2;  //min track transverse momentum cut
static const Double_t pTmax         = 30.0;  // max track pT
static const Double_t DcaCutMax      = 1.0;  //max single track DCA to PV
static const Double_t DcaCutMin      = 0.0;  //min single track DCA to PV
static const Int_t    RatioHit       = 0.52;   // TPC hits fitted/possible hits
static const Int_t    TpcCut         = 15;   // TPC hits fitted
static const Int_t    TpcMcCut       = 10;   // TPC hits fitted
static const Int_t    FlagMax        = 1000;   // track flag
static const Double_t EtaMax         = 1.0;  // max track pseudorapidity
static const Double_t zcut           = 30;    // zvertex cut
static const Double_t Chi2Cut        =100; //chi2 of the helix fit
static const Int_t refMultCutMin[3]={10,396,1}; //multiplicity cut for peripheral/central collisions
static const Int_t refMultCutMax[3]={31,2000, 2000};//multiplicity cut for peripheral/central collisions
//static const Double_t ProbKFCut      = 0.01; // probability of fit
//static const Double_t TrackLengthCut =40;   // min value for dEdxTrackLength

