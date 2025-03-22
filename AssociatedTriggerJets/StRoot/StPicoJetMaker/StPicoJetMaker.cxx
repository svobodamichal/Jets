#include "StPicoJetMaker.h"
#include <map>
#include <string>

ClassImp(StPicoJetMaker)

// _________________________________________________________
/*
StPicoJetMaker::StPicoJetMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName) :
        StMaker(name), mPicoDst(NULL), mPicoCuts(NULL), mBField(0.), mOutList(NULL), mMakerMode(StPicoJetMaker::kAnalyze),
        mMcMode(false), mOutputTreeName("picoJetTree"), mOutputFileBaseName(outputBaseFileName), mPicoDstMaker(picoMaker),
        mPicoEvent(NULL), mTree(NULL), mOutputFileTree(NULL), mOutputFileList(


	// _________________________________________________________                                                                                     */

// _________________________________________________________
StPicoJetMaker::StPicoJetMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName,  char const* inputHFListHFtree = "") : StMaker(name), mPicoDst(NULL), mBField(0.), mOutList(NULL), mMakerMode(StPicoJetMaker::kAnalyze), mMcMode(false),
  mOutputTreeName("picoHFtree"), mOutputFileBaseName(outputBaseFileName), mInputFileName(inputHFListHFtree), mPicoDstMaker(picoMaker), mPicoEvent(NULL), mTree(NULL), mOutputFileTree(NULL), mOutputFileList(NULL) {
  // -- constructor
}
                                                                                                                    

/*
StPicoJetMaker::StPicoJetMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName) : StPicoJetMaker(name)
{
  //   mInputFileName = name;
}
*/

// _________________________________________________________
StPicoJetMaker::~StPicoJetMaker() {
   // -- destructor 
  
  if (mPicoCuts)
    delete mPicoCuts;
  mPicoCuts = NULL;

}

// _________________________________________________________
Int_t StPicoJetMaker::Init() {
  // -- Inhertited from StMaker 
  //    NOT TO BE OVERWRITTEN by daughter class
  //    daughter class should implement InitHF()

  // -- check for cut class
  if (!mPicoCuts)
    mPicoCuts = new StPicoCuts();
  mPicoCuts->init();
  
  // -- file which holds list of histograms
  mOutputFileList = new TFile(Form("%s.%s.root", mOutputFileBaseName.Data(), GetName()), "RECREATE");
  mOutputFileList->SetCompressionLevel(1);

  if (mMakerMode == StPicoJetMaker::kWrite) {
    mOutputFileTree = new TFile(Form("%s.%s.root", mOutputFileBaseName.Data(), mOutputTreeName.Data()), "RECREATE");
    mOutputFileTree->SetCompressionLevel(1);
    mOutputFileTree->cd();

    // -- create OutputTree
    int BufSize = (int)pow(2., 16.);
    int Split = 1;
    if (!mTree) 
      mTree = new TTree("T", "T", BufSize);
    mTree->SetAutoSave(1000000); // autosave every 1 MBytes
  }

  // -- disable automatic adding of objects to file
  bool oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(false);

  // -- add list which holds all histograms  
  mOutList = new TList();
  mOutList->SetName(GetName());
  mOutList->SetOwner(true);

  // -- create event stat histograms
  initializeEventStats();

  // -- call method of daughter class
  InitJets();

  TH1::AddDirectory(oldStatus);

  // -- reset event to be in a defined state
  resetEvent();

  precomputeWeights();
  return kStOK;
}

// _________________________________________________________
Int_t StPicoJetMaker::Finish() {
  // -- Inherited from StMaker
  //    NOT TO BE OVERWRITTEN by daughter class
  //    daughter class should implement FinishHF()

  if (mMakerMode == StPicoJetMaker::kWrite) {
    mOutputFileTree->cd();
    mOutputFileTree->Write();
    mOutputFileTree->Close();
  }

  mOutputFileList->cd();
  mOutList->Write(mOutList->GetName(), TObject::kSingleKey);
	mOutList->Delete();
  
  // -- call method of daughter class
  FinishJets();

  mOutputFileList->Close();

  return kStOK;
}

// _________________________________________________________
void StPicoJetMaker::resetEvent() {
  // -- reset event
  mIdxPicoParticles.clear();
}

// _________________________________________________________
void StPicoJetMaker::Clear(Option_t *opt) {
  // -- Inherited from StMaker
  //    NOT TO BE OVERWRITTEN by daughter class
  //    daughter class should implement ClearHF()

  resetEvent();
}

// _________________________________________________________
Int_t StPicoJetMaker::Make() {
  // -- Inherited from StMaker
  //    NOT TO BE OVERWRITTEN by daughter class
  //    daughter class should implement MakeHF()
  // -- isPion, isKaon, isProton methods are to be 
  //    implemented by daughter class (
  //    -> methods of StHFCuts can and should be used

  if (!mPicoDstMaker) {
    LOG_WARN << " StPicoJetMaker - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  mPicoDst = mPicoDstMaker->picoDst();
  if (!mPicoDst) {
    LOG_WARN << " StPicoJetMaker - No PicoDst! Skip! " << endm;
    return kStWarn;
  }
  
  Int_t iReturn = kStOK;

  if (setupEvent()) {
    UInt_t nTracks = mPicoDst->numberOfTracks();

      int runNumber = mPicoDst->event()->runId();

      auto weightIt = mWeightMap.find(runNumber);
      if (weightIt != mWeightMap.end()) {
          double weightEVT = weightIt->second;

          // Fill histograms with the precomputed weight
          static_cast<TH1D*>(mOutList->FindObject("hrunId_weighted"))->Fill(runNumber, weightEVT);
      } else {
          std::cerr << "Warning: Precomputed weight not found for run number " << runNumber << std::endl;
      }

    // -- Fill vectors of particle types
    if (mMakerMode == StPicoJetMaker::kWrite || mMakerMode == StPicoJetMaker::kAnalyze) {
      for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack) {
        StPicoTrack* trk = mPicoDst->track(iTrack);

				if (trk->pMom().Perp() > 30) return kStOK; //throw out events with trk > 30 GeV/c
     // if (!trk || !mPicoCuts->isGoodTrack(trk)) continue;
		//good primary tracks only

			 if (!trk || !mPicoCuts->isGoodPrimaryTrack(trk)) continue;

        mIdxPicoParticles.push_back(iTrack);
      }
    }

    // -- call method of daughter class
    iReturn = MakeJets();
    //TODO: Fill good event histograms here - expand for other observables
    static_cast<TH1I*>(mOutList->FindObject("hevents_acc"))->Fill(1);
    static_cast<TH1I*>(mOutList->FindObject("hrefmult_acc"))->Fill(mPicoDst->event()->refMult());
    static_cast<TH1F*>(mOutList->FindObject("hzvertex_acc"))->Fill(mPrimVtx.z());
    static_cast<TH1D*>(mOutList->FindObject("hdeltaz_acc"))->Fill(mPrimVtx.z() - mPicoDst->event()->vzVpd());
    static_cast<TH2D*>(mOutList->FindObject("hz_refmult_acc"))->Fill(mPrimVtx.z(), mPicoDst->event()->refMult());
    static_cast<TH1I*>(mOutList->FindObject("hrunId_acc"))->Fill(mPicoDst->event()->runId());


  }
  
  // -- save information about all events, good or bad
  if (mMakerMode == StPicoJetMaker::kWrite)
    mTree->Fill();

  // -- fill basic event histograms - for all events
  //TODO: Fill all event histograms here - expand for other observables
  static_cast<TH1I*>(mOutList->FindObject("hevents"))->Fill(1);
  static_cast<TH1I*>(mOutList->FindObject("hrefmult"))->Fill(mPicoDst->event()->refMult());
  static_cast<TH1F*>(mOutList->FindObject("hzvertex"))->Fill(mPrimVtx.z());
  static_cast<TH1D*>(mOutList->FindObject("hdeltaz"))->Fill(mPrimVtx.z() - mPicoDst->event()->vzVpd());
  static_cast<TH2D*>(mOutList->FindObject("hz_refmult"))->Fill(mPrimVtx.z(), mPicoDst->event()->refMult());
  static_cast<TH1I*>(mOutList->FindObject("hrunId"))->Fill(mPicoDst->event()->runId());

  // -- reset event to be in a defined state
  resetEvent();
  return (kStOK && iReturn);
}

// _________________________________________________________
bool StPicoJetMaker::setupEvent() {

  // -- fill members from pico event, check for good events and fill event statistics

  mPicoEvent = mPicoDst->event();
  
  mBField = mPicoEvent->bField();
  mPrimVtx = mPicoEvent->primaryVertex();
  
  int aEventStat[mPicoCuts->eventStatMax()];
  
  bool bResult = mPicoCuts->isGoodEvent(mPicoDst, aEventStat);

  fillEventStats(aEventStat);

  return bResult;
}

// _________________________________________________________
void StPicoJetMaker::initializeEventStats() {
  // -- Initialize event statistics histograms

  //TODO: Add event histograms here - good and all events version

  int refmultbins = 130;
  float refmultmin = 0;
  float refmultmax = 650;

  int zbins = 50;
  float zmin = -50;
  float zmax = 50;


  const char *aEventCutNames[]   = {"all", "good run", "trigger", "#it{v}_{z}", "#it{v}_{z}-#it{v}^{VPD}_{z}", "refMult", "accepted"};

  mOutList->Add(new TH1F("hEventStat0","Event cut statistics 0;Event Cuts;Events", mPicoCuts->eventStatMax(), -0.5, mPicoCuts->eventStatMax()-0.5));
  TH1F *hEventStat0 = static_cast<TH1F*>(mOutList->Last());

  mOutList->Add(new TH1F("hEventStat1","Event cut statistics 1;Event Cuts;Events", mPicoCuts->eventStatMax(), -0.5, mPicoCuts->eventStatMax()-0.5));
  TH1F *hEventStat1 = static_cast<TH1F*>(mOutList->Last());

  for (unsigned int ii = 0; ii < mPicoCuts->eventStatMax(); ii++) {
    hEventStat0->GetXaxis()->SetBinLabel(ii+1, aEventCutNames[ii]);
    hEventStat1->GetXaxis()->SetBinLabel(ii+1, aEventCutNames[ii]);
  }

  //TODO: Add event ID histograms for AuAu Run 2016 and observables vs day information


  mOutList->Add(new TH1D("hrunId_weighted", "minimum bias events", 90913, 15076101, 15167014)); //15076101−15167014

  //All event histograms
  mOutList->Add(new TH1I("hevents", "number of events", 2, 0, 2));
  mOutList->Add(new TH1I("hrefmult", "Reference multiplicity", refmultbins, refmultmin, refmultmax));
  mOutList->Add(new TH1F("hzvertex", "z-position of primary vertex", zbins, zmin, zmax));
  mOutList->Add(new TH1D("hdeltaz", "zTPC-zVPD; #Delta [cm]", 80, -10, 10));
  mOutList->Add(new TH2D("hz_refmult", "zvertex vs refmult; z [cm]; refMult", zbins, zmin, zmax, refmultbins, refmultmin, refmultmax));

  mOutList->Add(new TH1I("hrunId", "runId", 90913, 15076101, 15167014));

  //Accepted event histograms
  mOutList->Add(new TH1I("hevents_acc", "number of events", 2, 0, 2));
  mOutList->Add(new TH1I("hrefmult_acc", "Reference multiplicity", refmultbins, refmultmin, refmultmax));
  mOutList->Add(new TH1F("hzvertex_acc", "z-position of primary vertex", zbins, zmin, zmax));
  mOutList->Add(new TH1D("hdeltaz_acc", "zTPC-zVPD; #Delta [cm]", 80, -10, 10));
  mOutList->Add(new TH2D("hz_refmult_acc", "zvertex vs refmult; z [cm]; refMult", zbins, zmin, zmax, refmultbins, refmultmin, refmultmax));
  
  mOutList->Add(new TH1I("hrunId_acc", "accepted events runId", 90913, 15076101, 15167014)); //15076101−15167014
}

//________________________________________________________________________
void StPicoJetMaker::fillEventStats(int *aEventStat) {
  // -- Fill event statistics 

  TH1F *hEventStat0 = static_cast<TH1F*>(mOutList->FindObject("hEventStat0"));
  TH1F *hEventStat1 = static_cast<TH1F*>(mOutList->FindObject("hEventStat1"));

  for (unsigned int idx = 0; idx < mPicoCuts->eventStatMax() ; ++idx) {
    if (!aEventStat[idx])
      hEventStat0->Fill(idx);
  }
  
  for (unsigned int idx = 0; idx < mPicoCuts->eventStatMax(); ++idx) {
    if (aEventStat[idx])
      break;
    hEventStat1->Fill(idx);
  }
}

//________________________________________________________________________
std::map<int, RunData> StPicoJetMaker::readDataFromFile(const std::string& filename) {
    std::map<int, RunData> dataMap;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return dataMap;
    }

    int runNumber, numberOfEvents;
    double sampledLuminosity, prescale, livetime;

    // Reading line-by-line and extracting data into the struct
    while (file >> runNumber >> numberOfEvents >> sampledLuminosity >> prescale >> livetime) {
        RunData runData = {runNumber, numberOfEvents, sampledLuminosity, prescale, livetime};
        dataMap[runNumber] = runData;
    }

    file.close();
    return dataMap;
}
//________________________________________________________________________
double StPicoJetMaker::calculateWeight(const RunData& htRunData, const RunData& mbRunData) {
    if (htRunData.sampledLuminosity == 0 || htRunData.prescale == 0 || mbRunData.sampledLuminosity == 0 || mbRunData.prescale == 0) {
            std::cerr << "Error: sampledLuminosity or prescale cannot be zero in weight calculation." << std::endl;
            return 0;}

    // Extract data for calculation
    double nEtvHT = 1.0;
    double psMB = mbRunData.prescale;
    double psHT = htRunData.prescale;
    double ltHT = htRunData.livetime;
    double ltMB = mbRunData.livetime;
    double nMBTotal = mbRunData.numberOfEvents;  // Assuming this is the total for MB
    double nHTTotal = htRunData.numberOfEvents;  // Assuming this is the total for HT
    double rMB = 0.88;   // Replace with actual rate if different
    double rHT = 0.98;    // Replace with actual rate if different

    cout << psMB << " " << psHT << " " << ltHT << " " << ltMB << " " << nMBTotal << " " << nHTTotal << " " << rMB << " " << rHT << endl;

    double ratioPS = psMB / psHT;
    double ratioLT = ltHT / ltMB;
    double nevtMB = nMBTotal * rMB;
    double nevtHT = nHTTotal * rHT;
    double ratioNevt = nevtMB / nevtHT;

    // Calculate the contribution for this run
    double weight = nEtvHT * ratioPS * ratioLT * ratioNevt;


    cout << "Weight: " << weight << endl;
    return weight;
}
//________________________________________________________________________
void StPicoJetMaker::precomputeWeights() {

    std::string file1 = "/gpfs01/star/pwg/svomich/Jets/AssociatedTriggerJets/StRoot/StPicoJetMaker/BHT2VPDMB-30_matched_cleaned.txt";
    std::string file2 = "/gpfs01/star/pwg/svomich/Jets/AssociatedTriggerJets/StRoot/StPicoJetMaker/VPDMB-30_matched_cleaned.txt";

    // Read data from both files into separate maps
    std::map<int, RunData> bhtRunDataMap = readDataFromFile(file1);
    std::map<int, RunData> vpdRunDataMap = readDataFromFile(file2);

    // Check if the maps have been successfully filled
    if (bhtRunDataMap.empty()) {
        std::cerr << "Warning: BHT file data is empty or could not be read." << std::endl;
    }
    if (vpdRunDataMap.empty()) {
        std::cerr << "Warning: VPD file data is empty or could not be read." << std::endl;
    }

    for (const auto& pair : bhtRunDataMap) {
        int runNumber = pair.first;
        const RunData& bhtRunData = pair.second;

        auto vpdIt = vpdRunDataMap.find(runNumber);
        if (vpdIt != vpdRunDataMap.end()) {
            const RunData& vpdRunData = vpdIt->second;
            mWeightMap[runNumber] = calculateWeight(bhtRunData, vpdRunData);
        } else {
            std::cerr << "Warning: Run number " << runNumber << " not found in VPD data map." << std::endl;
        }
    }
}
