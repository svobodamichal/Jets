#include <iostream>
#include <sstream>

#include <sstream>
#include <iomanip>
#include <cmath>
#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh" 
#include "fastjet/internal/base.hh"

using namespace std;  //robotmon
using namespace fastjet;

void fastjet(){
vector<PseudoJet> jetTracks;
		float eta = trk->gMom().pseudoRapidity();
		float pT = trk->gPt();
    float eta = trk->gMom().pseudoRapidity();
    float phi = trk->gMom().phi();
    float dca = (mPrimVtx - trk->origin()).mag();
		float charged = trk->charge();

		  PseudoJet inputParticle(trk->gMom().x(), trk->gMom().y(), trk->gMom().z(), trk->gMom().mag());
        jetTracks.push_back(inputParticle);
		return 0;
		}

