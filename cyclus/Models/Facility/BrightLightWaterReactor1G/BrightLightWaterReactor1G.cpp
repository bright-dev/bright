// BrightLightWaterReactor1G.cpp
// Implements the BrightLightWaterReactor1G class

#include "BrightLightWaterReactor1G.h"

#include "Logger.h"
#include "CycException.h"

using namespace std;

// stateic class members
BrightLightWaterReactor1G::refcount = 0;
BrightLightWaterReactor1G::engine = NULL;

/* --------------------
 * all MODEL classes have these members
 * --------------------
 */

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
BrightLightWaterReactor1G::BrightLightWaterReactor1G(){
  // get new engine
  if (engine != NULL){
    // start bright if it hasn't been started
    if (0 ==  bright::BRIGHT_DATA.length()) {
      bright::bright_start();
      bright::write_text = 0;
      bright::write_hdf5 = 0;
    }
    std::string libname = BRIGHT_DATA + "/LWR.h5";
    bright::load_track_nucs_hdf5(libname, std::string(""));
    engine = new bright::LightWaterReactor1G(libname);
  }
  refcount++;
};

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
BrightLightWaterReactor1G::~BrightLightWaterReactor1G(){
  refcount--;
  if (refcount == 0){
    delete engine;
    engine = NULL;
  };
};

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
void BrightLightWaterReactor1G::initModuleMembers(QueryEngine* qe) {
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
void BrightLightWaterReactor1G::cloneModuleMembersFrom(FacilityModel* src) {
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
std::string BrightLightWaterReactor1G::str() {
  return FacilityModel::str();
};

/* ------------------- */ 


/* --------------------
 * all COMMUNICATOR classes have these members
 * --------------------
 */

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void BrightLightWaterReactor1G::receiveMessage(msg_ptr msg) {}

/* ------------------- */ 


/* --------------------
 * all FACILITYMODEL classes have these members
 * --------------------
 */

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
vector<rsrc_ptr> BrightLightWaterReactor1G::removeResource(Transaction order){
  vector<rsrc_ptr> vrp = vector<rsrc_ptr>();
  return vrp;
}
    
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void BrightLightWaterReactor1G::addResource(Transaction trans, std::vector<rsrc_ptr> manifest){}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void BrightLightWaterReactor1G::handleTick(int time){}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void BrightLightWaterReactor1G::handleTock(int time){}

/* ------------------- */ 


/* --------------------
 * all MODEL classes have these members
 * --------------------
 */

extern "C" Model* constructBrightLightWaterReactor1G() {
  return new BrightLightWaterReactor1G();
}

extern "C" void destructBrightLightWaterReactor1G(Model* model) {
  delete model;
}

/* ------------------- */ 

