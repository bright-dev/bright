// BrightLightWaterReactor1G.cpp
// Implements the BrightLightWaterReactor1G class
#include "BrightLightWaterReactor1G.h"

using namespace std;

// static class members
unsigned int BrightLightWaterReactor1G::refcount = 0;
bright::LightWaterReactor1G * BrightLightWaterReactor1G::engine = NULL;

BrightLightWaterReactor1G::BrightLightWaterReactor1G(){
  // get new engine
  if (engine != NULL){
    // start bright if it hasn't been started
    if (0 ==  bright::BRIGHT_DATA.length()) {
      bright::bright_start();
      bright::write_text = 0;
      bright::write_hdf5 = 0;
    }
    std::string libname = bright::BRIGHT_DATA + "/LWR.h5";
    bright::load_track_nucs_hdf5(libname, std::string(""));
    engine = new bright::LightWaterReactor1G(libname, bright::lwr_defaults);
    engine.target_BU = 50.0;
  }
  refcount++;
  month_in_cycle_ = 1;
  cycle_length_ = 3;
};

BrightLightWaterReactor1G::~BrightLightWaterReactor1G(){
  refcount--;
  if (refcount == 0){
    delete engine;
    engine = NULL;
  };
};


/**
  TICK
  if stocks are empty, ask for a batch
  offer anything in the inventory
  if we're at the end of a cycle
     - begin the cycle
       - move currCore batch to inventory
       - move stocks batch to currCore
       - reset month_in_cycle clock
 
  TOCK
  advance month_in_cycle
  send appropriate materials to fill ordersWaiting.
 
  RECIEVE MATERIAL
  put it in stocks
 
  SEND MATERIAL
  pull it from inventory
  
 */

void BrightLightWaterReactor::init(xmlNodePtr cur) { 
  FacilityModel::init(cur);
  
  // move XML pointer to current model
  cur = XMLinput->get_xpath_element(cur,"model/BrightLightWaterReactor");

  // initialize ordinary objects
  capacity_ = strtod(XMLinput->get_xpath_content(cur,"capacity"), NULL);
  //cycle_length_ = strtod(XMLinput->get_xpath_content(cur,"cycletime"), NULL);
  lifetime_ = strtol(XMLinput->get_xpath_content(cur,"lifetime"), NULL, 10);
  startConstrYr_ = strtol(XMLinput->get_xpath_content(cur,"startConstrYear"), NULL, 10);
  startConstrMo_ = strtol(XMLinput->get_xpath_content(cur,"startConstrMonth"), NULL, 10);
  startOpYr_ = strtol(XMLinput->get_xpath_content(cur,"startOperYear"), NULL, 10);
  startOpMo_ = strtol(XMLinput->get_xpath_content(cur,"startOperMonth"), NULL, 10);
  licExpYr_ = strtol(XMLinput->get_xpath_content(cur,"licExpYear"), NULL, 10);
  licExpMo_ = strtol(XMLinput->get_xpath_content(cur,"licExpMonth"), NULL, 10);
  state_ = XMLinput->get_xpath_content(cur,"state");
  typeReac_ = XMLinput->get_xpath_content(cur,"typeReac");
  CF_ = strtod(XMLinput->get_xpath_content(cur,"elecCF"), NULL);

  // all facilities require commodities - possibly many
  string recipe_name;
  xmlNodeSetPtr nodes = XMLinput->get_xpath_elements(cur, "fuelpair");

  // for each fuel pair, there is an in and an out commodity
  for (int i = 0; i < nodes->nodeNr; i++){
    xmlNodePtr pair_node = nodes->nodeTab[i];

    // get commods
    in_commod = XMLinput->get_xpath_content(pair_node,"incommodity");
    out_commod = XMLinput->get_xpath_content(pair_node,"outcommodity");
  };

  stocks_ = deque<Fuel>();
  currCore_ = deque<Fuel>();
  inventory_ = deque<Fuel>();
  ordersWaiting_ = deque<msg_ptr>();
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
void BrightLightWaterReactor::copy(BrightLightWaterReactor* src) {

  FacilityModel::copy(src);

  capacity_ = src->capacity_;
  cycle_length_ = src->cycle_length_;
  lifetime_ = src->lifetime_;
  month_in_cycle_ = src->month_in_cycle_;
  startConstrYr_ = src->startConstrYr_;
  startOpYr_ = src->startOpYr_;
  startOpMo_ = src->startOpMo_;
  licExpYr_ = src->licExpYr_;
  licExpMo_ = src->licExpMo_;
  state_ = src->state_;
  typeReac_ = src->typeReac_;
  CF_ = src->CF_;
  in_commod = src->in_commod;
  out_commod = src->out_commod;

  stocks_ = deque<Fuel>();
  currCore_ = deque< pair<string, mat_rsrc_ptr > >();
  inventory_ = deque<Fuel>();
  ordersWaiting_ = deque<msg_ptr>();
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
void BrightLightWaterReactor::copyFreshModel(Model* src) {
  copy(dynamic_cast<BrightLightWaterReactor*>(src));
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
std::string BrightLightWaterReactor::str() { 
  std::string s = FacilityModel::str(); 
  s += "    converts commodity '"
    + in_commod
    + "' into commodity '"
    + out_commod
    + "'.";
  return s;
};

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
void BrightLightWaterReactor::beginCycle() {
  if ( !stocks_.empty() ) {
    // move stocks batch to currCore
    string batchCommod = stocks_.front().first;
    mat_rsrc_ptr batchMat = stocks_.front().second;
    stocks_.pop_front();
    Fuel inBatch;
    inBatch = make_pair(batchCommod, batchMat);
    LOG(LEV_DEBUG2, "BLWR1G") << "Adding a new batch to the core";
    currCore_.push_back(inBatch);
    // reset month_in_cycle_ clock
    month_in_cycle_ = 1;
  } else {
    LOG(LEV_DEBUG3, "BLWR1G") << "Beginning a cycle with an empty core. Why??";
    // wait for a successful transaction to fill the stocks.
    // reset the cycle month to zero 
    month_in_cycle_=0;
  }
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
void BrightLightWaterReactor::endCycle() {
  LOG(LEV_DEBUG2, "BLWR1G") << "Ending a cycle.";
  month_in_cycle_ = 0;
  if (currCore_.size() == 0) {
    LOG(LEV_DEBUG3, "BLWR1G") << "Ended a cycle with an empty core. Why??";
    return;
  }

  // move a batch out of the core 
  string batchCommod = currCore_.front().first;
  mat_rsrc_ptr batchMatFeed = currCore_.front().second;
  currCore_.pop_front();

  // Get a pyne composition for bright to use
  map<int, double>::iterator i;
  map<int, double> incomp = map<int, double>();
  map<int, double> batch_feed = batchMatFeed.isoVector().comp().map();
  for (map<int, double>::iterator i = batch_feed.begin(); i != batch_feed.end(); i++)
    incomp[pyne::nucname::zzaaam(i->first)] = i->second;

  // do the calculation
  pyne::Material mat_prod = engine.calc(incomp);
  double mass_prod = mat_prod.mass;
  map<int, double> outcomp = map<int, double>();
  CompMapPtr outComp(new CompMap(MASS));
  for (i = mat_prod.comp.begin(); i != mat_prod.comp.end(); i++)
    (*outComp)[(i->first)/10] = (i->second) * mass_prod;

  // change the composition to the compositon of the spent fuel type
  mat_rsrc_ptr batchMatProd = mat_rsrc_ptr(new Material(outComp));

  // move converted material into Inventory
  Fuel outBatch;
  outBatch = make_pair(out_commod, batchMatProd);
  inventory_.push_back(outBatch);
};

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
void BrightLightWaterReactor::receiveMessage(msg_ptr msg) {
  // is this a message from on high? 
  if(msg->trans().supplier()==this){
    // file the order
    ordersWaiting_.push_front(msg);
    LOG(LEV_INFO5, "BLWR1G") << name() << " just received an order.";
  }
  else {
    throw CycException("BrightLightWaterReactor is not the supplier of this msg.");
  }
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
vector<rsrc_ptr> BrightLightWaterReactor::removeResource(Transaction order) {
  double newAmt = 0;

  mat_rsrc_ptr m;
  mat_rsrc_ptr newMat;
  mat_rsrc_ptr toAbsorb;

  // start with an empty manifest
  vector<rsrc_ptr> toSend;

  // pull materials off of the inventory stack until you get the trans amount
  while (order.resource()->quantity() > newAmt && !inventory_.empty() ) {
    for (deque<Fuel>::iterator iter = inventory_.begin(); 
        iter != inventory_.end(); 
        iter ++){
      // be sure to get the right commodity
      if (iter->first == order.commod()) {
        m = iter->second;

        // start with an empty material
        newMat = mat_rsrc_ptr(new Material());

        // if the inventory obj isn't larger than the remaining need, send it as is.
        if (m->quantity() <= (order.resource()->quantity() - newAmt)) {
          newAmt += m->quantity();
          newMat->absorb(m);
          inventory_.pop_front();
        } else { 
          // if the inventory obj is larger than the remaining need, split it.
          toAbsorb = m->extract(order.resource()->quantity() - newAmt);
          newAmt += toAbsorb->quantity();
          newMat->absorb(toAbsorb);
        }
        toSend.push_back(newMat);
        LOG(LEV_DEBUG2, "BLWR1G") <<"BrightLightWaterReactor "<< ID()
          <<"  is sending a mat with mass: "<< newMat->quantity();
      }
    }
  }    
  return toSend;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
void BrightLightWaterReactor::addResource(Transaction trans, std::vector<rsrc_ptr> manifest) {
  // grab each material object off of the manifest
  // and move it into the stocks.
  for (vector<rsrc_ptr>::iterator thisMat=manifest.begin();
       thisMat != manifest.end();
       thisMat++) {
    LOG(LEV_DEBUG2, "BLWR1G") <<"BrightLightWaterReactor " << ID() << " is receiving material with mass "
        << (*thisMat)->quantity();
    stocks_.push_front(make_pair(trans.commod(), boost::dynamic_pointer_cast<Material>(*thisMat)));
  }
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
void BrightLightWaterReactor::handleTick(int time) {
  LOG(LEV_INFO3, "BLWR1G") << name() << " is ticking {";

  // if at beginning of cycle, beginCycle()
  // if stocks are empty, ask for a batch
  // offer anything in the inventory
  
  // BEGIN CYCLE
  if(month_in_cycle_ == 1){
    LOG(LEV_INFO4, "BLWR1G") << " Beginning a new cycle";
    this->beginCycle();
  };

  makeRequests();
  makeOffers();
  LOG(LEV_INFO3, "BLWR1G") << "}";
}

void BrightLightWaterReactor::makeRequests(){
  // MAKE A REQUEST
  if(this->stocksMass() != 0) {
    return;
  }

  // It chooses the next in/out commodity pair in the preference lineup
  Recipe request_commod_pair;
  Recipe offer_commod_pair;
  request_commod_pair = fuelPairs_.front().first;
  offer_commod_pair = fuelPairs_.front().second;
  string in_commod = request_commod_pair.first;
  IsoVector in_recipe = request_commod_pair.second;

  // It can accept only a whole batch
  double requestAmt;
  double minAmt = in_recipe.mass();
  // The Recipe Reactor should ask for a batch if there isn't one in stock.
  double sto = this->stocksMass(); 
  // subtract sto from batch size to get total empty space. 
  // Hopefully the result  if (space <= 0) { is either 0 or the batch size 
  double space = minAmt - sto; // KDHFLAG get minAmt from the input ?
  // this will be a request for free stuff
  double commod_price = 0;

  if (space <= 0) {
    // don't request anything
    return;
  } else if (space <= minAmt) {
    // if empty space is less than monthly acceptance capacity
    requestAmt = space;
  } else if (space >= minAmt) {
    // empty space is more than monthly acceptance capacity
    // upper bound is the batch size minus the amount in stocks.
    requestAmt = capacity_ - sto;
  }

  LOG(LEV_INFO4, "BLWR1G") << " making requests {";

  MarketModel* market = MarketModel::marketForCommod(in_commod);
  Communicator* recipient = dynamic_cast<Communicator*>(market);

  // request a generic resource
  gen_rsrc_ptr request_res = gen_rsrc_ptr(new GenericResource(in_commod, "kg", requestAmt));

  // build the transaction and message
  Transaction trans(this, REQUEST);
  trans.setCommod(in_commod);
  trans.setMinFrac( minAmt/requestAmt );
  trans.setPrice(commod_price);
  trans.setResource(request_res);

  LOG(LEV_INFO5, "BLWR1G") << name() << " has requested " << request_res->quantity()
                           << " kg of " << in_commod << ".";
  sendMessage(recipient, trans);
  LOG(LEV_INFO4, "BLWR1G") << "}";
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
void BrightLightWaterReactor::sendMessage(Communicator* recipient, Transaction trans){
      msg_ptr msg(new Message(this, recipient, trans)); 
      msg->sendOn();
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
void BrightLightWaterReactor::makeOffers(){
  LOG(LEV_INFO4, "BLWR1G") << " making offers {";
  // MAKE OFFERS
  // decide how much to offer

  // there is no minimum amount a null facility may send
  double min_amt = 0;
  // this will be an offer of free stuff
  double commod_price = 0;

  // there are potentially many types of batch in the inventory stack
  double inv = this->inventoryMass();
  // send an offer for each material on the stack 
  string commod;
  Communicator* recipient;
  double offer_amt;
  for (deque<pair<string, mat_rsrc_ptr > >::iterator iter = inventory_.begin(); 
       iter != inventory_.end(); 
       iter ++){
    // get commod
    commod = iter->first;
    MarketModel* market = MarketModel::marketForCommod(commod);
    // decide what market to offer to
    recipient = dynamic_cast<Communicator*>(market);
    // get amt
    offer_amt = iter->second->quantity();

    // make a material to offer
    mat_rsrc_ptr offer_mat = mat_rsrc_ptr(new Material(iter->second));
    offer_mat->print();
    offer_mat->setQuantity(offer_amt);

    // build the transaction and message
    Transaction trans(this, OFFER);
    trans.setCommod(commod);
    trans.setMinFrac(min_amt/offer_amt);
    trans.setPrice(commod_price);
    trans.setResource(offer_mat);

    LOG(LEV_INFO5, "BLWR1G") << name() << " has offered " << offer_mat->quantity()
                             << " kg of " << commod << ".";

    sendMessage(recipient, trans);
  }
  LOG(LEV_INFO4, "BLWR1G") << "}";
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
void BrightLightWaterReactor::handleTock(int time) {
  LOG(LEV_INFO3, "BLWR1G") << name() << " is tocking {";
  // at the end of the cycle
  if (month_in_cycle_ > cycle_length_){
    this->endCycle();
  };

  // check what orders are waiting, 
  while(!ordersWaiting_.empty()){
    msg_ptr order = ordersWaiting_.front();
    order->trans().approveTransfer();
    ordersWaiting_.pop_front();
  };
  month_in_cycle_++;
  LOG(LEV_INFO3, "BLWR1G") << "}";
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    

int BrightLightWaterReactor::cycleLength() {
  return cycle_length_;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
void BrightLightWaterReactor::setCycleLength(int length) {
  cycle_length_ = length;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
double BrightLightWaterReactor::capacity() {
  return capacity_;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
void BrightLightWaterReactor::setCapacity(double cap) {
  capacity_ = cap;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
double BrightLightWaterReactor::inventorySize() {
  return inventory_size_;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
void BrightLightWaterReactor::setInventorySize(double size) {
  inventory_size_ = size;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
int BrightLightWaterReactor::facLife() {
  return lifetime_;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
void BrightLightWaterReactor::setFacLife(int lifespan) {
  lifetime_ = lifespan;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
double BrightLightWaterReactor::capacityFactor() {
  return CF_;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
void BrightLightWaterReactor::setCapacityFactor(double cf) {
  CF_ = cf;
}


//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
double BrightLightWaterReactor::inventoryMass(){
  double total = 0;

  // Iterate through the inventory and sum the amount of whatever
  // material unit is in each object.

  for (deque< pair<string, mat_rsrc_ptr> >::iterator iter = inventory_.begin(); 
       iter != inventory_.end(); 
       iter ++){
    total += iter->second->quantity();
  }

  return total;
}

double BrightLightWaterReactor::stocksMass(){
  double total = 0;

  // Iterate through the stocks and sum the amount of whatever
  // material unit is in each object.
  if(!stocks_.empty()){
    for (deque< pair<string, mat_rsrc_ptr> >::iterator iter = stocks_.begin(); 
         iter != stocks_.end(); 
         iter ++){
        total += iter->second->quantity();
    };
  };
  return total;
}


extern "C" Model* constructBrightLightWaterReactor1G() {
  return new BrightLightWaterReactor1G();
}

extern "C" void destructBrightLightWaterReactor1G(Model* model) {
  delete model;
}
