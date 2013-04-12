// BrightLightWaterReactor1G.h
#ifndef _BRIGHT_CYCLUS_LIGHTWATERREACTOR1G_
#define _BRIGHT_CYCLUS_LIGHTWATERREACTOR1G_

#include "Logger.h"
#include "FacilityModel.h"
#include "QueryEngine.h"

#include "light_water_reactor1g.h"

/**
   @class BrightLightWaterReactor1G 
    
   This FacilityModel is a cyclus wrapper for bright LightWaterReactor1G class.
    
   @section intro Introduction 
   Place an introduction to the model here. 
    
   @section modelparams Model Parameters 
   Place a description of the required input parameters which define the 
   model implementation. 
    
   @section optionalparams Optional Parameters 
   Place a description of the optional input parameters to define the 
   model implementation. 
    
   @section detailed Detailed Behavior 
   Place a description of the detailed behavior of the model. Consider 
   describing the behavior at the tick and tock as well as the behavior 
   upon sending and receiving materials and messages. 
 */
class BrightLightWaterReactor1G : public FacilityModel  {
 public:
  BrightLightWaterReactor1G();
  virtual ~BrightLightWaterReactor1G();

  static unsigned int refcount;
  static bright::LightWaterReactor1G * engine;

    
  /**
     Initialize members related to derived module class

     @param qe a pointer to a QueryEngine object containing initialization data
   */
  virtual void initModuleMembers(QueryEngine* qe);
  
  /**
     A verbose printer for the BrightLightWaterReactor1G
   */
   virtual std::string str();

  /**
     Initializes a BrightLightWaterReactor1G object by copying the members of another.

     @param src the model from which to copy initialized members
    */
   virtual void cloneModuleMembersFrom(FacilityModel* src);

  /**
     Transacted resources are extracted through this method 
      
     @param order the msg/order for which resource(s) are to be prepared 
     @return list of resources to be sent for this order 
      
   */ 
  virtual std::vector<rsrc_ptr> removeResource(Transaction order);

  /**
     Transacted resources are received through this method 
      
     @param trans the transaction to which these resource objects belong 
     @param manifest is the set of resources being received 
   */ 
  virtual void addResource(Transaction trans,
                              std::vector<rsrc_ptr> manifest);

/* ------------------- */ 


/* --------------------
 * all COMMUNICATOR classes have these members
 * --------------------
 */
 public:
  /**
     The BrightLightWaterReactor1G should ignore incoming messages 
   */
  virtual void receiveMessage(msg_ptr msg);

/* -------------------- */


/* --------------------
 * all FACILITYMODEL classes have these members
 * --------------------
 */
 public:
  /**
     The handleTick function specific to the BrightLightWaterReactor1G. 
      
     @param time the time of the tick 
   */
  virtual void handleTick(int time);

  /**
     The handleTick function specific to the BrightLightWaterReactor1G. 
      
     @param time the time of the tock 
   */
  virtual void handleTock(int time);

/* ------------------- */ 


/* --------------------
 * _THIS_ FACILITYMODEL class has these members
 * --------------------
 */

/* ------------------- */ 

};

#endif

