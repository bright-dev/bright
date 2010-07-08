// MassStream.h
// Class that defines material flows for the Fuel Cycle.  Input into most FCComp objects.
// Effectively an normalized isotopic linked list with associated mass.
// Contains other functions for mixing mass streams and setting up sub streams.

#if !defined(_Bright_MassStream_)
#define _Bright_MassStream_

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h> 

#include "H5Cpp.h"

#include "bright.h"
#include "isoname.h"

// Set Type Definitions
typedef int Iso;
typedef double Weight;
typedef std::map<Iso, Weight> CompDict;
typedef CompDict::iterator CompIter;

class MassStream
{
// Parent Class for Handling Mass Streams"
private:
    double getArrayNthEntry(H5::DataSet *, int );

protected:
    //MassStream function definitions
    double get_comp_sum ();

public:
    void norm_comp_dict ();

    void load_from_hdf5 (std::string, std::string, int = -1);
    void load_from_text (char *);

    //Fundemental mass stream data
    CompDict comp;
    double mass;
    std::string name;

    //MassStream Constructors
    MassStream ();
    MassStream (CompDict, double = -1.0, std::string = "");
    MassStream (char *, double = -1.0, std::string = "");
    MassStream (std::string, double = -1.0, std::string = "");

    //Get Functions
    CompDict    get_comp() const {return comp;};
    double      get_mass() const {return mass;};
    std::string get_name() const {return name;};

    //Set Functions
    void set_comp(CompDict cd)   {comp = cd;};
    void set_mass(double d)      {mass = d; };
    void set_name(std::string s) {name = s; };

    //MassStream function definitions
    void Print (); 
    void Normalize ();
    CompDict multByMass();

    //Sub-Stream Computation
    MassStream getSubStream (std::set<int>,  std::string n = "");
    MassStream getSubStream (std::set<std::string>,  std::string n = "");
    MassStream getU   (std::string = "");
    MassStream getPU  (std::string = "");
    MassStream getLAN (std::string = "");
    MassStream getACT (std::string = "");
    MassStream getTRU (std::string = "");
    MassStream getMA  (std::string = "");
    MassStream getFP  (std::string = "");

};

//std::ostream& operator<< (std::ostream& os, MassStream& ms);
std::ostream& operator<< (std::ostream& os, MassStream ms);

//Overloaded Operators
MassStream operator+ (MassStream, double);
MassStream operator+ (double, MassStream);
MassStream operator+ (MassStream, MassStream);

MassStream operator* (MassStream, double);
MassStream operator* (double, MassStream);
MassStream operator/ (MassStream, double);


//Exceptions
class HDF5BoundsError: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Index of point is out of bounds.  Cannot read in from HDF5 File.";
  }
};

#endif
