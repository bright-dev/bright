.. _usersguide_custom:

*************************
Writing Custom Components
*************************
One of the main purposes of the bright framework in addition to providing a suite of fuel cycle components
is to allow the user or developer to easily create their own custom components.  The natural way to leverage 
the existing functionality is to subclass existing components.  This may be done *either* on the C++ level
or on the Python level.

============
FCComp Mixin
============
The abstract base class for bright is the fuel cycle component :class:`FCComp <bright.fccomp.FCComp>`.
This class handles much of the logic for initialization, I/O, and contains member variables which
are common to all components.  To successfully subclass :class:`FCComp <bright.fccomp.FCComp>`, 
the child class should override the following attributes and methods:

* **track_params**: A set of strings used as keys in the params_prior_calc and params_after_calc
  maps.  This denotes what parameters are fundementally important for this component to follow
  and determines which parameters are written to output.
* **calc([input])**: A method which computes and returns the output material mat_prod given an input
  material mat_feed.  If input is not supplied as an argument to this function, the mat_feed material
  which is currently on the component is used.
* **calc_params()**: A method which fills params_prior_calc and params_after_calc with the 
  relevant values.  Should only be called after calc() has executed.

The exact semantics of subclassing depend on which language is used.  For conciseness Python will be
used here.

-------
Example
-------
Suppose that a new type of cladding material which can only be produced in outer space has just 
been discovered!  This rare substance has been dubbed Gundanium alloy and is comprised of silver 
and neodymium in equal parts (AgNd).  A very expensive fabrication facility has been constructed.
As it turns out, the amount of Gundanium that may be produced is a strong function of the number 
of g-forces.  As a plucky young fuel cycle student, you have decided to model this facility.

You start out by subclassing :class:`FCComp <bright.fccomp.FCComp>`::

    pass


