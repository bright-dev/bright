.. _usersguide_custom:

*************************
Writing Custom Components
*************************
One of the main purposes of the bright framework in addition to providing a suite of fuel cycle components
is to allow the user or developer to easily create their own custom components.  The natural way to leverage 
the existing functionality is to subclass existing components.  This may be done *either* on the C++ level
or on the Python level.

==================
Subclassing FCComp
==================
The abstract base class for bright is the fuel cycle component :class:`FCComp <bright.fccomp.FCComp>`.
This class handles much of the logic for initialization, I/O, and contains member variables which
are common to all components.  To successfully subclass :class:`FCComp <bright.fccomp.FCComp>`, 
the child class should override the following attributes and methods:

* **track_params**: A set of strings used as keys in the params_prior_calc and params_after_calc
  maps.  This denotes what parameters are fundamentally important for this component to follow
  and determines which parameters are written to output.
* **calc([input])**: A method which computes and returns the output material mat_prod given an input
  material mat_feed.  If input is not supplied as an argument to this function, the mat_feed material
  which is currently on the component is used.
* **calc_params()**: A method which fills params_prior_calc and params_after_calc with the 
  relevant values.  Should only be called after calc() has executed.

The exact semantics of subclassing depend on which language is used.  For conciseness Python will be
used here.

------------------------
Example: Gundanium Alloy
------------------------
Suppose that a new type of cladding material which can only be produced in outer space has just 
been discovered!  This rare substance has been dubbed **Gundanium alloy** and is comprised of silver 
and neodymium in equal parts (AgNd).  A very expensive fabrication facility has been constructed.
As it turns out, the amount of Gundanium that may be produced is a strong function of the number 
of g-forces.  As a plucky young fuel cycle student, you have decided to model this facility.

You start out by subclassing :class:`FCComp <bright.fccomp.FCComp>`::

    from bright import bright_conf
    from bright.fccomp import FCComp
    from pyne.material import Material

    class GundaniumFab(FCComp):

        # Override constructor to provide a default g-force
        # Also, set track_params for the class here
        def __init__(self, g=1E-3, name=""):
            track_params = set(['g', 'mass_AgNd'])
            super(GundaniumFab, self).__init__(params=track_params, name=name)
            self.g = g

        # Override calc() method to create Gundanium mass 
        def calc(self, input=None):
            """Removes silver and neodymium from a material to make Gundanium alloy."""
            # Check the input first
            if input is None:
                pass
            elif isinstance(input, Material):
                self.mat_feed = input
            else:
                self.mat_feed = Material(input)

            feed = self.mat_feed

            ag = feed['Ag':'Cd']
            nd = feed['Nd':'Pm']

            if ag.mass <= nd.mass:
                nd.mass = ag.mass
            else:
                ag.mass = nd.mass

            agnd = ag + nd
            agnd.mass = agnd.mass * 10**(-self.g/9.8)
            agnd.name = "Gundanium Alloy"
            agnd.atoms_per_mol = 2

            self.mat_prod = agnd
            return agnd

        # Override the calc_params() to set the appropriate parameter values
        def calc_params(self):
            """Calculate fabrication parameters."""
            self.params_prior_calc['g'] = self.g
            self.params_after_calc['g'] = self.g

            self.params_prior_calc['mass_AgNd'] = 0.0
            self.params_after_calc['mass_AgNd'] = self.mat_prod.mass


    if __name__ == '__main__':
        # Init the nuclides
        bright_conf.track_nucs = set(['Ag107', 'B10', 'ND144'])

        # Create an instance of the sub-class and some material
        gf = GundaniumFab(0.98, "Shangri-La")
        mat = Material({'Ag107': 10.0, 'B10': 42.0, 'ND144': 65.0})

        # Calculate the product produced
        prod = gf.calc(mat)
        print prod

        # Set the parameters and display output
        gf.calc_params()
        gf.write()


----------------
Other Subclasses
----------------
Any of the other daughter classes of :class:`FCComp <bright.fccomp.FCComp>` may be subclassed
and their behavior altered.  More sophisticated components may require additional methods or 
attributes to be specified. In all cases, the three attributes above must be implemented.  

An example of subclassing that has become part of the bright suite is the 
:class:`OrigenReactorMG <bright.origen_reactormg.OrigenReactorMG>` component.  This class
inherits from the standard :class:`ReactorMG <bright.reactormg.ReactorMG>` class and swaps 
out the parent's transmutation methods with the an ORIGEN 2.2 based approach.  Please refer
to the source code for more implementation details.

----------
Adaptation
----------
Another powerful feature of this subclassing approach is the ability to adapt the existing 
classes to new use cases.  Suppose a pricing model (based on the mass of the output) is desired.
Thus all components should have an associated price() method.  Thin subclasses which mix
an adapter and the base classes can easily be defined.  For example::

    from bright.api import *

    # Adapter class

    class PriceAdapter(object):
        """I am useless on my own."""
        def price(self):
            return self.mat_prod.mass * 42.0


    # Adapted classes

    class PricedFCComp(FCComp, PriceAdapter):
        pass

    class PricedEnrichment(Enrichment, PriceAdapter):
        pass

    class PricedStorage(Storage, PriceAdapter):
        pass

