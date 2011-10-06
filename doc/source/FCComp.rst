************
FCComp Class
************
This is the central class to all Bright fuel cycle components.  The premise of `FCComp` is that all components take a 
mass stream input, perform some operation on the input, and most return a mass stream as output.  This class, therefore,  
defines the elements common to all fuel cycle components.  As a base class it is not meant to be called directly, except 
in a child-class's initializer.

The following documentation is provided as needed information to produce your own components and as a description of 
how other fuel cycle objects behave.

.. currentmodule:: bright

.. class:: FCComp([paramlist, name])

    Base Fuel Cycle Component Class.

    Args:
        * `paramlist` (string list): A set of parameter names (str) that the component will track.
        * `name` (str): The name of the fuel cycle component instance.


    Note that this automatically calls the protected :meth:`initialize` C function.


.. _FCComp_Attributes:

=================
FCComp Attributes
=================

.. attribute:: FCComp.ms_feed

    :attr:`ms_feed` is a :doc:`MassStream <MassStream>` object that represents the flow
    of material into this component for this pass.  

    This attribute may be accessed and altered directly (public).


.. attribute:: FCComp.ms_prod

    :attr:`ms_prod` is a :doc:`MassStream <MassStream>` object that represents the flow
    of material out of this component for this pass.  
        
    This attribute may be accessed and altered directly (public).

    .. note:: Calling :meth:`FCComp.calc` should calculate :attr:`ms_prod` from the :attr:`ms_feed` values.


.. attribute:: FCComp.params_prior_calc

    This is a dictionary (or C map) that represents component specific parameters at input for this pass.  
    The keys are restricted to strings while their associated values are floats (doubles).  For example, 
    A reprocessing component may record the total mass input using a "Mass" key.  
        
    This attribute may be accessed and altered directly (public).

.. attribute:: FCComp.params_after_calc

    This is a dictionary (or C map) that represents component specific parameters at output for this pass.  
    The keys are restricted to strings while their associated values are floats (doubles).  For example, 
    A reactor components have a "BUd" key that represents the discharge burnup that the input
    fuel achieved.

    This attribute may be accessed and altered directly (public).

        
    .. note::

        The :attr:`params_prior_calc` and :attr:`params_after_calc` attributes do not have meaningful values until :meth:`FCComp.calc_params` is
        called.  This should be done after :meth:`FCComp.calc` but prior to output.

.. attribute:: FCComp.pass_num

    An integer representing the number of passes this component has been cycled through.  
    It starts at zero and is incremented by one each time ``FCComp.write_ms_pass()`` is called.

    This attribute may be accessed and altered directly (public).

.. attribute:: FCComp.name

    The string identifier for the component.  Defaults to an empty string.

    This attribute may be accessed and altered directly (public).

.. attribute:: FCComp.track_params

    A list (C set) of strings that holds the keys of :attr:`params_prior_calc` and :attr:`params_after_calc`.
    Every component type has its own set of parameters it is able to track.  This is why 
    :attr:`track_params` is a component-specific attribute, while :func:`track_nucs` is a module-level
    object.

    This attribute is set during initialization and is protected.


.. _FCComp_Methods:

==============
FCComp Methods
==============

.. method:: FCComp.calc([input])

    This method is used to determine a component's output isotopics from its input isotopics.
    Therefore, this is typically where the bulk of a fuel cycle component's algorithm lies.
    As each component type has a distinct methodology, the :meth:`calc` method  needs 
    to be overridden child classes.

    This method should return :attr:`ms_prod` so that component calculations may be easily 
    daisy-chained together.

    Args:
        * `input` (dict or MassStream): If input is present, it set as the component's 
          :attr:`ms_feed`.  If input is a isotopic dictionary (zzaaam keys, float values), this
          dictionary is first converted into a MassStream before being set as :attr:`ms_feed`.

    Returns:
        * `output` (MassStream): :attr:`ms_prod`.


.. method:: FCComp.initialize(paramlist[, name])

    This performs basic component initialization and in called by all ``FCComp`` constructors.
    As such this function is protected and will not show up in the Python bindings of Bright.

    The :meth:`initialize` method sets :attr:`track_params` and :attr:`name` as defined in the 
    function call and sets :attr:`pass_num` to zero.  Additionally, it creates and prepares isotopic 
    and parameter files for output.  These files are named "``{FCComp.name}Isos.txt``" and 
    "``{FCComp.name}Params.txt``" respectively.

    In the future, there may be an option to output in HDF5.  File initialization would take place
    here.

    .. warning:: This function will write over the isotopic and parameter files that already exist 
       in the current working directory.  Either rename the files or the components if you 
       want to ensure that data is not lost.

.. method:: FCComp.calc_params()

    By calling this method, all parameter values are calculated and set for the fuel cycle component.
    This should be done following a :meth:`calc` calculation but before data is written out.
    If a component has important parameters associated with it, this function must be overridden and called.

    Note that this is called first thing when `write_params_pass` is called.  For example, reprocessing only 
    has a "Mass" parameter.  Translated into Python, :meth:`calc_params` here looks like the following::

        def calc_params(self):
            self.params_prior_calc["Mass"]  = self.ms_feed.mass
            self.params_after_calc["Mass"] = self.ms_prod.mass
            return


.. method:: FCComp.write_ms_pass()

    This method is responsible for adding a new pass to the output text file 
    "``{FCComp.name}Isos.txt``" for this component.  Further calculations should
    not be performed after :meth:`write_ms_pass` has been called.

    This function has one very important subtlety: it does not write out mass streams data.
    Rather, input columns are given as normalized isotopic vectors (:attr:`MassStream.MassStream.comp`).
    As weight fractions, input columns are in units of ``[kgInIso/kgms_feed.mass]``.
    Moreover, the output columns are given in terms relative to the mass of the input mass, 
    ``[kgOutIso/kgms_feed.mass]``.  These are calculated via the following expressions.

    .. math::

        \mbox{inpcol[iso]} = \mbox{ms\_feed.comp[iso]}

        \mbox{outcol[iso]} = \mbox{ms\_prod.comp[iso]} \times \frac{\mbox{ms\_prod.mass}}{\mbox{ms\_feed.mass}}

    Because of the units of these two columns, total mass flow data may often only be recovered via the 
    a "Mass" parameter in the "``{FCComp.name}Params.txt``" file.  Here is a sample ``LWRIsos.txt`` file for a
    light water reactor for the first pass::

        Isotope 1in             1out	
        H1      0.000000E+00    0.000000E+00
        H3      0.000000E+00    8.568522E-08
        HE4     0.000000E+00    4.421615E-07
        B10     0.000000E+00    0.000000E+00
        B11     0.000000E+00    0.000000E+00
        C14     0.000000E+00    4.015091E-11
        O16     0.000000E+00    0.000000E+00
        SR90    0.000000E+00    8.221283E-04
        TC99    0.000000E+00    1.112580E-03
        CS137   0.000000E+00    1.821226E-03
        U234    0.000000E+00    2.807466E-06
        U235    4.773292E-02    8.951725E-03
        U236    0.000000E+00    6.155297E-03
        U237    0.000000E+00    1.719458E-05
        U238    9.522671E-01    9.211956E-01
        U239    0.000000E+00    6.953862E-07
        NP237   0.000000E+00    8.057270E-04
        PU238   0.000000E+00    2.842232E-04
        PU239   0.000000E+00    5.353362E-03
        PU240   0.000000E+00    2.114728E-03

.. method:: FCComp.write_params_pass()

    What :meth:`write_ms_pass` does for a component's input and output isotopics, 
    this function does for the components parameters.  To ensure that meaningful 
    data is available, :meth:`write_params_pass` first must have :meth:`calc_params` 
    called elsewhere in the program.  Note that to get the pass numbering correct, 
    :attr:`pass_num` should always be incremented prior to this method.  The 
    following is an example of "``{FCComp.name}Params.txt``" for a light water 
    reactor spent fuel reprocessing facility::

        Param   1in             1out	
        Mass    9.985828E-01    9.975915E-01

.. method:: FCComp.write_text()

    This method calls :meth:`write_ms_pass` and then, if available, calls 
    :meth:`write_params_pass`.  This is convience function for producing 
    text-based output.  However, using :meth:`write` is recommended.

.. method:: FCComp.write_hdf5()

    This method writes out the isotopic pass data to an HDF5 file. 
    Then, if available, it also writes parameter data as well.  
    Using :meth:`write` instead is recommended.

.. method:: FCComp.write()

    This is a convenience function that first increments up :attr:`pass_num`.
    Then, it checks to see if there are any parameters for this component.
    If there are, it sets the current values using :meth:`calc_params`.

    If :attr:`bright.write_hdf5` is set, then :meth:`write_hdf5` is called.

    If :attr:`bright.write_text` is set, then :meth:`write_text` is called.

    This is what is most often used to write Bright output.  Therefore it is
    seen as the last step for every component in each pass.  
