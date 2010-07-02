************
FCComp Class
************
This is the central class to all Bright fuel cycle components.  The premise of `FCComp` is that all components take a 
mass stream input, perform some operation on the input, and most return a mass stream as output.  This class, therefore,  
defines the elements common to all fuel cycle components.  As a base class it is not meant to be called directly, except 
in a child-class's initializer.

The following documentation is provided as needed information to produce your own components and as a description of 
how other fuel cycle objects behave.

.. currentmodule:: BriPy

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

.. attribute:: FCComp.IsosIn

    :attr:`IsosIn` is a :doc:`MassStream <MassStream>` object that represents the flow
    of material into this component for this pass.  

    This attribute may be accessed and altered directly (public).


.. attribute:: FCComp.IsosOut

    :attr:`IsosOut` is a :doc:`MassStream <MassStream>` object that represents the flow
    of material out of this component for this pass.  
        
    This attribute may be accessed and altered directly (public).

    .. note:: Calling :meth:`FCComp.doCalc` should calculate :attr:`IsosOut` from the :attr:`IsosIn` values.


.. attribute:: FCComp.ParamsIn

    This is a dictionary (or C map) that represents component specific parameters at input for this pass.  
    The keys are restricted to strings while their associated values are floats (doubles).  For example, 
    A reprocessing component may record the total mass input using a "Mass" key.  
        
    This attribute may be accessed and altered directly (public).

.. attribute:: FCComp.ParamsOut

    This is a dictionary (or C map) that represents component specific parameters at output for this pass.  
    The keys are restricted to strings while their associated values are floats (doubles).  For example, 
    A reactor components have a "BUd" key that represents the discharge burnup that the input
    fuel achieved.

    This attribute may be accessed and altered directly (public).

        
    .. note::

        The :attr:`ParamsIn` and :attr:`ParamsOut` attributes do not have meaningful values until :meth:`FCComp.setParams` is
        called.  This should be done after :meth:`FCComp.doCalc` but prior to output.

.. attribute:: FCComp.PassNum

    An integer representing the number of passes this component has been cycled through.  
    It starts at zero and is incremented by one each time ``FCComp.writeIsoPass()`` is called.

    This attribute may be accessed and altered directly (public).

.. attribute:: FCComp.name

    The string identifier for the component.  Defaults to an empty string.

    This attribute may be accessed and altered directly (public).

.. attribute:: FCComp.params2track

    A list (C set) of strings that holds the keys of :attr:`ParamsIn` and :attr:`ParamsOut`.
    Every component type has its own set of parameters it is able to track.  This is why 
    :attr:`params2track` is a component-specific attribute, while :func:`isos2track` is a module-level
    object.

    This attribute is set during initialization and is protected.


.. _FCComp_Methods:

==============
FCComp Methods
==============

.. method:: FCComp.doCalc([input])

    This method is used to determine a component's output isotopics from its input isotopics.
    Therefore, this is typically where the bulk of a fuel cycle component's algorithm lies.
    As each component type has a distinct methodology, the :meth:`doCalc` method  needs 
    to be overridden child classes.

    This method should return :attr:`IsosOut` so that component calculations may be easily 
    daisy-chained together.

    Args:
        * `input` (dict or MassStream): If input is present, it set as the component's 
          :attr:`IsosIn`.  If input is a isotopic dictionary (zzaaam keys, float values), this
          dictionary is first converted into a MassStream before being set as :attr:`IsosIn`.

    Returns:
        * `output` (MassStream): :attr:`IsosOut`.


.. method:: FCComp.initialize(paramlist[, name])

    This performs basic component initialization and in called by all ``FCComp`` constructors.
    As such this function is protected and will not show up in the Python bindings of Bright.

    The :meth:`initialize` method sets :attr:`params2track` and :attr:`name` as defined in the 
    function call and sets :attr:`PassNum` to zero.  Additionally, it creates and prepares isotopic 
    and parameter files for output.  These files are named "``{FCComp.name}Isos.txt``" and 
    "``{FCComp.name}Params.txt``" respectively.

    In the future, there may be an option to output in HDF5.  File initialization would take place
    here.

    .. warning:: This function will write over the isotopic and parameter files that already exist 
       in the current working directory.  Either rename the files or the components if you 
       want to ensure that data is not lost.

.. method:: FCComp.setParams()

    By calling this method, all parameter values are calculated and set for the fuel cycle component.
    This should be done following a :meth:`doCalc` calculation but before data is written out.
    If a component has important parameters associated with it, this function must be overridden and called.

    Note that this is called first thing when `writeParamPass` is called.  For example, reprocessing only 
    has a "Mass" parameter.  Translated into Python, :meth:`setParams` here looks like the following::

        def setParams(self):
            self.ParamsIn["Mass"]  = self.IsosIn.mass
            self.ParamsOut["Mass"] = self.IsosOut.mass
            return


.. method:: FCComp.writeIsoPass()

    This method is responsible for adding a new pass to the output file "``{FCComp.name}Isos.txt``"
    for this component.  Additionally, this is where :attr:`PassNum` is incremented by one.  This 
    function is therefore seen as the last step for any component in every pass.  Further calculations 
    should not be performed after :meth:`writeIsoPass` has been called.

    This function has one very important subtlety: it does not write out mass streams data.
    Rather, input columns are given as normalized isotopic vectors (:attr:`MassStream.MassStream.comp`).
    As weight fractions, input columns are in units of ``[kgInIso/kgIsosIn.mass]``.
    Moreover, the output columns are given in terms relative to the mass of the input mass, 
    ``[kgOutIso/kgIsosIn.mass]``.  These are calculated via the following expressions.

    .. math::

        \mbox{inpcol[iso]} = \mbox{IsosIn.comp[iso]}

        \mbox{outcol[iso]} = \mbox{IsosOut.comp[iso]} \times \frac{\mbox{IsosOut.mass}}{\mbox{IsosIn.mass}}

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

.. method:: FCComp.writeParamPass()

    What :meth:`writeIsoPass` does for a component's input and output isotopics, this function does for 
    the components parameters.  To ensure that meaningful data is available, :meth:`writeParamPass` first 
    calls :meth:`setParams`.  Note that to get the pass numbering correct, :meth:`writeIsoPass` should
    always be called prior to this method.  The following is an example of "``{FCComp.name}Params.txt``"
    for a light water reactor spent fuel reprocessing facility::

        Param   1in             1out	
        Mass    9.985828E-01    9.975915E-01

.. method:: FCComp.writeout()

    This is a convenience function that calls :meth:`writeIsoPass`, followed by :meth:`writeParamPass` if it is
    available.  This is what is most often used to write Bright output.  

    In a future version of Bright, this may be replaced or extended to write to an HDF5 database.
