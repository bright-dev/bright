*******************
Mass Stream Package
*******************
This package provides a simple, standard way to represent fuel cycle mass flows.  The ``MassStream`` objects
represent the flows (arrows) between fuel cycle component objects.  A mass stream has three main attributes:

 * ``MassStream.mass``: The current mass of the flow in the units of the problem.
 * ``MassStream.comp``: A dictionary (map) of isotopes to their normalized fractional weights.  
   Isotopic keys are given in :ref:`zzaaam <isoform>` (integer) form.
 * ``MassStream.name``: (Optional) The label or name  of the stream.

Furthermore, mass streams are supposed to represent a new (if highly specialized) numeric type. 
Thus, some of the often-used numeric 

Presented below is the MassStream package, which contains the MassStream object and its member functions and attributes.  
For a listing of isonames raw docstrings (helpful if you need to know C++ signatures), please see :doc:`MassStream_raw`.

=======================================
:mod:`MassStream` -- Mass Stream Module
=======================================
.. currentmodule:: MassStream

All functionality may be found in the ``MassStream`` package::

 import MassStream 

This package holds the important ``MassStream`` object.

.. class:: MassStream([compdict[, mass[, name]]])

   MassStream fuel cycle flow object.

   Args:
      * `compdict` (dict or str): This is the input isotopic component dictionary.
        This dictionary need not be normalized; ``MassStream`` initialization will
        automatically renormalize the stream.  Thus compdict simply is a dictionary
        of relative weights.  The keys of compdict must be integers representing
        isotopes in :ref:`zzaaam <isoform>` form.  The values are floats for 
        each isotope's weight fraction.

        If a string is provided instead of a dictionary, then ``MassStream`` will
        read in the compdict vector from a file at the string's location using 
        the :meth:`load_from_text` function.  

        If no compdict is provided, an empty MassStream object is constructed.

   Keyword Args:
      * `mass` (float): This is the mass of the new stream. If the mass provided
        is negative (default -1.0) then the mass of the new stream is calculated from 
        the sum of compdict's components before normalization.  If the mass here
        is positive or zero, then this mass overrides the calculated one.
      * `name` (str):  A string label for the stream.  Helpful for large numbers of 
        streams. Default ''.


.. method:: MassStream.load_from_hdf5(filename, groupname[, row=N-1])

    A :class:`MassStream` object may be initialized from an HDF5 file.
    The HDF5 representation of a :class:`MassStream` is a group that holds several 
    extendable array datasets.  One array is entitled "Mass" while the other datasets
    are nuclide names in LLAAAM form ("U235", "NP237", *etc*).  For example::

        File.h5 (file)
            |-- MassStream (group)
                |-- Mass (array)
                |-- H1 (array)
                |-- O16 (array)
                |-- U235 (array)
                |-- PU239 (array)
                |-- ...

    The arrays are all of length N, where each row typically represents a different 
    fuel cycle pass.  The sum of all of the nuclide arrays should sum to one, like 
    :attr:`MassStream.comp`. 

    Args:
        * `filename`  (str): Path to HDF5 file that contains the data to read in.    
        * `groupname` (str): Path to HDF5 group that represents the data. 
            In the above example, groupname = "/MassStream".    

    Keyword Args:
        * `row` (int): The index of the arrays from which to read the data.  This 
            ranges from 0 to N-1.  Defaults to the last element of the array.
            Negative indexing is allowed (row[-N] = row[0]).

    Usage:
        This function loads data into a pre-existing :class:`MassStream`.  
        Initialization is therefore a two-step process::

            ms = MassStream()
            ms.load_from_hdf5("afile.h5", "/foo/bar/ms", -3)


.. method:: MassStream.load_from_text(filename)

    A :class:`MassStream` object may be initialized from a simple text file.
    The text representation of MassStreams are nuclide identifiers in the 
    first column and mass or weight values in the second column.  For example, 
    for natural uranium::

        922340  0.000055
        U235    0.00720
        92238   0.992745

    Data in this file must be whitespace separated.  Any valid nuclide naming
    scheme may be used for any isotope.

    Args:
        * `filename`  (str): Path to HDF5 file that contains the data to read in.    

    Usage:
        This function loads data into a pre-existing :class:`MassStream`.  
        Initialization is therefore a two-step process::

            ms = MassStream()
            ms.load_from_text("natu.h5")

        This function is most often called implicitly by the
        :class:`MassStream` constructor.

.. _MassStream_Attributes:

----------------------
MassStream Attributes
----------------------

.. attribute:: MassStream.comp

    The isotopic component vector. The comp attribute is *always* normalized!  In Python, 
    this is simply a dictionary of zzaaam isotopes to float-type weight fractions. In C++, 
    comp is a ``map<int, double>``. 

.. attribute:: MassStream.mass

    Float mass value of the stream.

.. attribute:: MassStream.name

    String name of the mass stream.

Continuing with the natural uranium example::

    #Natural Uranium MassStream
    nu = MassStream.MassStream("natural_uranium.txt", 10.0, "NU")
    print(nu.comp[922350])
    nu.mass = nu.mass * 0.5
    print(nu.name, nu.mass)
    print(nu.comp[922350])

    #Output of above code
    "0.00720"
    "NU 5.0"
    "0.00720"


.. _MassStream_Methods:

------------------
MassStream Methods
------------------

.. method:: MassStream.Normalize()

    This convenience function normalizes the mass stream by setting its mass = 1.0.

.. method:: MassStream.Print()

    This prints a string representation of the MassStream to stdout.  Print is 
    particularly useful in C++.  In Python, this method simply duplicates 
    the functionality you would get from the built-in ``str()`` function.

.. method:: MassStream.multByMass()

    This function multiplies :attr:`comp` by :attr:`mass` and returns the resultant isotopic vector.

    Returns:
        * `isovec` (dict): For a MassStream ``ms``, 

          .. math:: \mbox{isovec[iso]} = \mbox{ms.comp[iso]} \times \mbox{ms.mass}


Once again, we'll use the natural uranium stream from before as an example::

    nu.Normalize()
    print( str(nu) )
    #Mass Stream: NU
    #        Mass: 1
    #        ---------
    #        U234    5.5e-05
    #        U235    0.0072
    #        U238    0.992745

    nu.mass = nu.mass * 30.0
    nu.Print()
    #Mass Stream: NU
    #        Mass: 30
    #        ---------
    #        U234    5.5e-05
    #        U235    0.0072
    #        U238    0.992745

    unnorm_nu = nu.multByMass()
    print("Unnormalized = ", unnorm_nu )
    #Unnormalized = {922340: 0.00165, 922350: 0.216, 922380: 29.782350000000001}


.. _SubStream_Methods:

-----------------
SubStream Methods
-----------------

.. method:: MassStream.getSubStreamInt(iso_list [, name])

    Grabs a subset of the mass streams and returns a new stream comprised of only
    the specified nuclides.  The elements or isotopes included in the new substream
    are determined by iso_list, which is must only contain integers. 

    The input here is seen as a suggestion and so no error is raised if a nuclide 
    is asked for via iso_list that is not present in the original mass stream.

    Args:
        * `isoname` (list): Elements and isotopes to be taken from current stream.
          Members of this list must be integers.  For example, ``[92, 942390]``
          would take all uranium atoms and Pu-239.  
        * `name` (str): The name of the substream.

     Returns:
        * `substream` (MassStream): A new mass stream object that only 
          has the members given in iso_list.  The mass of the substream
          is calculated based on the weight fraction composition and mass
          of the original mass stream.


.. method:: MassStream.getSubStreamStr(iso_list [, name])

    Grabs a subset of the mass streams and returns a new stream comprised of only
    the specified nuclides.  The elements or isotopes included in the substream
    are determined by iso_list, which is must only contain strings. 

    The input here is seen as a suggestion and so no error is raised if a nuclide 
    is asked for via iso_list that is not present in the original mass stream.

    Args:
        * `isoname` (list): Elements and isotopes to be taken from current stream.
           Members of this list must be strings.  For example, ``['U', 'PU239']``
           would take all uranium atoms and Pu-239.  
        * `name` (str): The name of the substream.

    Returns:
        * `substream` (MassStream): A new mass stream object that only 
          has the members given in iso_list.  The mass of the substream
          is calculated based on the weight fraction composition and mass
          of the original mass stream.

For example, if we just wanted to grab U-234 and Pu-239 from the natural uranium stream::

    nu.getSubStreamStr(["U234", "PU239"], "U-234")

    #Note that no Pu is present, so we are left with only U-234
    #Mass Stream: U-234
    #        Mass: 5.5e-05
    #        ---------
    #        U234    1



.. method:: MassStream.getU([name])

    Convenience method that gets the Uranium portion of a mass stream.

    Args:
        * `name` (str): The name of the substream.

    Returns:
        * `substream` (MassStream): A new mass stream object that only 
          has Uranium members.  


.. method:: MassStream.getPU([name])

    Convenience method that gets the Plutonium portion of a mass stream.

    Args:
        * `name` (str): The name of the substream.

    Returns:
        * `substream` (MassStream): A new mass stream object that only 
          has Plutonium members.  


.. method:: MassStream.getLAN([name])

    Convenience method that gets the Lanthanide portion of a mass stream.

    Args:
        * `name` (str): The name of the substream.

    Returns:
        * `substream` (MassStream): A new mass stream object that only 
          has Lanthanide members.  


.. method:: MassStream.getACT([name])

    Convenience method that gets the Actinide portion of a mass stream.

    Args:
        * `name` (str): The name of the substream.

    Returns:
        * `substream` (MassStream): A new mass stream object that only 
          has Actinide members.  


.. method:: MassStream.getTRU([name])

    Convenience method that gets the Transuranic portion of a mass stream.

    Args:
        * `name` (str): The name of the substream.

    Returns:
        * `substream` (MassStream): A new mass stream object that only 
          has Transuranic members.  


.. method:: MassStream.getMA([name])

    Convenience method that gets the Minor Actinide portion of a mass stream.

    Args:
        * `name` (str): The name of the substream.

    Returns:
        * `substream` (MassStream): A new mass stream object that only 
          has Minor Actinide members.  


.. method:: MassStream.getFP([name])

    Convenience method that gets the Fission Product portion of a mass stream.

    Args:
        * `name` (str): The name of the substream.

    Returns:
        * `substream` (MassStream): A new mass stream object that only 
          has Fission Product members. 

.. toctree::
   :hidden:

   MassStream_raw
