***********************
Isotopic Naming Package
***********************
This package is used to convert between various isotopic naming schema.  Currently only 
three naming conventions are supported. 

.. _isoform:

 #. **zzaaam**: This type places the charge of the nucleus out front, then has three 
    spaces for the atomic mass number, and ends with a metastable flag (0 ground, 1 first excitied state).
    Uranium-235 here would be expressed as '922350'.
 #. **LLAAAM**: This is the more common and more readable notation.  The chamical symbol (one or two characters long)
    is first, followed by the atomic weight.  Lastly if the nuclide is metasstable, the letter *M* is concatanated 
    to the end.  For example, 'H-1' and 'Am242M' are both valid.  Note that isoname will always return LLAAAM form with
    the dash removed and all letters uppercase.
 #. **MCNP**: The MCNP format for entering nuclides is unfortunately non-standard.  In most ways it is similar 
    to zzaaam form, except that it lacks the metastable flag.  For information on how metastable isotopes are named, 
    please consult the MCNPX documentation for more information.

Presented below are the various isoname functions and attributes.  For a listing of isonames raw docstrings 
(helpful if you need to know C++ signatures), please see :doc:`isoname_raw`.

================================
:mod:`isoname` -- Isoname Module
================================
.. currentmodule:: isoname

All functionality may be found in the ``isoname`` package::

 import isoname 

This contains several zzaaam, LLAAAM, and MCNP converter function as well as other helpful module attributes.


-----------------------
Conversion Dictionaries
-----------------------

.. attribute:: isoname.LLzz

   Dictionary that is used to convert an elemental symbol (str) to its charge Z-number (int).
   For example::

      isoname.LLzz["HE"] = 2
      isoname.LLzz["U"]  = 92


.. attribute:: isoname.zzLL

   Dictionary that is used to convert a charge Z-number (int) to its elemental symbol (str).
   For example::

      isoname.LLzz[1]  = "H"
      isoname.LLzz[94] = "PU"



-------------------
Element Groups (LL)
-------------------
Element groups for the Lanthanides, Actinides, Transuranics, Minor Actinides, and Fission Products.

.. attribute:: isoname.LAN

.. attribute:: isoname.ACT

.. attribute:: isoname.TRU

.. attribute:: isoname.MA

.. attribute:: isoname.FP

The groups are defined as follows::

   isoname.LAN = ['CE', 'DY', 'ER', 'EU', 'GD', 'HO', 'LA', 'LU', 'ND', 'PM', 'PR', 'SM', 'TB', 'TM', 'YB']

   isoname.ACT = ['AC', 'AM', 'BH', 'BK', 'CF', 'CM', 'DB', 'DS', 'ES', 'FM', 'HS', 'LR', 'MD', 'MT', 'NO', 
                        'NP', 'PA', 'PU', 'RF', 'RG', 'SG', 'TH', 'U']

   isoname.TRU = ['AM', 'BH', 'BK', 'CF', 'CM', 'DB', 'DS', 'ES', 'FM', 'HS', 'LR', 'MD', 'MT', 'NO', 'NP', 
                        'PU', 'RF', 'RG', 'SG']

   isoname.MA  = ['AM', 'BH', 'BK', 'CF', 'CM', 'DB', 'DS', 'ES', 'FM', 'HS', 'LR', 'MD', 'MT', 'NO', 'NP', 
                        'RF', 'RG', 'SG']

   isoname.FP  = ['AG', 'AL', 'AR', 'AS', 'AT', 'AU', 'B',  'BA', 'BE', 'BI', 'BR', 'C',  'CA', 'CD', 'CE', 
                        'CL', 'CO', 'CR', 'CS', 'CU', 'DY', 'ER', 'EU', 'F',  'FE', 'FR', 'GA', 'GD', 'GE', 
                        'H',  'HE', 'HF', 'HG', 'HO', 'I',  'IN', 'IR', 'K',  'KR', 'LA', 'LI', 'LU', 'MG', 
                        'MN', 'MO', 'N',  'NA', 'NB', 'ND', 'NE', 'NI', 'O',  'OS', 'P',  'PB', 'PD', 'PM', 
                        'PO', 'PR', 'PT', 'RA', 'RB', 'RE', 'RH', 'RN', 'RU', 'S',  'SB', 'SC', 'SE', 'SI', 
                        'SM', 'SN', 'SR', 'TA', 'TB', 'TC', 'TE', 'TI', 'TL', 'TM', 'V',  'W',  'XE',  'Y', 
                        'YB', 'ZN', 'ZR']


-------------------
Element Groups (zz)
-------------------
Element groups for the Lanthanides, Actinides, Transuranics, Minor Actinides, and Fission Products.

.. attribute:: isoname.lan

.. attribute:: isoname.act

.. attribute:: isoname.tru

.. attribute:: isoname.ma

.. attribute:: isoname.fp

The groups are defined as follows::

   isoname.lan = [57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71]

   isoname.act = [89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111]

   isoname.tru = [93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111]

   isoname.ma  = [93, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111]

   isoname.fp  = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 
                     29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 
                     54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 
                     79, 80, 81, 82, 83, 84, 85, 86, 87, 88]


---------------------
From zzaaam functions
---------------------

.. function:: isoname.zzaaam_2_LLAAAM(nuc)
 
   Takes a zzaaam nuclide and returns its LLAAAM form.

   Args:
      * `nuc` (int): Input nuclide of zzaaam form.

   Returns:
      * `newnuc` (str): Output nuclide of LLAAAM form.

.. function:: isoname.zzaaam_2_MCNP(nuc)
 
   Takes a zzaaam nuclide and returns its MCNP form.

   Args:
      * `nuc` (int): Input nuclide in zzaaam form.

   Returns:
      * `newnuc` (int): Output nuclide in MCNP form.


.. autofunction:: isoname.zzaaam_2_LLAAAM_List

.. autofunction:: isoname.zzaaam_2_MCNP_List



---------------------
From LLAAAM functions
---------------------

.. function:: isoname.LLAAAM_2_zzaaam(nuc)
 
   Takes a LLAAAM nuclide and returns its zzaaam form.

   Args:
      * `nuc` (str): Input nuclide of LLAAAM form.

   Returns:
      * `newnuc` (int): Output nuclide of zzaaam form.

.. function:: isoname.LLAAAM_2_MCNP(nuc)
 
   Takes a LLAAAM nuclide and returns its MCNP form.

   Args:
      * `nuc` (str): Input nuclide in LLAAAM form.

   Returns:
      * `newnuc` (int): Output nuclide in MCNP form.


.. autofunction:: isoname.LLAAAM_2_zzaaam_List

.. autofunction:: isoname.LLAAAM_2_MCNP_List



-------------------
From MCNP functions
-------------------

.. function:: isoname.MCNP_2_zzaaam(nuc)
 
   Takes a MCNP nuclide and returns its zzaaam form.

   Args:
      * `nuc` (int): Input nuclide in MCNP form.

   Returns:
      * `newnuc` (int): Output nuclide in zzaaam form.

.. function:: isoname.MCNP_2_LLAAAM(nuc)
 
   Takes a MCNP nuclide and returns its LLAAAM form.

   Args:
      * `nuc` (int): Input nuclide of MCNP form.

   Returns:
      * `newnuc` (str): Output nuclide of LLAAAM form.


.. autofunction:: isoname.MCNP_2_zzaaam_List

.. autofunction:: isoname.MCNP_2_LLAAAM_List


--------------------
From mixed functions
--------------------

.. function:: isoname.mixed_2_zzaaam(nuc)
 
   Takes an arbirary nuclide and returns its zzaaam form.

   Args:
      * `nuc` (int or str): Input nuclide.

   Returns:
      * `newnuc` (int): Output nuclide in zzaaam form.

.. function:: isoname.mixed_2_LLAAAM(nuc)
 
   Takes an arbirary nuclide and returns its LLAAAM form.

   Args:
      * `nuc` (int or str): Input nuclide.

   Returns:
      * `newnuc` (str): Output nuclide in LLAAAM form.

.. function:: isoname.mixed_2_MCNP(nuc)
 
   Takes an arbirary nuclide and returns its MCNP form.

   Args:
      * `nuc` (int or str): Input nuclide.

   Returns:
      * `newnuc` (int): Output nuclide in MCNP form.


.. autofunction:: isoname.mixed_2_zzaaam_List

.. autofunction:: isoname.mixed_2_LLAAAM_List

.. autofunction:: isoname.mixed_2_MCNP_List


-------------------------------------
Isotope Vector Key Swapping Functions
-------------------------------------

Isotopic vectors are really just dictionaries or maps. 
The following functions provide an easy way to alter change the type
of all of the keys in such a vector while maintaining the value.
Unfortunately, this process involves a deep copy.

.. autofunction:: isoname.isovec_keys_2_zzaaam

.. autofunction:: isoname.isovec_keys_2_LLAAAM

.. autofunction:: isoname.isovec_keys_2_MCNP

-----------------------
Other isoname functions
-----------------------

.. function:: isoname.CurrentForm(nuc)

   Determines the form of a nuclide from the options zzaaam, LLAAAM, and MCNP.

   Args:
      * `nuc` (int or str): Input nuclide.

   Returns:
      * `FormFlag` (str): The form string from ["zzaaam", "LLAAAM", "MCNP"].
      * If the form cannot be determined, CurrentForm raises a RuntimeError 
        in Python and throws an IndeterminateNuclideForm exception in C++.

.. autofunction:: isoname.RearRemoveDuplicates




.. toctree::
   :hidden:

   isoname_raw
