==================================
Welcome to the Char documentation!
==================================
Cross-sections Have Awesome Rates (CHAR) is a Python program that computes neutron cross-sections
as a function of time interned in a nuclear power reactor.  

The underlying model that is used is the neutron transport code `Serpent <http://montecarlo.vtt.fi/>`_ 
for all available isotopes.  If a nuclide is not present in Serpent's ACE database, physical 
cross section models are used instead.  

.. warning::
    Unfortunately, the present version of Serpent  is not sufficient to calculate the group-to-group 
    scattering cross-section for a single material.  To rectify this, the author has a development 
    patch to Serpent that Char requires to run correctly.  This patch is available upon request.  
    Additionally, the patch has been submitted upstream so hopefully future versions of Serpent will 
    include this capability.


Char currently has the following dependencies:
   #. `NumPy <http://numpy.scipy.org/>`_
   #. `SciPy <http://www.scipy.org/>`_
   #. `PyTables <http://www.pytables.org/>`_
   #. `PyNE <http://pyne.github.com/pyne/>`_
   #. `Serpent <http://montecarlo.vtt.fi/>`_
   #. `Enthought Tool Suite <http://code.enthought.com/projects/index.php>`_ (Optional, for UI)

The source code for char may be found at the
`GitHub project site <http://github.com/scopatz/char>`_.
Or you may simply clone from the master using git::

    git clone git://github.com/scopatz/char.git

No Windows builds are currently available because much of the design assumes a posix environment.  
However, with some refactoring, a Windows version should be possible.  Please contact me 
if you are interested in this feature set.

--------
Contents
--------

.. toctree::
   :maxdepth: 1

   runtime
   configuration
   contact

=============
Helpful Links
=============
	
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
