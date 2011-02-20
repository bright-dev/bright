==================================
Welcome to bright's documentation!
==================================
bright (*Br-eye Pie*), or Bright/Python, is a set of Python bindings for the Bright nuclear fuel cycle model.
Bright is a pure C++ library that models many canonical components such as reactors, 
storage facilities, and more.  These components are then linked to one another using 
a Mass Stream object.  Lastly, an isotopic naming module is available that conveniently 
converts between several standard nuclide naming schemes. 

Arbitrarily hooking together fuel cycle objects in C-code is usually an unwieldy task that 
has low computational overhead.  Bright, therefore, is simply a collection of object models that 
allow another program to connect them.  The bright bindings enable this hooking to be done 
in Python.  Thus by using bright, Python itself becomes the fuel cycle interpreter while allowing the 
heavy lifting of the models to be performed by the faster C++ code.

Since Bright fuel cycle objects will usually be called via Python, this documentation 
serves for both the Bright objects and bright interface.

Bright currently has the following dependencies:
   #. `HDF5 <http://www.hdfgroup.org/HDF5/>`_
   #. `Boost <http://www.boost.org/>`_

bright requires these additional dependencies:
   #. `Boost::Python <http://www.boost.org/doc/libs/release/libs/python/doc/>`_
   #. `NumPy <http://numpy.scipy.org/>`_
   #. `SciPy <http://www.scipy.org/>`_
   #. `MatPlotLib <http://matplotlib.sourceforge.net/>`_
   #. `numexpr <http://code.google.com/p/numexpr/>`_ (now needed for Pytables).
   #. `PyTables <http://www.pytables.org/>`_

The source code for Bright and bright may be found at the 
`GitHub project site <http://github.com/bright-dev>`_.
Or you may simply branch from the trunk using git::

    git clone git://github.com/bright-dev/bright.git

Additionally, a `Windows build is provided here <http://nukestar.me.utexas.edu/scopatz/bright/BriPy-0.23.win32-py2.6.msi>`_.

Lastly, you can join up as a Bright developer at 
`the Google groups page <http://groups.google.com/group/bright-dev>`_.

--------
Contents
--------

.. toctree::
    :maxdepth: 2

    Tutorial
    isoname   
    MassStream
    BriPy
    DevNotes

=============
Helpful Links
=============
	
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
