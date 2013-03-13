==================================
Welcome to Bright's documentation!
==================================
Bright is a nuclear fuel cycle code.  Bright is a C++ library that models many canonical 
components such as reactors, storage facilities, and more. These components are then linked to one another using 
the Material object from PyNE.  

Arbitrarily hooking together fuel cycle objects in C-code is usually an unwieldy task that 
has low computational overhead.  Bright, therefore, is simply a collection of object models that 
allow another program to connect them.  The bright bindings enable this hooking to be done 
in Python or C++.  Thus by using bright, Python itself becomes the fuel cycle interpreter while allowing the 
heavy lifting of the models is performed by the faster C++ code.

Since bright fuel cycle objects will usually be called via Python, this documentation 
discusses the Python level, though it serves for both interfaces.

Bright currently has the following external C++ dependencies:

   #. `HDF5 <http://www.hdfgroup.org/HDF5/>`_
   #. `PyNE <http://pyne.github.com/pyne/>`_

From python, Bright may require these additional dependencies:

   #. `Cython <http://cython.org/>`_
   #. `NumPy <http://numpy.scipy.org/>`_
   #. `SciPy <http://www.scipy.org/>`_
   #. `ETS <http://code.enthought.com/>`_ (for GUI)
   #. `PyTables <http://www.pytables.org/>`_

The source code for Bright and bright may be found at the 
`GitHub project site <http://github.com/bright-dev/bright>`_.
Or you may simply clone the development branch using git::

    git clone git://github.com/bright-dev/bright.git

Lastly, you can join up as a Bright developer at 
`the Google groups page <http://groups.google.com/group/bright-dev>`_.

--------
Contents
--------

.. toctree::
    :maxdepth: 1

    usersguide/tutorial
    usersguide/index
    libref/index

=============
Helpful Links
=============
	
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
