********************
Notes for Developers
********************

.. _win_devnotes:

============================
Notes for Windows Developers
============================

.. warning:: Windows development is deprecated. Enter at your own risk, ye brave souls!

The following is a recipe for getting a bright development environment setup in a Windows environment.
First, make sure that you have the following programs installed with the appropriate versions:

    #. `Python 2.7 <http://www.python.org/download/>`_
    #. A unix compatible editor, like `Editra <http://editra.org/download>`_
    #. `(Git, if you want to connect to the repository) <http://code.google.com/p/msysgit/>`_
    #. Microsoft Visual Studio (2008, not 2010)
    #. `HDF5 <http://www.hdfgroup.org/HDF5/>`_ 
    #. `NumPy <http://numpy.scipy.org/>`_
    #. `SciPy <http://scipy.org/>`_
    #. `MatPlotLib <http://matplotlib.sourceforge.net/>`_
    #. `PyNE <http://pyne.github.com/pyne/>`_
    #. `PyTables <http://www.pytables.org/>`_

You will need to add HDF5 location information to the environment. Thus you'll need to set the 
``LIB`` and ``INCLUDE`` environmental variables to something similar to the following::

    LIB     += C:\hdf5\lib
    INCLUDE += C:\hdf5\include

