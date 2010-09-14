********************
Notes for Developers
********************

.. _Win_DevNotes:

============================
Notes for Windows Developers
============================

The following is a recipe for getting a bright development environment setup in a Windows environment.
First, make sure that you have the following programs installed with the appropriate versions:

    #. `Python 2.6 <http://www.python.org/download/>`_
    #. A unix compatible editor, like `Editra <http://editra.org/download>`_
    #. `(Git, if you want to connect to the repository) <http://code.google.com/p/msysgit/>`_
    #. Microsoft Visual Studio (2008, not 2010)
    #. Boost C++ Libraries (Should be available from your distribution's package manager.  
       `You can find window's binaries here. <http://www.boostpro.com/download/>`_)
    #. `NumPy <http://numpy.scipy.org/>`_
    #. `SciPy <http://scipy.org/>`_
    #. `MatPlotLib <http://matplotlib.sourceforge.net/>`_
    #. `MetaSci <http://nukestar.me.utexas.edu/scopatz/metasci/>`_
    #. `PyTables <http://www.pytables.org/>`_

Note that Boost puts its binaries in a non-standard place.  (Is there anything really standard about Windows?!)  However, 
to compile, you will have to point to them.  Thus you'll need to set the ``LIB`` and ``INCLUDE`` environmental variables
to somethign similar to the following::

    LIB     = C:\Program Files (x86)\boost\boost_1_44\lib
    INCLUDE = C:\Program Files (x86)\boost\boost_1_44;%INCLUDE%

Enjoy!
