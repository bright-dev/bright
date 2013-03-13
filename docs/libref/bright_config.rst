.. _bright_bright_config:

********************
Bright Configuration
********************
This module provides tools for manipulating bright-namespace level configuration 
variables which may affect all fuel cycle components.

.. currentmodule:: bright.bright_config

All functionality may be found in the ``bright_config`` module::

    import bright.bright_config

====================
Configuration Object
====================

.. autoclass:: BrightConf()

    .. autoattribute:: BRIGHT_DATA
    .. autoattribute:: track_nucs
    .. autoattribute:: track_nucs_order
    .. autoattribute:: verbosity
    .. autoattribute:: write_hdf5
    .. autoattribute:: write_text
    .. autoattribute:: write_hdf5
    .. autoattribute:: output_filename

    
================
Helper Functions
================

.. autofunction:: bright_start()
.. autofunction:: load_track_nucs_hdf5(filename, datasetname="", clear=False)
.. autofunction:: load_track_nucs_text(filename, clear=False)
.. autofunction:: sort_track_nucs()


===========
Helper Data
===========

.. data:: lib

    Path to pure C++ library directory.

.. data:: includes

    Path to C++ and Cython header file directory.

.. data:: bright_data

    Path to bright data directory.
