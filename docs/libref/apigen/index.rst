.. _libref_apigen:

=================
API Generation
=================
At its core, Bright is a collection of C++ classes.  However, while one *could*
certainly write a fuel cycle simulator in pure C++ by linking these classes 
together, this is generally not the simplest or most elegant simulator engine.

Thus a major area of development for Bright is an automatic wrapper generator
layer which backends to other simulators.  The two main ones are Python and
Cyclus.  The infrastructure for creating these API is described here.

.. toctree::
    :maxdepth: 1

    typesystem
    autodescribe
    cythongen
    main
    utils
