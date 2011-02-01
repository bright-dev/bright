==================
Char Configuration
==================
Char is highly configurable, which enables easy small and large changes to the reactor of interest.  
All of char's configuration parameters are stored in a file, which is by default called ``defchar.py``.
However, the path to this definition file is given as a runtime argument, so the file may have any name.

As its name implies, ``defchar.py`` is a python module that serves as a cookie jar for all of the
parameters used.  Because it is a python file and is executed at runtime, the underlying variables 
that char uses may be calculated in as complex or as simple a way as the user desires.  Here, 
only simple examples will be demonstrated.

Char also comes with a variety of default parameters.  In some cases, these defaults may be imported
directly into ``defchar.py``.  In other instances, defaults are selected automatically if a value is 
not provided.

What follows is an explanation of all of the input paramters that may be specified in a ``defchar.py``
file.
