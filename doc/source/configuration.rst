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


---------------------
General Specification
---------------------
* **reactor** (str): A string name for the reactor.  This is used as an identified throughout
  much of char.  Only letters and underscores are valid characters.  Exmples include
  ``"lwr"``, ``"fast_reactor"``, and ``"LWR_more_perts"``.
* **burn_regions** (int): This specifies the number of annular burn regions in a fuel pin.
  The absolut minimum number is 1, however three to four are probably required for decent statistics.
  The serpent manual recomends ten. 
* **burn_time** (int or float): This is the total time, in days, that all burnup calculations should 
  be run for. As such, it is a surrogate for the fluence.  For a light water reactor, a value of 
  ``4200`` is a decent guess.  Generally, you want to overshoot rather than undershoot this value.
* **time_step** (int or float): Time step by which to increment ``burn_time``.  The code generally
  works best when ``(0 == burn_time % time_step)``.  It is usually better to guess smaller values
  for this parameter rather than larger ones.
* **email** (str): An email address to send runtime updates to, "char@zeon.com".
* **scheduler** (str): If present, this value specifices which scheduler to use.  Accepted values include 
  an empty string ``''`` to indicate no scheduler and ``"PBS"`` to use the torque scheduler.
* **number_cpus** (int): This is the number of CPUs to run the transport code on if parallel processing 
  is available.
* **cpus_per_node** (int): If being run on a cluster, this value indicates the number of processors 
  available per node.
* **verbosity** (int): Indicates how much information should be printed to stdout.  The default level 
  is zero. You may set this arbitrarily high.

