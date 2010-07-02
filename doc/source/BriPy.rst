*********************
Bright/Python Package
*********************
The BriPy package is the central place for accessing all of Bright's features::

    from BriPy import *

This command implicitly imports the three extension modules: `isoname`, `MassStream`, and `FCComps`.
The first two are given separate packages in Python as a convenience to non-Bright programs.  
Switching between isotopic names is a common task.  A mass stream model, like the one here, is useful
in a variety of situations.  Partitioning `isoname` and `MassStream` out of `BriPy` allows for them to 
be easily imported elsewhere as needed.

However, the fuel cycle component library `FCComps` is essentially synonymous with Bright.  Therefore, the 
way to access all of the fuel cycle components is explicitly through `BriPy`.  

While `BriPy` is primarily a collection of components (under `FCComps`), there are a couple of module 
level functions and attributes that must be set for Bright to run successfully.  These are explained here.
Meanwhile, because of their potential for complexity, each fuel cycle component object is given its own page.  

.. toctree::
   :maxdepth: 3

   FCComps   
   FCComps_raw   


====================================
:mod:`BriPy` -- Bright/Python Module
====================================
.. currentmodule:: BriPy

.. function:: BrightStart()

    This function primarily serves to grab the environmental variable `BRIGHT_DATA` and set its value 
    within the extension module.  In Python, though not in C, `BRIGHT_DATA` is set automatically 
    to be the directory of the Bright/Python installation (usually within the site-packages directory).
    Also when you `import BriPy`, this function is called.  Therefore, the typical BriPy user should
    never need to explicitly use this function.  However, if you wish to override the default behavior 
    you may do so with the following::

        import os
        import BriPy

        os.putenv('BRIGHT_DATA', 'some/other/dir/')
        BriPy.BrightStart()

        #Continue with the fuel cycle code...

    The need for a `BrightStart()`-like function is that it ensures (in a language irrelevant way) that Bright
    can find its common data libraries that hold information on common cross-section data, half-lives, *etc*.


.. function:: isos2track([isolist])

    This function gets/sets the Bright global variable that determines which nuclides are tracked through
    the fuel cycle components.  Historically, each component was allowed to track its own unique set of 
    isotopes.  With more complex fuel cycles, this methodology quickly becomes redundant.  Therefore, 
    all components now refer to a common isotopic set.

    In posix, this set is empty by default.  In Windows, it only contains U-235.  Therefore 
    calling `isos2track()` becomes necessary at the beginning of most Bright codes so that the components 
    may be executed appropriately.  Moreover, this isotopic list is given and received in zzaaam form.  
    However, using `isoname` functions you can specify isotopic lists in a more natural way::

        import BriPy

        #Set isos2track
        isolist = BriPy.mixed_2_zzaaam_List(["U235", 94239, "H1", 80160])
        BriPy.isos2track(isolist)

        #Get isos2track
        BriPy.isos2track()  #returns the list [922350, 942390, 10010, 80160]

    Args:
        * `isolist` (zzaaam list): Integer list of isotopes (in zzaaam form) for fuel cycle components to track.

    Returns:
        * `None` (None):           Returns nothing if used to set isotopes.
        * `curriso` (zzaaam list): Returns current list of isotopes that are being tracked if `isolist` 
          was not set.


.. function:: verbosity([v])

    The verbosity level is a global parameter that determines how much information is displayed to stdout
    at runtime.  This is particularly useful in debugging.  By default, the verbosity is zero (corresponding
    to no extra output).  This function can be used to get/set the current level to any other positive 
    integer.  Increased levels will typically yield more output::

        import BriPy

        BriPy.verbosity()    #returns 0, need more...

        BriPy.verbosity(100) #Let me see everything!

    Args:
        * `v` (int): New verbosity level.

    Returns:
        * `None` (None): Returns nothing if `v` was set.
        * `currv` (int): Returns the current verbosity level if `v` was not set.
