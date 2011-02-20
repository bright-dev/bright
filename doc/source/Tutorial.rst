******************
The bright Tutorial
******************
bright (Br-eye Pie), or Bright/Python, is a set of Python bindings for the Bright nuclear fuel cycle model. Bright is a pure C++ library that models many canonical components such as reactors, storage facilities, and more. These components are then linked to one another using a Mass Stream object. Lastly, an isotopic naming module is available that conveniently converts between several standard nuclide naming schemes.

bright contains three modules. This modules are isoname, MassStream, and FCComps.
The command::

    import bright

will import the three modules. However, the isoname and MassStream are also provided as separate packages for convenience to non-Bright programs::

    import isoname
    import MassStream

==============
Isoname module
==============
This package is used to convert between various isotopic naming schema. Type the following command to get all the functionality of the Isoname package::

    import isoname

The following are the only three naming conventions currently supported by this package.

.. _isoform:

 #. **zzaaam**: This type places the charge of the nucleus out front, then has three 
    spaces for the atomic mass number, and ends with a metastable flag (0 ground, 1 first excited state).
    Uranium-235 here would be expressed as '922350'.
 #. **LLAAAM**: This is the more common and more readable notation.  The chemical symbol (one or two characters long)
    is first, followed by the atomic weight.  Lastly if the nuclide is metastable, the letter *M* is concatenated 
    to the end.  For example, 'H-1' and 'Am242M' are both valid.  Note that isoname will always return LLAAAM form with
    the dash removed and all letters uppercase.
 #. **MCNP**: The MCNP format for entering nuclides is unfortunately non-standard.  In most ways it is similar 
    to zzaaam form, except that it lacks the metastable flag.  For information on how metastable isotopes are named, 
    please consult the MCNPX documentation for more information.

For more information about the different attributes and methods of the Isoname module refer to its documentation.

=================
MassStream module
=================
This package provides a simple, standard way to represent fuel cycle mass flows.  The ``MassStream`` objects
represent the flows (arrows) between fuel cycle component objects.  A mass stream has three main attributes:

 * ``MassStream.mass``: The current mass of the flow in the units of the problem.
 * ``MassStream.comp``: A dictionary (map) of isotopes to their normalized fractional weights.  
   Isotopic keys are given in :ref:`zzaaam <isoform>` (integer) form.
 * ``MassStream.name``: (Optional) The label or name  of the stream.

Type the following command to get all the functionality of the MassStream package::

    import MassStream

For more information about the different attributes and methods of the MassStream module refer to its documentation.

==============
FCComps module
==============
This module contains all of the fuel cycle components, such as the Enrichment model and Light the Water Reaction model. For more information about the different components refer to the FCComps documentation.

===================================================
Example: Once-Through Nuclear Fuel Cycle Simulation
===================================================
In this section, we present a small program that models a nuclear fuel cycle. This will help in getting familiar with the different modules and methods of bright.

First, we need to get the functionality of bright and also some operating system functionally provided by the os module. This is done by the following code::

    import bright
    import os

To produce a nuclear reaction, Uranium is needed. Uranium is obtained from the Earth's crust in its natural form. Natural Uranium has an isotope concentration of 0.0055% U-234, 0.72% U-235, and 99.2745% U-238. Thus, the nuclear fuel cycle is started with a stream of natural Uranium. To represent natural uranium as a stream using bright, the zzaaam representation of the isotopes are put together with their corresponding concentrations in a dictionary::

    nu = bright.MassStream({922340 : 0.000055, 922350 : 0.0072, 922380: 0.992745})

However, a reaction will not take place unless the Uranium has a bigger concentration of the isotope U-235. Thus, the natural Uranium needs to be enriched. The enrichment process is accomplished by an Enrichment component. In this particular example, the Uranium is enriched to a concentration of 3.6% U-235::

    enrd = bright.UraniumEnrichmentDefaults()
    enrd.xP_j = 0.036
    enr = bright.Enrichment(enrd, "Enrichment")

The first line in the code snippet above gets the default parameters of the Enrichment component. The parameter xP_j needs to be changed to 3.6% since its default value is 5%. Finally, the Enrichment component is created with the parameters wanted and its is given the name "Enrichment." Now that the Enrichment component is created and has the parameters wanted, it is time to enrich the natural Uranium. This is done with the method calc, which calculates the output of the Enrichment component (or any component) from the input MassStream::

    calc(nu)

Now that we have the Uranium with the concentration of U-235 needed to produce a reaction, it is time to feed it to the Light Water Reactor. In order to create a Light Water Reactor, we need the path to the LWR HDF5 data library. This is done in the following line of code::

    data_dir = os.getenv("BRIGHT_DATA")
    lwr_data = data_dir + "/LWR.h5"

Also, parameter data needs to be provided in order to initialize the Light Water Reactor. In this example, the parameters BUt is set to 40::

    lwrd = bright.LWRDefaults()
    lwrd.BUt = 40.0

The Light Water Reactor is instantiated with the following line of code::

    lwr = bright.LightWaterReactor1G(lwr_data, lwrd, "LWR")

The MassStream that is produced by the Enrichment component can now be feed to the Light Water Reactor::

    lwr.calc(enr.ms_prod)

It is important to know that ms_prod is calculated after calc is called. Finally, we feed instantiate a Storage component and feed the output MassStream of the Light Water Reactor to it::

    st = bright.Storage("Storage")
    st.calc(lwr.ms_prod)

Lastly, every fuel cycle component contains a ``write()`` method that is used for outputting 
data to the hard disk in either text or HDF5 format. 

The complete program of this nuclear fuel cycle simulation is provided below::

    import bright
    import os

    # Set-up pointer to reactor database
    data_dir = os.getenv("BRIGHT_DATA")
    lwr_data = data_dir + "/LWR.h5"

    # Customize output
    bright.write_text(False)
    bright.write_hdf5(True)
    bright.load_track_isos_hdf5(lwr_data)
    
    # Enrichment Calculation
    nu = bright.MassStream({922340 : 0.000055, 922350 : 0.0072, 922380: 0.992745})
    enrd = bright.UraniumEnrichmentDefaults()
    enrd.xP_j = 0.036
    enr = bright.Enrichment(enrd, "Enrichment")
    enr.calc(nu)
    enr.write()

    # Reactor Calculation
    lwrd = bright.LWRDefaults()
    lwrd.BUt = 35.0
    lwrd.batches = 3
    lwr = bright.LightWaterReactor1G(lwr_data, lwrd, "LWR")
    lwr.calc(enr.ms_prod)
    lwr.write()

    # Storage Calculation
    st = bright.Storage("Storage")
    st.decay_time = 5.0 * 365.25 * 24.0 * 3600.0
    st.calc(lwr.ms_prod)
    st.write()
