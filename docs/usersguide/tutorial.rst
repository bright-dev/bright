.. _usersguide_tutorial:

*******************
The Bright Tutorial
*******************
Bright contains many nuclear fuel cycle components, such as a multi-component enrichment model and 
a light water reactor model. For more information about the individual components please refer to the 
:ref:`libref`.


===================================================
Example: Once-Through Nuclear Fuel Cycle Simulation
===================================================
The following walks through the creation of a small program which models a once-thorugh nuclear fuel cycle. 
This will aid in becoming familiar with the various modules and methods of bright.

First, we need to get the functionality of bright, the underlying pyne library, and the operating system 
This is done by the following code::

    from bright.api import *
    from pyne.material import Material
    import os

To produce power from a standard light-water reactor, uranium is needed. Uranium is obtained from 
the Earth's crust in its natural form. Natural Uranium has an isotope concentration of 0.0055% U-234, 
0.72% U-235, and 99.2745% U-238. Thus, the nuclear fuel cycle is started with natural uranium material. 
To represent natural uranium material using pyne, the isotope names are put 
together with their corresponding concentrations in a dictionary::

    nu = Material({"U234": 0.000055, "U235": 0.0072, "U238": 0.992745})

However, natural uranium requires a larger concentration of the U-235 isotope to make a valid light-water reactor
fuel. Increasing this concentration is known as enrichment.  The enrichment process is accomplished by an 
Enrichment component. In this particular example, the natural uranium is enriched to a concentration of 
3.6% U-235::

    leu = uranium_enrichment_defaults()
    leu.xP_j = 0.036
    enr = Enrichment(leu, "enrich")

The first line in the code snippet above gets the default parameters for uranium enrichment. 
The parameter xP_j is then changed to 3.6% since the default value is 5%. Finally, the Enrichment 
component is created with the 'low-enriched uranium' (leu) parameters desired and its is given the name "enrich." 
Now that the Enrichment component has been instantiated, the natural uranium may be enriched with the calc()
method.  This calculates the output of the Enrichment component from the input material::

    enr.calc(nu)

Now that we have the low enriched uranium with the concentration of U-235 needed to fuel a reactor, 
we need to feed it to a light water reactor (LWR). In order to create such a component, we need the path 
to the LWR data library.  Bright comes bundled with a default one-energy group LWR library.  The path 
to this file may found with the following::

    data_dir = os.getenv("BRIGHT_DATA")
    lwr_data = data_dir + "/LWR.h5"

Additionally, LWR parameter data needs to be provided in order to initialize the component. In this example, 
the parameters are taken from a default case but the target burnup (BUt) is instead set to a value 
of 40 [MWd/kgIHM]::

    lwrd = lwr_defaults()
    lwrd.BUt = 40.0

The LightWaterReactor1G component is thus instantiated::

    lwr = bright.LightWaterReactor1G(lwr_data, lwrd, "LWR")

The material that is produced by the Enrichment component can now be feed to the LightWaterReactor1G::

    lwr.calc(enr.mat_prod)

It is important to note that mat_prod is calculated only after calc() is called.  Finally, we instantiate 
a Storage component and feed it the output material of the LWR::

    st = Storage("Storage")
    st.decay_time = 5.0 * 365.25 * 24.0 * 3600.0
    st.calc(lwr.mat_prod)

Lastly, every fuel cycle component has a ``write()`` method that is used for outputting 
data to the hard disk in either text or HDF5 format. 

Thus the complete program of this nuclear fuel cycle simulation is provided below::

    from bright.api import *
    from pyne.material import Material
    import os

    # Set-up pointer to reactor database
    data_dir = os.getenv("BRIGHT_DATA")
    lwr_data = data_dir + "/LWR.h5"

    # Customize output
    bright_conf.write_text = False
    bright_conf.write_hdf5 = True
    load_track_nucs_hdf5(lwr_data)

    # Enrichment Calculation
    nu = Material({"U234": 0.000055, "U235": 0.0072, "U238": 0.992745})
    leu = uranium_enrichment_defaults()
    leu.xP_j = 0.036
    enr = Enrichment(enrich_params=leu, name="enrich")
    enr.calc(nu)
    enr.write()

    # Reactor Calculation
    lwrd = lwr_defaults()
    lwrd.BUt = 35.0
    lwrd.batches = 3
    lwr = LightWaterReactor1G(lwr_data, lwrd, "LWR")
    lwr.calc(enr.mat_prod)
    lwr.write()

    # Storage Calculation
    st = Storage("Storage")
    st.decay_time = 5.0 * 365.25 * 24.0 * 3600.0
    st.calc(lwr.mat_prod)
    st.write()

