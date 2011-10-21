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
