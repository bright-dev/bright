# General sets
mt_tallies = {"sigma_t": 1, 
              "sigma_e": 2, 
              "sigma_i": 4, 
              "sigma_2n": 16, 
              "sigma_3n": 17, 
              "sigma_f": 18, 
              "sigma_a": 27, 
              "sigma_gamma": 102, 
              "sigma_proton": 103,
              "sigma_d": 104,
              "sigma_t": 105,
              "sigma_He3": 106,
              "sigma_alpha": 107,
              }
mt_map.update({"sigma_i{0}".format(i): i+50 for i in range(1, 41)})


# MCNP sets
mcnp_tallies = {}
mcnp_tallies.update(mt_tallies)
mcnp_tallies.update({
    "sigma_t": -1, 
    "sigma_f": -2, 
    "nubar": -3, 
    "chi": -4, 
    "sigma_a": -5,
    })

mcnp_no_scattering = set(["sigma_t", "sigma_f", "nubar", "chi", "sigma_a"])
mcnp_basic = set(["sigma_t", "sigma_f", "nubar", "chi", "sigma_a", "sigma_e", "sigma_i"])
mcnp_advanced = set(["sigma_t", "sigma_f", "nubar", "chi", "sigma_a", "sigma_e", "sigma_i", 
                     "sigma_2n", "sigma_3n", "sigma_gamma", "sigma_proton", "sigma_alpha"])

# Serpent Sets
serpent_tallies = {}
serpent_tallies.update(mt_tallies)
serpent_tallies.update({
    "sigma_t": -1, 
    "sigma_a": -2, 
    "sigma_f": -6, 
    "nubar_sigma_f": -7,  
    "chi": None,
    })

serpent_no_scattering = set(["sigma_t", "sigma_f", "nubar_sigma_f", "sigma_a"])
serpent_basic = set(["sigma_t", "sigma_f", "nubar_sigma_f", "sigma_a", "sigma_e", "sigma_i"])
serpent_advanced = set(["sigma_t", "sigma_f", "nubar_sigma_f", "sigma_a", "sigma_e", "sigma_i", 
                        "sigma_2n", "sigma_3n", "sigma_gamma", "sigma_proton", "sigma_alpha"])

serpent_default = set([
    "sigma_t", 
    "sigma_f", 
    "nubar_sigma_f",  
#   "sigma_a", 
    "sigma_e", 
    "sigma_i1", 
    "sigma_i2",
    "sigma_i3", 
    "sigma_i4", 
    "sigma_i5",
    "chi",
    ])
