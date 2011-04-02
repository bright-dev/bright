# General sets

total_only = {"sigma_t": 1}

no_scattering = {"sigma_t": 1, "sigma_f": 18, "sigma_a": 27}

basic = {"sigma_t": 1, "sigma_f": 18, "sigma_a": 27, "sigma_e": 2, "sigma_i": 4}

advanced = {"sigma_t": 1, "sigma_f": 18, "sigma_a": 27, "sigma_e": 2, "sigma_i": 4, 
            "sigma_2n": 16, "sigma_3n": 17, "sigma_gamma": 102, "sigma_proton": 103, "sigma_alpha": 107}

# MCNP sets

mcnp_no_scattering = {"sigma_t": -1, "sigma_f": -2, "nubar": -3, "chi": -4, "sigma_a": -5}

mcnp_basic = {"sigma_t": -1, "sigma_f": -2, "nubar": -3, "chi": -4, "sigma_a": -5, "sigma_e": 2, "sigma_i": 4}

mcnp_advanced = {"sigma_t": -1, "sigma_f": -2, "nubar": -3, "chi": -4, "sigma_a": -5, "sigma_e": 2, "sigma_i": 4, 
            "sigma_2n": 16, "sigma_3n": 17, "sigma_gamma": 102, "sigma_proton": 103, "sigma_alpha": 107}

# Serpent Sets

serpent_no_scattering = {"sigma_t": -1, "sigma_f": -6, "nubar_sigma_f": -7, "sigma_a": -2}

serpent_basic = {"sigma_t": -1, "sigma_f": -6, "nubar_sigma_f": -7,  "sigma_a": -2, "sigma_e": 2, "sigma_i": 4}

serpent_advanced = {"sigma_t": -1, "sigma_f": -6, "nubar_sigma_f": -7, "sigma_a": -2, "sigma_e": 2, "sigma_i": 4, 
            "sigma_2n": 16, "sigma_3n": 17, "sigma_gamma": 102, "sigma_proton": 103, "sigma_alpha": 107}

#serpent_default = {"sigma_f": 18, "nubar": -7,  "sigma_a": 27, "sigma_e": 2, "sigma_i": 4}

serpent_default = {
            "sigma_t": -1, 
            "sigma_f": -6, 
            "nubar_sigma_f": -7,  
#            "sigma_a": -2, 
            "sigma_e": 2, 
            "sigma_i1": 51, 
            "sigma_i2": 52,
            "sigma_i3": 53, 
            "sigma_i4": 54, 
            "sigma_i5": 55,
            }
