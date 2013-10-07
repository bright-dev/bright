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
              "sigma_deut": 104,
              "sigma_trit": 105,
              "sigma_He3": 106,
              "sigma_alpha": 107,
              }
mt_tallies.update({"sigma_i{0}".format(i): i+50 for i in range(1, 41)})

# Tallies, which if they exist, sum to make the absorption XS
sigma_a_tallies = set(["sigma_f", "sigma_gamma", "sigma_gamma_x", "sigma_proton", 
                       "sigma_deut", "sigma_trit", "sigma_He3", "sigma_alpha"])


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
#    "sigma_f": -6, 
    "nubar_sigma_f": -7,  
    "sigma_s_gh": None,
    "chi": None,
    "sigma_gamma_x": None,
    "sigma_2n_x": None,
    })

serpent_no_scattering = set(["sigma_t", "sigma_f", "nubar_sigma_f", "sigma_a"])
serpent_basic = set(["sigma_t", "sigma_f", "nubar_sigma_f", "sigma_a", "sigma_e", "sigma_i"])
serpent_advanced = set(["sigma_t", "sigma_f", "nubar_sigma_f", "sigma_a", "sigma_e", "sigma_i", 
                        "sigma_2n", "sigma_3n", "sigma_gamma", "sigma_proton", "sigma_alpha"])

serpent_default = set([
    "sigma_t", 
    "sigma_f", 
    "nubar_sigma_f",
    "sigma_a", 
    "sigma_e", 
    "sigma_i1", 
    "sigma_i2",
    "sigma_i3", 
    "sigma_i4", 
    "sigma_i5",
    "sigma_s_gh",
    "chi",
    "sigma_gamma", 
    "sigma_2n", 
    "sigma_3n", 
    "sigma_alpha", 
    "sigma_proton", 
    "sigma_deut",
    "sigma_trit",
    "sigma_He3",
    "sigma_gamma_x",
    "sigma_2n_x",
    ]) 



# The following defines tallies that are not valid (in serpent)
# even thought they appear as proper MT numbers in the ACE file.
#   key = (zzaaam, temp_flag)
#   value = set(mt)
restricted_tallies = {
    (80160, '06c'): set([16, 103, 104, 105]), 
    (110230, '06c'): set([16, 17]), 
    (340790, '06c'): set([17, 106]), 
    (360850, '06c'): set([17]), 
    (380890, '06c'): set([17, 104, 105, 106]), 
    (380900, '06c'): set([17, 105]), 
    (390910, '06c'): set([17]), 
    (400930, '06c'): set([17]), 
    (400950, '06c'): set([17]), 
    (410940, '06c'): set([17]), 
    (410950, '06c'): set([17, 106]), 
    (430990, '06c'): set([17]), 
    (441060, '06c'): set([17, 105]), 
    (461070, '06c'): set([17]), 
    (501230, '06c'): set([17, 104, 105]), 
    (501250, '06c'): set([17]), 
    (501260, '06c'): set([17, 104, 105]), 
    (511240, '06c'): set([17]), 
    (511250, '06c'): set([17]), 
    (511260, '06c'): set([17]), 
    (531290, '06c'): set([17]), 
    (551340, '06c'): set([17, 106]), 
    (551350, '06c'): set([17]), 
    (551360, '06c'): set([17]), 
    (551370, '06c'): set([17]), 
    (561330, '06c'): set([17]), 
    (561400, '06c'): set([17, 106]), 
    (611470, '06c'): set([17, 106]), 
    (621480, '06c'): set([17]), 
    (621510, '06c'): set([17]), 
    (631520, '06c'): set([17]), 
    (631540, '06c'): set([17]), 
    (631550, '06c'): set([17]), 
    (631560, '06c'): set([17]), 
    (822060, '06c'): set([17]), 
    (822070, '06c'): set([17]), 
    (822080, '06c'): set([17]), 
    (832090, '06c'): set([17]), 
    (882260, '06c'): set([17]), 
    (892270, '06c'): set([17]), 
    (902280, '06c'): set([17]), 
    (902290, '06c'): set([17]), 
    (902300, '06c'): set([17]), 
    (902320, '06c'): set([17]), 
    (912310, '06c'): set([17]), 
    (922320, '06c'): set([17, 38]), 
    (922330, '06c'): set([17]), 
    (922340, '06c'): set([17, 38]), 
    (922350, '06c'): set([17]), 
    (922360, '06c'): set([17, 21, 38]), 
    (922370, '06c'): set([17]), 
    (922380, '06c'): set([17]), 
    (922390, '06c'): set([17]), 
    (932350, '06c'): set([17]), 
    (932360, '06c'): set([17]), 
    (932370, '06c'): set([17]), 
    (932380, '06c'): set([17]), 
    (932390, '06c'): set([17]), 
    (942360, '06c'): set([17]), 
    (942370, '06c'): set([17]), 
    (942380, '06c'): set([17]), 
    (942390, '06c'): set([17]), 
    (942400, '06c'): set([17, 21]), 
    (942410, '06c'): set([17]), 
    (942420, '06c'): set([17]), 
    (942430, '06c'): set([17]), 
    (942440, '06c'): set([17]), 
    (942460, '06c'): set([17]), 
    (952410, '06c'): set([17]), 
    (952420, '06c'): set([17]), 
    (952421, '06c'): set([17]), 
    (952430, '06c'): set([17]), 
    (952440, '06c'): set([17]), 
    (952441, '06c'): set([17]), 
    (962420, '06c'): set([17]), 
    (962410, '06c'): set([17, 20]), 
    (962430, '06c'): set([17, 21]), 
    (962440, '06c'): set([17]), 
    (962450, '06c'): set([17]), 
    (962460, '06c'): set([17, 21]), 
    (962470, '06c'): set([17]), 
    (962480, '06c'): set([17]), 
    (962490, '06c'): set([17]), 
    (962500, '06c'): set([17]), 
    (972490, '06c'): set([17]), 
    (982490, '06c'): set([17]), 
    (982500, '06c'): set([17]), 
    (982510, '06c'): set([17]), 
    (982520, '06c'): set([17]), 
    }
