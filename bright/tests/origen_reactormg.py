from bright import ReactorMG, bright_config
import metasci.nuke.origen as msno
import subprocess
import numpy as np


class OrigenReactorMG(ReactorMG):
    """A multi-group reactor that replaces the transmutation step with origen."""

    def __init__(self, *args, **kwargs):
        super(OrigenReactorMG, self).__init__(*args, **kwargs)
        self._nearest_neighbors = []

    def burnup_core(self):
        """Overrides the burnup core functio so that we may intercept certain methods."""

        # prep the core
        self.init_core()

        # Loop through all time steps
        #for s in range(1):
        for s in range(self.S):
            # Set the current time
            self.bt_s = s
            self.burn_time = self.burn_times[s]

            if (2 <= bright_config.verbosity):
                print "Time step {0} = {1} days".format(s, self.burn_times[s])

            # Find the nearest neightbors for this time.
            self.calc_nearest_neighbors()
            self._nearest_neighbors.append(np.array(self.nearest_neighbors))

            # Interpolate cross section in preparation for 
            # criticality calculation.
            self.interpolate_cross_sections()

            # Fold the mass weights for this time step
            self.calc_mass_weights()
            #import pdb; pdb.set_trace()
            self.fold_mass_weights()

            # Preform the criticality calulation
            self.assemble_multigroup_matrices()
            self.calc_criticality()

            # Preform the burnup calulation
            self.assemble_transmutation_matrices()

            if (s != self.S - 1):
                self.calc_transmutation()


    def assemble_transmutation_matrices(self):
        """Reassemble transmutation matrices for ORIGEN."""
        K = self.K
        s = self.bt_s
        phi_t = self.phi_t
        phi_tg = self.phi_tg
        norm_phi = phi_tg[s] / phi_t[s]

        rx_map = {'SNG':   'sigma_gamma_itg', 
                  'SN2N':  'sigma_2n_itg', 
                  'SN3N':  'sigma_3n_itg', 
                  'SNF':   'sigma_f_itg', 
                  'SNA':   'sigma_a_itg', 
                  'SNP':   'sigma_proton_itg', 
                  'SNGX':  'sigma_gamma_x_itg', 
                  'SN2NX': 'sigma_2n_x_itg',
                  }

        xs = {key: dict() for key in rx_map.keys()}
        for rx, sig_attr in rx_map.items():
            sig = getattr(self, sig_attr)
            for iso in K:
                xs[rx][iso] = (sig[iso][s] * norm_phi).sum()

        self._xs = xs


    def calc_transmutation(self):
        """Use ORIGEN as a backend to perform the transmutation calculation."""
        K = self.K
        s = self.bt_s
        T_it = self.T_it

        # Make cross section library
        msno.write_tape9(name_org='template.tape9', **self._xs)

        # Make input mass stream
        isovec = {iso: 1E3 * T_it[iso][s] for iso in K}
        msno.write_tape4(isovec)

        # Make origen input file
        ir_type = 'IRP'
        ir_time = self.burn_times[s+1] - self.burn_times[s]
        ir_value = self.specific_power
        nlb = (219, 220, 221)
        nes = (True, False, False)
        cutoff = 1E-300
        otn = [5]
        msno.write_tape5_irradiation(ir_type, ir_time, ir_value, nlb, cut_off=cutoff, out_table_nes=nes, out_table_num=otn)

        # Run origen
        rtn = subprocess.check_call("o2_therm_linux.exe", shell=True)

        # Parse origen output 
        res = msno.parse_tape6()
        outvec = {key: sum([v[-1] for v in value]) * 1E-3 for key, value in res['table_5']['nuclide']['data'].items() if (key in K) and not np.isnan(value).any()}
        nullvec = {k: 0.0 for k in K if k not in outvec}
        outvec.update(nullvec)

        # update the transmutation matrix
        sp1 = s + 1
        for iso in K:
            T_it[iso][sp1] = outvec[iso]
        self.T_it = T_it

        # update the burnup
        BU_t = self.BU_t
        deltaBU = ir_time * ir_value
        BU_t[sp1] = BU_t[s] + deltaBU
        self.BU_t = BU_t
        print "   BU_t = {0}".format(BU_t[sp1])
