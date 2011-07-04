from bright import ReactorMG, bright_config


class OrigenReactorMG(ReactorMG):
    """A multi-group reactor that replaces the transmutation step with origen."""

    def __init__(self, *args, **kwargs):
        super(OrigenReactorMG, self).__init__(*args, **kwargs)

    def burnup_core(self):
        """Overrides the burnup core functio so that we may intercept certain methods."""

        # prep the core
        self.init_core()

        # Loop through all time steps
        for s in range(self.S):
            # Set the current time
            self.bt_s = s
            self.burn_time = self.burn_times[s]

            if (2 <= bright_config.verbosity):
                print "Time step {0} = {1} days".format(s, self.burn_times[s])

            # Find the nearest neightbors for this time.
            self.calc_nearest_neighbors()

            # Interpolate cross section in preparation for 
            # criticality calculation.
            self.interpolate_cross_sections()

            # Fold the mass weights for this time step
            self.calc_mass_weights()
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
        s = self.bt_s
        phi_t = self.phi_t
        phi_tg = self.phi_tg
        norm_phi = phi_tg[s] / phi_t[s]
        K = self.K

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
        
