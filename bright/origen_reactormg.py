import subprocess

import numpy as np
from pyne import origen22 
from pyne.material import Material

from bright import bright_conf
from bright.reactormg import ReactorMG


class OrigenReactorMG(ReactorMG):
    """A multi-group reactor that replaces the transmutation step with origen.  Daughter class of ReactorMG.

    Parameters
    ----------
    tape9 : str or dict
        A tape9 dictionary or a path to a TAPE9.INP file to be used as the default tape9 file.
    reactor_parameters : ReactorParameters or None, optional 
        A special data structure that contains information on how to setup and run the reactor.
    track_params : set of str or None, optional 
        A set of strings that represents what parameter data the reactor should store and set.  
        Different reactor types may have different characteristic parameters that are of interest.
    name : str, optional 
        The name of the reactor fuel cycle component instance.

    """

    def __init__(self, *args, **kwargs):
        tape9 = kwargs.pop('tape9', None)
        super(OrigenReactorMG, self).__init__(*args, **kwargs)
        self._nearest_neighbors = []

        if tape9 is None:
            raise ValueError("must supply a default TAPE9 file to OrigenReactorMG.")
        elif isinstance(tape9, basestring):
            self._tape9 = origen22.parse_tape9(tape9)
        elif hasattr(tape9, 'keys'):
            self._tape9 = tape9
        else:
            # assume sequence of tape9s
            self._tape9 = origen22.merge_tape9(tape9)


    def burnup_core(self):
        """Overrides the burnup_core() method to appropriately intercept helper methods."""

        # prep the core
        self.init_core()

        # Loop through all time steps
        #for s in range(1):
        for s in range(self.S):
            # Set the current time
            self.bt_s = s
            self.burn_time = self.burn_times[s]

            if (2 <= bright_conf.verbosity):
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
        """Reassemble transmutation matrices as ORIGEN TAPE9 input."""
        K = self.K
        s = self.bt_s
        phi_t = self.phi_t
        phi_tg = self.phi_tg
        norm_phi = phi_tg[s] / phi_t[s]

        rx_map = {'sigma_gamma':   'sigma_gamma_itg', 
                  'sigma_2n':      'sigma_2n_itg', 
                  'sigma_3n':      'sigma_3n_itg', 
                  'sigma_f':       'sigma_f_itg', 
                  'sigma_alpha':   'sigma_a_itg', 
                  'sigma_p':       'sigma_proton_itg', 
                  'sigma_gamma_x': 'sigma_gamma_x_itg', 
                  'sigma_2n_x':    'sigma_2n_x_itg',
                  }

        # Make an overlay for tape9 called xs which contains 
        # the same cross-section stricture as self._tape9
        xs = {nlb: dict() for nlb in self._tape9.keys() if self._tape9[nlb]['_type'] == 'xsfpy'}
        for nlb, val in xs.items():
            val['_type'] = self._tape9[nlb]['_type']
            val['_subtype'] = self._tape9[nlb]['_subtype']
            for rx in rx_map.keys():
                if rx not in self._tape9[nlb]:
                    continue
                val[rx] = {}
                sig = getattr(self, rx_map[rx])
                for nuc in self._tape9[nlb][rx]:
                    if nuc in K:
                        val[rx][nuc] = (sig[nuc][s] * norm_phi).sum()

        self._xs = xs


    def calc_transmutation(self):
        """Use ORIGEN as a backend to perform the transmutation calculation."""
        K = self.K
        s = self.bt_s
        T_it = self.T_it

        # Make origen cross section library
        t9 = origen22.merge_tape9([self._xs, self._tape9])
        origen22.write_tape9(t9)

        # Make input mass stream
        mat = Material({iso: 1E3 * T_it[iso][s] for iso in K})
        origen22.write_tape4(mat)

        # Make origen input file
        irr_type = 'IRP'
        irr_time = self.burn_times[s+1] - self.burn_times[s]
        irr_value = self.specific_power
        t5kw = {
            'decay_nlb': sorted(nlb for nlb in t9.keys() if t9[nlb]['_type'] == 'decay')[:3],
            'xsfpy_nlb': sorted(nlb for nlb in t9.keys() if t9[nlb]['_type'] == 'xsfpy')[:3],
            'out_table_nes': (True, False, False),
            'cut_off': 1E-300, 
            'out_table_num': [5],
            }
        origen22.write_tape5_irradiation(irr_type, irr_time, irr_value, **t5kw)

        # Run origen
        rtn = subprocess.check_call("o2_therm_linux.exe", shell=True)

        # Parse origen output 
        res = origen22.parse_tape6()
        #outvec = {key: sum([v[-1] for v in value]) * 1E-3 for key, value in res['table_5']['nuclide']['data'].items() if (key in K) and not np.isnan(value).any()}
        outvec = res['materials'][-1].comp
        nullvec = {k: 0.0 for k in K if k not in outvec}
        outvec.update(nullvec)

        # update the transmutation matrix
        sp1 = s + 1
        for nuc in K:
            T_it[nuc][sp1] = outvec[nuc]
        self.T_it = T_it

        # update the burnup
        BU_t = self.BU_t
        deltaBU = irr_time * irr_value
        BU_t[sp1] = BU_t[s] + deltaBU
        self.BU_t = BU_t
        print "   BU_t = {0}".format(BU_t[sp1])
