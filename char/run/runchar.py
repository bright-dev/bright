from __future__ import print_function

import os
import time
import subprocess

import numpy as np
import tables as tb

from pyne.material import Material

from pyne.utils import message, failure

class RunChar(object):
    """A controller to run char very generally."""

    def __init__(self, n_code, env):
        """Args:
            * n_code: a neutron transport model
            * env: the environment to execute the model in.
        """
        self.n_code = n_code
        self.env = env


    #
    # Controleer functions
    # 

    def init_h5(self):
        """Inits the char library."""

        # Make a new HDF5 file.
        if (self.env['options'].MAKE_INPUT):
            self.n_code.init_h5()
            
        # the cross-section generation
        if self.env['options'].RUN_BURNUP:
            self.n_code.init_h5_burnup()

        # Make Cross-sections as a separate step from the burnup calculation
        if self.env['options'].RUN_XS_GEN:
            self.n_code.init_h5_xs_gen()

            #if self.env['xs_models_needed']:
            #    self.n_code.init_h5_flux_g()
            self.n_code.init_h5_flux_g()

        # Run initial nuclide sensitivity calculation
        if self.env['options'].RUN_DELTAM:
            self.n_code.init_h5_deltam()


    def burnup(self, idx):
        """Runs the burnup portion of char.

        idx : a list of perturbation indices that 
              could be supplied to range() or slice().
        """
        # Make sure we only run with the right strides
        ridx = idx[:2] + [self.n_code.ntimes]

        # run the burnup steps
        for n in range(*ridx):
            res, dep = self.n_code.run_burnup_pert(n)
            self.n_code.write_burnup(n, res, dep)


    def xs_gen(self, idx, nucs):
        """Runs the cross-section generation portion of char.

        idx : a list of perturbation indices that 
              could be supplied to range() or slice().
        nucs : a set of nuclides to run (zzaaam-form).
        """
        nucs_in_serpent = (nucs & set(self.env['core_transmute_in_serpent']))
        nucs_not_in_serpent = (nucs & set(self.env['core_transmute_not_in_serpent']))

        # Loop over the perturbation steps
        for n in range(*idx):
            # Grab the Material at this time.
            ms_n = Material()
            ms_n.from_hdf5(self.env['reactor'] + ".h5", "/Ti0", n, protocol=0)

            # Calc restricted mass streams
            ms_n_in_serpent = ms_n[self.env['core_transmute_in_serpent']]
            ms_n_not_in_serpent = ms_n[self.env['core_transmute_not_in_serpent']]

            # Read in some common parameters from the data file
            with tb.openFile(self.env['reactor'] + ".h5", 'r') as  rx_h5:
                E_g = np.array(rx_h5.root.energy[n][::-1])
                E_n = np.array(rx_h5.root.hi_res.energy.read()[::-1])
                phi_n = np.array(rx_h5.root.hi_res.phi_g[n][::-1])

            # Run and write the high resolution flux
            if (phi_n < 0.0).all():
                res, det = self.n_code.run_flux_g_pert(n, ms_n_in_serpent)
                self.n_code.write_flux_g(n, res, det)
                with tb.openFile(self.env['reactor'] + ".h5", 'r') as  rx_h5:
                    phi_n = np.array(rx_h5.root.hi_res.phi_g[n][::-1])

            #
            # Loop over all output nuclides...
            #
            # ...that are valid in serpent
            for nuc in nucs_in_serpent:
                res, det = self.n_code.run_xs_gen_pert(nuc, n, ms_n_in_serpent, E_n, E_g, phi_n)
                self.n_code.write_xs_gen(nuc, n, res, det)

            # ...that are NOT valid in serpent
            for nuc in nucs_not_in_serpent:
                xsd = self.n_code.run_xs_mod_pert(nuc, n, E_n, E_g, phi_n)
                self.n_code.write_xs_mod(nuc, n, xsd)


    def deltam(self, idx, nucs, sidx):
        """Runs the nuclide sensitivity study.

        idx : a list of perturbation indices that 
              could be supplied to range() or slice().
        nucs : a set of nuclides to run (zzaaaam-form).
        sidx : a list of sensitivity indices that 
              could be supplied to range() or slice().
        """
        # Make sure we only run with the right strides
        ridx = idx[:2] + [self.n_code.ntimes]

        # Loop over all perturbations.
        for n in range(*ridx):
            # Loop over all nuclides
            for nuc_zz in nucs:
                # Calulate this nuclides new values of IHM concentration
                nuc_fracs = self.env['deltam'] * self.ihm_stream.comp[nuc_zz]

                # Skip nuclides that would be pertubed over 1 kgIHM
                if (1.0 < nuc_fracs).any():
                    continue

                # Loop over all nuclide sesnitivities
                for s in range(*sidx):
                    res, dep = self.n_code.run_deltam_pert(nuc_zz, n, s, nuc_fracs)
                    self.n_code.write_deltam(nuc_zz, n, s, nuc_fracs, res, dep)

