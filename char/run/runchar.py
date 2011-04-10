from __future__ import print_function

import os
import time
import subprocess

from metasci.colortext import message, failure


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

    def run_init(self):
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
            self.n_code.init_h5_flux_g()

        # Run initial isotope sensitivity calculation
        if self.env['options'].RUN_DELTAM:
            self.n_code.init_h5_deltam()


    def run_burnup(self, idx):
        """Runs the burnup portion of char.

        idx : a list of indeces that could be supplied 
              to range() or slice().
        """
        # Make sure we only run with the right strides
        ridx = idx[:2] + [self.n_code.ntimes]

        # run the burnup steps
        for n in range(*ridx):
            res, dep = self.n_code.run_burnup_pert(n)
            self.n_code.write_burnup(n, res, dep)

