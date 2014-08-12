"""Plugin that runs burnup-criticality calculation.
"""
from __future__ import print_function
import os

from bright.xsgen.plugins import Plugin
from bright.xsgen.utils import RunControl, NotSpecified
from bright.xsgen.openmc_origen import OpenMCOrigen

SOLVER_ENGINES = {'openmc+origen': OpenMCOrigen}

class XSGenPlugin(Plugin):

    requires = ('bright.xsgen.pre',)

    defaultrc = RunControl(
        solver=NotSpecified,
        openmc_cross_sections=NotSpecified,
        )

    rcdocs = {
        'openmc_cross_sections': 'Path to the cross_sections.xml file for OpenMC',
        'solver': ('The physics codes that are used to solve the '
                   'burnup-criticality problem and compute cross sections and '
                   'transmutation matrices.'),
        }

    def update_argparser(self, parser):
        parser.add_argument('--solver', dest='solver', help=self.rcdocs['solver'])
        parser.add_argument("--openmc-cross-sections", dest="openmc_cross_sections",
            help=self.rcdocs['openmc_cross_sections'])

    def setup(self, rc):
        self._ensure_omcxs(rc)

        # do after all other values have been setup
        if rc.solver is NotSpecified:
            raise ValueError('a solver type must be specified')
        rc.engine = SOLVER_ENGINES[rc.solver](rc)


    def execute(self, rc):
        for state in rc.states:
            lib = rc.engine.generate(state)
            for writer in rc.writers:
                writer.write(state, lib)

    #
    # ensure functions
    #

    def _ensure_omcxs(self, rc):
        if rc.openmc_cross_sections is not NotSpecified:
            rc.openmc_cross_sections = os.path.abspath(rc.openmc_cross_sections)
        elif 'CROSS_SECTIONS' in os.environ:
            rc.openmc_cross_sections = os.path.abspath(os.environ['CROSS_SECTIONS'])
            