from __future__ import print_function
import os
import sys
import subprocess

from pyne import nucname

# templates are from openmc/examples/lattice/simple

SETTINGS_TEMPLATE = """<?xml version="1.0"?>
<settings>
  <eigenvalue>
    <batches>{k_cycles}</batches>
    <inactive>{k_cycles_skip}</inactive>
    <particles>{k_particles}</particles>
  </eigenvalue>
  <source>
    <space type="box">
      <parameters>0 0 0 {unit_cell_height} {unit_cell_height} {unit_cell_height}</parameters>
    </space>
  </source>
</settings>
"""

MATERIALS_TEMPLATE = """<?xml version="1.0"?>
<materials>
  <default_xs>71c</default_xs>
  <material id="1">
    <density value="{fuel_density}" units="g/cc" />
    {_fuel_nucs}
  </material>
  <material id="2">
    <density value="{cool_density}" units="g/cc" />
    {_cool_nucs}
    <sab name="HH2O" xs="71t" />
  </material>
  <material id="3">
    <density value="{clad_density}" units="g/cc" />
    {_clad_nucs}
  </material>
</materials>
"""

class OpenMCOrigen(object):
    """An that combines OpenMC for k-code calculations and ORIGEN for 
    transmutation.
    """

    def __init__(self, rc):
        self.rc = rc
        self.statelibs = {}
        self.builddir = 'build-' + rc.reactor
        if not os.path.isdir(self.builddir):
            os.makedirs(self.builddir)

    def pwd(self, state):
        return os.path.join(self.builddir, str(hash(state)), 'omc')

    def context(self, state):
        rc = self.rc
        ctx = dict(rc._dict)
        ctx.update(zip(rc.perturbation_params, state))
        return ctx

    def generate(self, state):
        """Generates a library for a given state."""
        if state in self.statelibs:
            return self.statelibs[state]
        if state.burn_times != 0.0:
            raise ValueError("Burn must start at t=0.")
        self.openmc(state)
        return self.statelibs[state]

    def openmc(self, state):
        pwd = self.pwd(state)
        if not os.path.isdir(pwd):
            os.makedirs(pwd)
        self._make_omc_input(state)
        # run openmc
        # parse results

    def _make_omc_input(self, state):
        pwd = self.pwd(state)
        ctx = self.context(state)
        # settings
        settings = SETTINGS_TEMPLATE.format(**ctx)
        with open(os.path.join(pwd, 'settings.xml'), 'w') as f:
            f.write(settings)
        # materials
        ctx['_fuel_nucs'] = _mat_to_nucs(rc.fuel_material)
        ctx['_clad_nucs'] = _mat_to_nucs(rc.clad_material)
        ctx['_cool_nucs'] = _mat_to_nucs(rc.cool_material)
        materials = MATERIALS_TEMPLATE.format(**ctx)
        with open(os.path.join(pwd, 'materials.xml'), 'w') as f:
            f.write(materials)

def _mat_to_nucs(mat):
    nucs = []
    template = '<nuclide name="{nuc}" wo="{mass}" />'
    for nuc, mass in mat.comp:
        nucs.append(template.format(nuc=nucname.serpent(nuc), mass=mass))
    nucs = "\n    ".join(nucs)
    return nucs
