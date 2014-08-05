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

GEOMETRY_TEMPLATE = """<?xml version="1.0"?>
<geometry>
  <cell id="1" fill="5" surfaces="1 -2 3 -4" />
  <cell id="101" universe="1" material="1" surfaces="-5" />
  <cell id="102" universe="1" material="void" surfaces="5 -6" />
  <cell id="103" universe="1" material="3" surfaces="6 -7" />
  <cell id="104" universe="1" material="2" surfaces="7" />
  <cell id="201" universe="2" material="2" surfaces="8 -8" />
  <cell id="302" universe="3" material="3" surfaces="9 -9" />
  <lattice id="5">
    <type>rectangular</type>
    <dimension>{_latt_shape0} {_latt_shape1}</dimension>
    <lower_left>0.0 0.0</lower_left>
    <width>{_latt_x_pitch} {_latt_y_pitch}</width>
    <universes>
      {lattice}
    </universes>
  </lattice>
  <surface id="1" type="x-plane" coeffs="0.0" boundary="reflective" />
  <surface id="2" type="x-plane" coeffs="{_latt_x_pitch}" boundary="reflective" />
  <surface id="3" type="y-plane" coeffs="0.0" boundary="reflective" />
  <surface id="4" type="y-plane" coeffs="{_latt_y_pitch}" boundary="reflective" />
  <surface id="5" type="z-cylinder" coeffs="0.0 0.0 {fuel_cell_radius}" />
  <surface id="6" type="z-cylinder" coeffs="0.0 0.0 {void_cell_radius}" />
  <surface id="7" type="z-cylinder" coeffs="0.0 0.0 {clad_cell_radius}" />
  <surface id="8" type="z-cylinder" coeffs="0.0 0.0 0.0" />
  <surface id="9" type="z-cylinder" coeffs="0.0 0.0 0.0" />
</geometry>
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
        rc = self.rc
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
        # geometry
        ctx['lattice'] = ctx['lattice'].strip().replace('\n', '\n      ')
        ctx['_latt_shape0'] = ctx['lattice_shape'][0]
        ctx['_latt_shape1'] = ctx['lattice_shape'][1]
        ctx['_latt_x_pitch'] = ctx['unit_cell_pitch'] * ctx['lattice_shape'][0]
        ctx['_latt_y_pitch'] = ctx['unit_cell_pitch'] * ctx['lattice_shape'][1]
        geometry = GEOMETRY_TEMPLATE.format(**ctx)
        with open(os.path.join(pwd, 'geometry.xml'), 'w') as f:
            f.write(geometry)

def _mat_to_nucs(mat):
    nucs = []
    template = '<nuclide name="{nuc}" wo="{mass}" />'
    for nuc, mass in mat.comp.items():
        nucs.append(template.format(nuc=nucname.serpent(nuc), mass=mass))
    nucs = "\n    ".join(nucs)
    return nucs
