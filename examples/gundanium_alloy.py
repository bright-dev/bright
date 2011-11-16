from bright import bright_conf
from bright.fccomp import FCComp
from pyne.material import Material

class GundaniumFab(FCComp):

    # Override constructor to provide a default g-force
    # Also, set track_params for the class here
    def __init__(self, g=1E-3, name=""):
        track_params = set(['g', 'mass_AgNd'])
        super(GundaniumFab, self).__init__(params=track_params, name=name)
        self.g = g

    # Override calc() method to create Gundanium mass 
    def calc(self, input=None):
        """Removes silver and neodymium from a material to make Gundanium alloy."""
        # Check the input first
        if input is None:
            pass
        elif isinstance(input, Material):
            self.mat_feed = input
        else:
            self.mat_feed = Material(input)

        feed = self.mat_feed

        ag = feed['Ag':'Cd']
        nd = feed['Nd':'Pm']

        if ag.mass <= nd.mass:
            nd.mass = ag.mass
        else:
            ag.mass = nd.mass

        agnd = ag + nd
        agnd.mass = agnd.mass * 10**(-self.g/9.8)
        agnd.name = "Gundanium Alloy"
        agnd.atoms_per_mol = 2

        self.mat_prod = agnd
        return agnd

    # Override the calc_params() to set the appropriate parameter values
    def calc_params(self):
        """Calculate fabrication parameters."""
        self.params_prior_calc['g'] = self.g
        self.params_after_calc['g'] = self.g

        self.params_prior_calc['mass_AgNd'] = 0.0
        self.params_after_calc['mass_AgNd'] = self.mat_prod.mass


if __name__ == '__main__':
    # Init the nuclides
    bright_conf.track_nucs = set(['Ag107', 'B10', 'ND144'])

    # Create an instance of the sub-class and some material
    gf = GundaniumFab(0.98, "Shangri-La")
    mat = Material({'Ag107': 10.0, 'B10': 42.0, 'ND144': 65.0})

    # Calculate the product produced
    prod = gf.calc(mat)
    print prod

    # Set the parameters and display output
    gf.calc_params()
    gf.write()
