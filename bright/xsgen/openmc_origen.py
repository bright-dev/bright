from __future__ import print_function

class OpenMCOrigen(object):
    """An that combines OpenMC for k-code calculations and ORIGEN for 
    transmutation.
    """

    def __init__(self, rc):
        self.rc = rc

    def generate(self, state):
        pass