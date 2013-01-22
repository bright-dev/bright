calc_ds = \
"""calc(input=None)
This method performs the relatively simply task of multiplying the current 
input stream by the SE to form a new output stream::

    incomp  = self.mat_feed.mult_by_mass()
    outcomp = {}
    for iso in incomp.keys():
        outcomp[iso] = incomp[iso] * sepeff[iso]
    self.mat_prod = Material(outcomp)
    return self.mat_prod

Parameters
----------
input : dict or Material or None, optional 
    If input is present, it set as the component's mat_feed.  If input is a 
    isotopic dictionary (zzaaam keys, float values), this dictionary is first 
    converted into a Material before being set as mat_feed.

Returns
-------
output : Material
    mat_prod

"""


desc = {
    'docstrings': {
        'methods': {
            'calc': calc_ds,
            }
        }
    }
