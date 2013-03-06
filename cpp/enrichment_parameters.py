mod_ds = """Python wrapper for enrichment parameters."""

class_ds = \
"""This class is a collection of values that mirror the attributes in 
Enrichment that are required for the cascade model to run. Like 
ReactorParameters, this class takes no arguments on initialization.  An 
empty ErichmentParameters instance has all values set to zero.
"""

desc = {
    'docstrings': {
        'module': mod_ds,
        'class': class_ds,
        'attrs': {},
        },
    'attrs': {
        'j': 'nucid',
        'k': 'nucid',
        },
    }

desc['docstrings']['attrs']['alpha_0'] = \
"""The :math:`\\alpha_0` attribute specifies the overall stage separation factor
for the cascade.  This should be set on initialization.  Values should be
greater than one.  Values less than one represent de-enrichment.
""" 

desc['docstrings']['attrs']['Mstar_0'] = \
"""The :math:`M^*_0` represents a first guess at what the `Mstar` should be.
The value of Mstar_0 on initialization should be in the ballpark
of the optimized result of the Mstar attribute.  However, :math:`M^*_0` must
always have a value between the weights of the j and k key components.
"""

desc['docstrings']['attrs']['j'] = \
"""This is an integer in zzaaam-form that represents the jth key component.
This nuclide is preferentially enriched in the product stream.
For standard uranium cascades j is 922350 (ie U-235).
"""
desc['docstrings']['attrs']['k'] = \
"""This is an integer in zzaaam-form that represents the kth key component.
This nuclide is preferentially enriched in the waste stream.
For standard uranium cascades k is 922380 (ie U-238).
"""

desc['docstrings']['attrs']['N0'] = \
"""This is the number of enriching stages initially guessed by the user."""

desc['docstrings']['attrs']['M0'] = \
"""This is the number of stripping stages initially guessed by the user."""

desc['docstrings']['attrs']['xP_j'] = \
"""This is the target enrichment of the jth isotope in the
product stream mat_prod.  The :math:`x^P_j` value is set by 
the user at initialization or run-time.  For typical uranium 
vectors, this value is about U-235 = 0.05.
"""

desc['docstrings']['attrs']['xW_j'] = \
"""This is the target enrichment of the jth isotope in the
waste stream ms_tail.  The :math:`x^W_j` value is set by the 
user at initialization or runtime.  For typical uranium vectors,
this value is about U-235 = 0.0025.
"""

