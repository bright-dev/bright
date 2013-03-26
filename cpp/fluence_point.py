desc = {
    'docstrings': {
        'class': ("This class holds three simple data points that represent "
                  "a fluence point."),
        'attrs': {
            'f': "Index (int) of fluence immediately lower than the value of F.",
            'F': ("Fluence value itself (float).  In units of [neutrons/kilobarn], "
                  "abbr [n/kb]."),
            'm': ("The slope (float) dBU/dF between points f and f+1.  Has the odd "
                  "units of [MWd kb / kgIHM n]"),
            },
        },
    ('FluencePoint', ('f', 'int32', '0'), ('F', 'float64', '0.0'), 
                     ('m', 'float64', '0.0')): None,
    }


mod = {'FluencePoint': desc,
       'docstring': "Python wrapper for the fluence point.",}
