from bright.apigen import typesystem as ts

from nose.tools import assert_equal


def check_canon(t, exp):
    obs = ts.canon(t)
    assert_equal(obs, exp)
    return obs

def test_canon():
    ts.refined_types['comp_map'] = ['map', 'nucid', 'float64']
    cases = (
        ('str', 'str'),
        (['str'], ('str', 0)),
        ('f4', 'float32'),
        ('nucid', ('int32', 'nucid')),
        (['nucid'], (('int32', 'nucid'), 0)),
        (['set', 'complex'], ('set', 'complex128', 0)),
        (['map', 'nucid', 'float'], ('map', ('int32', 'nucid'), 'float64', 0)),
        ('comp_map', (('map', ('int32', 'nucid'), 'float64', 0), 'comp_map')),
        (['char', '*'], ('char', '*')),
        (['char', 42], ('char', 42)),
        (['map', 'nucid', ['set', 'nucname']], 
            ('map', ('int32', 'nucid'), ('set', ('str', 'nucname'), 0), 0)),
        (['intrange', 1, 2], ('int32', ('intrange', 
                                            ('low', 'int32', 1), 
                                            ('high', 'int32', 2)))), 
    )
    for t, exp in cases:
        yield check_canon, t, exp            # Check that the case works,
        yield check_canon, ts.canon(t), exp  # And that it actually is canonical.

    del ts.refined_types['comp_map']
