import warnings
import nose
import numpy as np
import tables as tb

from numpy.testing import assert_array_equal, assert_array_almost_equal

from pyne import nucname

from tally_types import sigma_a_tallies

np.seterr(divide='ignore')

rx_h5 = None
nucs = None

npert = None
G = None

#
# Helper funcs
#

def _run_tests(path):
    """Runs tests on a library located at path"""
    global rx_h5, nucs, npert, G
    rx_h5 = tb.openFile(path, 'r')

    nucs = rx_h5.root.transmute_nucs_LL[:]
    npert = len(rx_h5.root.perturbations)
    G = len(rx_h5.root.energy[0]) - 1

    nose.runmodule(__name__, argv=[__file__])

    rx_h5.close()


def is_data_group(grp):
    name = grp._v_name
    return any(['sigma' in name, 
                'chi' == name, 
                'nubar' == name, 
                'Ti0' == name, 
                'hi_res' == name, 
                 ])


def read_array(grp, name):
    h5arr = getattr(grp, name)
    nparr = np.array(h5arr)
    return h5arr, nparr




#
# Test everything
#

def check_isnan(arr):
    assert not np.isnan(arr).any()


def check_le(arr1, arr2, names=None):
    cond = (arr1 <= arr2).all()
    if not cond:
        try:
            q = arr1 / arr2
            mask = np.isnan(q)
            assert_array_almost_equal(1.0, q[mask])
            assert (q[mask-True] <= 1.0000001).all()
            return
        except AssertionError:
            pass

        if names is None:
            names = ['arr1', 'arr2']
        print names[0] + ' = ' + repr(arr1)
        print names[1] + ' = ' + repr(arr2)
        msg = 'not ({0} <= {1})'.format(*names)
        print (arr1 <= arr2)
        print msg
        raise AssertionError(msg)



def check_eq(arr1, arr2, names=None):
    cond = (arr1 == arr2).all()
    if not cond:
        if names is None:
            names = ['arr1', 'arr2']
        print names[0] + ' = ' + repr(arr1)
        print names[1] + ' = ' + repr(arr2)
        msg = '{0} != {1}'.format(*names)
        print msg
        raise AssertionError(msg)


def check_array_eq(arr1, arr2, names=None):
    try:
        assert_array_equal(arr1, arr2)
    except AssertionError as e:
        msg = '{0} != {1}'.format(*names)
        print msg
        raise e


def check_array_almost_eq(arr1, arr2, names=None, decimal=6):
    try:
        assert_array_almost_equal(arr1, arr2, decimal)
    except AssertionError as e:
        try:
            q = arr2 / arr1
            assert_array_almost_equal(1.0, q, decimal)
            return
        except AssertionError:
            pass

        msg = '{0} != {1}'.format(*names)
        print msg
        raise e


def check_shape(arr, npert, G, name=None):
    cond = (arr.shape in [(npert, ), (npert, G), (npert, G, G)])
    if not cond:
        if name is None:
            name = 'arr'
        print name + '.shape = ' + repr(arr.shape)
        msg = '{0} not in {1}'.format(name, [(npert, ), (npert, G), (npert, G, G)])
        print msg
        raise AssertionError(msg)


def test_basics():
    for grp in rx_h5.root:
        if is_data_group(grp):
            for arr in grp:
                a = np.array(arr)
                yield check_isnan, a
                yield check_le, np.array(0.0), a, ['zero', arr._v_pathname]

                if ('hi_res' == grp._v_name):
                    hi_G = len(grp.energy) - 1
                    if (arr._v_name != 'energy'):
                        yield check_shape, a, npert, hi_G, arr._v_pathname
                else:
                    yield check_shape, a, npert, G, arr._v_pathname


def test_phi():
    phi_arr, phi = read_array(rx_h5.root, 'phi')
    phi_g_arr, phi_g = read_array(rx_h5.root, 'phi_g')
    yield check_shape, phi, npert, G, phi_arr._v_pathname
    yield check_shape, phi_g, npert, G, phi_g_arr._v_pathname
    yield check_array_almost_eq, 1.0, phi / phi_g.sum(axis=-1), [phi_arr._v_pathname, 'sum(' + phi_g_arr._v_pathname + ')'], 5

    if hasattr(rx_h5.root, 'hi_res'):
        phi_arr, phi = read_array(rx_h5.root.hi_res, 'phi')
        phi_g_arr, phi_g = read_array(rx_h5.root.hi_res, 'phi_g')
        yield check_array_almost_eq, 1.0, phi / phi_g.sum(axis=-1), [phi_arr._v_pathname, 'sum(' + phi_g_arr._v_pathname + ')'], 5


def test_root_misc():
    time0_arr, time0 = read_array(rx_h5.root, 'time0')
    yield check_shape, time0, npert, G, time0_arr._v_pathname

    BU0_arr, BU0 = read_array(rx_h5.root, 'BU0')
    yield check_shape, BU0, npert, G, BU0_arr._v_pathname

    Phi_arr, Phi = read_array(rx_h5.root, 'Phi')
    yield check_shape, Phi, npert, G, Phi_arr._v_pathname

    pert_arr, pert = read_array(rx_h5.root, 'perturbations')
    yield check_shape, pert, npert, G, pert_arr._v_pathname


#
# Test Cross sections
#


def test_sigma_f():
    if not hasattr(rx_h5.root, 'sigma_f'):
        raise nose.SkipTest

    for nuc in nucs:
        nuc_zz = nucname.zzaaam(nuc)

        sig_t_arr, sig_t = read_array(rx_h5.root.sigma_t, nuc)
        sig_f_arr, sig_f = read_array(rx_h5.root.sigma_f, nuc)
        nu_sig_f_arr, nu_sig_f = read_array(rx_h5.root.nubar_sigma_f, nuc)

        yield check_le, sig_f, sig_t, [sig_f_arr._v_pathname, sig_t_arr._v_pathname]

        if 86 <= (nuc_zz/10000):
            mask = (sig_f != 0.0)
            nu = nu_sig_f[mask] / sig_f[mask]
            yield check_le, 1.0, nu, ['1.0', 'nu(' + sig_f_arr._v_pathname + ')']
            yield check_le, nu, 5.5, ['nu(' + sig_f_arr._v_pathname + ')', '5.5']
        else:
            yield check_eq, 0.0, sig_f, ['0.0', sig_f_arr._v_pathname]
            yield check_eq, 0.0, nu_sig_f, ['0.0', nu_sig_f_arr._v_pathname]


def test_chi():
    if not hasattr(rx_h5.root, 'chi'):
        raise nose.SkipTest

    for nuc in nucs:    
        nuc_zz = nucname.zzaaam(nuc)

        chi_arr, chi = read_array(rx_h5.root.chi, nuc)
        sig_f_arr, sig_f = read_array(rx_h5.root.sigma_f, nuc)

        if 86 <= (nuc_zz/10000):
            yield check_array_almost_eq, 1.0, chi.sum(axis=1), ['1.0', 'sum(' + chi_arr._v_pathname + ')']
        else:
            yield check_eq, 0.0, chi, ['0.0', chi_arr._v_pathname]


def test_sigma_s():
    if not hasattr(rx_h5.root, 'sigma_s'):
        raise nose.SkipTest

    for nuc in nucs:    
        nuc_zz = nucname.zzaaam(nuc)

        sig_t_arr, sig_t = read_array(rx_h5.root.sigma_t, nuc)
        sig_s_arr, sig_s = read_array(rx_h5.root.sigma_s, nuc)
        sig_s_gh_arr, sig_s_gh = read_array(rx_h5.root.sigma_s_gh, nuc)

        yield check_le, sig_s, sig_t, [sig_s_arr._v_pathname, sig_t_arr._v_pathname]
        yield check_array_almost_eq, sig_s, sig_s_gh.sum(axis=-2), [sig_s_arr._v_pathname, 'sum(' + sig_s_gh_arr._v_pathname + ')']


def test_sigma_a():
    #raise nose.SkipTest
    if not hasattr(rx_h5.root, 'sigma_a'):
        raise nose.SkipTest

    for nuc in nucs:    
        nuc_zz = nucname.zzaaam(nuc)

        sig_t_arr, sig_t = read_array(rx_h5.root.sigma_t, nuc)
        sig_a_arr, sig_a = read_array(rx_h5.root.sigma_a, nuc)

        yield check_le, sig_a, sig_t, [sig_a_arr._v_pathname, sig_t_arr._v_pathname]

        tot_sig_a = np.zeros(sig_a.shape, dtype=float)

        for tally in sigma_a_tallies:
            if hasattr(rx_h5.root, tally):
                sig_tally_arr, sig_tally = read_array(getattr(rx_h5.root, tally), nuc)
                yield check_le, sig_tally, sig_a, [sig_tally_arr._v_pathname, sig_a_arr._v_pathname]
                tot_sig_a += sig_tally

        yield check_le, tot_sig_a, sig_a, ['sum(sig_a_parts)', sig_a_arr._v_pathname]


