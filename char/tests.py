import nose
import numpy as np
import tables as tb

from numpy.testing import assert_array_equal

import isoname

rx_h5 = None
isos_LL = None


#
# Helper funcs
#

def _run_tests(path):
    """Runs tests on a library located at path"""
    global rx_h5, isos_LL
    rx_h5 = tb.openFile(path, 'r')
    isos_LL = [iso_LL for iso_LL in rx_h5.root.transmute_isos_LL]
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
        if names is None:
            names = ['arr1', 'arr2']
        print names[0] + ' = ' + repr(arr1)
        print names[1] + ' = ' + repr(arr2)
        msg = 'not ({0} <= {1})'.format(*names)
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


def test_basics():
    raise nose.SkipTest
    for grp in rx_h5.root:
        if is_data_group(grp):
            for arr in grp:
                a = np.array(arr)
                yield check_isnan, a
                yield check_le, np.array(0.0), a, ['zero', arr._v_pathname]

#
# Test Cross sections
#


def test_sigma_f():
    if not hasattr(rx_h5.root, 'sigma_f'):
        return

    for iso_LL in isos_LL:
        iso_zz = isoname.LLAAAM_2_zzaaam(iso_LL)

        sig_t_arr, sig_t = read_array(rx_h5.root.sigma_t, iso_LL)
        sig_f_arr, sig_f = read_array(rx_h5.root.sigma_f, iso_LL)
        nu_sig_f_arr, nu_sig_f = read_array(rx_h5.root.nubar_sigma_f, iso_LL)

        yield check_le, sig_f, sig_t, [sig_f_arr._v_pathname, sig_t_arr._v_pathname]

        if 89 <= (iso_zz%10000):
            mask = (sig_f != 0.0)
            nu = nu_sig_f[mask] / sig_f[mask]
            yield check_le, 1.0, nu, ['1.0', 'nu(' + sig_f_arr._v_pathname + ')']
            yield check_le, nu, 5.0, ['nu(' + sig_f_arr._v_pathname + ')', '5.0']
        else:
            yield check_eq, 0.0, sig_f, ['0.0', sig_f_arr._v_pathname]
            yield check_eq, 0.0, nu_sig_f, ['0.0', nu_sig_f_arr._v_pathname]
