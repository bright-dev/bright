from nose.tools import assert_equal


from char.char import utils


def test_load_nuc_file():
    expected = [10010, 80160, 954241]
    observed = utils.load_nuc_file('nuc_set.txt')
    assert_equal(expected, observed)
