from nose.tools import assert_equal, assert_raises


from char import utils


def test_load_nuc_file():
    expected = [10010, 80160, 952421]
    observed = utils.load_nuc_file('nuc_set.txt')
    assert_equal(expected, observed)


def test_temp_flag():
    assert_equal(utils.temperature_flag(600), '06c')
    assert_equal(utils.temperature_flag(1800.0), '18c')
    assert_raises(ValueError, utils.temperature_flag, 601)
    assert_raises(ValueError, utils.temperature_flag, -600)
    assert_raises(ValueError, utils.temperature_flag, 30000000)
