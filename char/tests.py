import nose
import numpy as np
import tables as tb


f = None


def _run_tests(path):
    """Runs tests on a library located at path"""
    global f
    f = tb.openFile(path, 'r')
#    nose.main()
    print __file__
#    nose.run(argv=[__file__])
    nose.runmodule(__name__, argv=[__file__])
    f.close()


def test_a():
    print "I am here"
    print dir(f)
    raise TypeError
