import sys

_imported = False
_paths = ['v3_2',]
_glbs = globals()
_locs = locals()

for _path in _paths:
    try:
        _cindex = __import__('.{0}.cindex'.format(_path), _glbs, _locs, [], -1)
        _imported = True
    except _cindex.LibclangError:
        print "Failed to import " + _path
    if _imported:
        break

_currmod = sys.modules[__name__]
for a in dir(_cindex):
    if not a.startswith('_'):
        setattr(a, _currmod, getattr(a, _currmod))
