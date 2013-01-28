"""Top-level automatic API generators for Bright."""
import os
import argparse
from pprint import pprint
try:
    import cPickle as pickle
except ImportError:
    import pickle

from bright.apigen import typesystem as ts
from bright.apigen.cythongen import gencpppxd, genpxd, genpyx
from bright.apigen.autodescribe import describe, merge_descriptions

CLASSES = [
    # classname, base filename, make cython bindings, make cyclus bindings
    ('FCComp', 'fccomp', False, False),
    ('Reprocess', 'reprocess', True, True),
    ]

class DescriptionCache(object):
    """A quick persistent cache for descriptions from files.  
    The keys are (classname, filename) tuples.  The values are 
    (hashes-of-the-file, description-dictionary) tuples."""

    def __init__(self, cachefile=os.path.join('build', 'desc.cache')):
        self.cachefile = cachefile
        if os.path.isfile(cachefile):
            with open(cachefile, 'r') as f:
                self.cache = pickle.load(f)
        else:
            self.cache = {}

    def isvalid(self, classname, filename):
        """Boolean on whether the cach value for a (classname, filename)
        tuple matches the state of the file on the system."""
        key = (classname, filename)
        if key not in self.cache:
            return False
        cachehash = self.cache[key][0]
        with open(filename, 'r') as f:
            filestr = f.read()
        currhash = hash(filestr)
        return cachehash == currhash

    def __getitem__(self, key):
        return self.cache[key][1]  # return the description only

    def __setitem__(self, key, value):
        classname, filename = key
        with open(filename, 'r') as f:
            filestr = f.read()
        currhash = hash(filestr)
        self.cache[key] = (currhash, value)

    def __delitem__(self, key):
        del self.cache[key]

    def dump(self):
        """Writes the cache out to the filesystem."""
        with open(self.cachefile, 'w') as f:
            pickle.dump(self.cache, f, pickle.HIGHEST_PROTOCOL)


# singleton
cache = DescriptionCache()


def describe_class(classname, filename, verbose=False):
    """Returns a description dictionary for a class (called classname) 
    living in a file (called filename)."""
    # C++ description
    cppfilename = filename + '.cpp'
    if cache.isvalid(classname, cppfilename):
        cppdesc = cache[classname, cppfilename]
    else:
        cppdesc = describe(cppfilename, classname=classname, verbose=verbose)
        cache[classname, cppfilename] = cppdesc

    # python description
    if os.path.isfile(filename + '.py'):
        glbs = globals()
        locs = {}
        execfile(filename + '.py', glbs, locs)
        if 'desc' not in locs:
            pydesc = {}
        elif callable(locs['desc']):
            pydesc = eval('desc()', glbs, locs)
        else:
            pydesc = locs['desc']
    else:
        pydesc = {}

    desc = merge_descriptions([cppdesc, pydesc])
    basefilename = os.path.split(filename)[-1]
    desc['cpp_filename'] = '{0}.cpp'.format(basefilename)
    desc['header_filename'] = '{0}.h'.format(basefilename)
    desc['metadata_filename'] = '{0}.py'.format(basefilename)
    desc['pxd_filename'] = '{0}.pxd'.format(basefilename)
    desc['pyx_filename'] = '{0}.pyx'.format(basefilename)
    desc['cpppxd_filename'] = 'cpp_{0}.pxd'.format(basefilename)
    return desc

# Classes and type to preregister with the typesyetem prior to doing any code 
# generation.  
PREREGISTER_KEYS = ['name', 'cython_c_type', 'cython_cimport', 'cython_cy_type',
                    'cython_cyimport', 'cython_c2py', 'cython_py2c']
PREREGISTER_CLASSES = [
    ('Material', 'cpp_material.Material', ('pyne', 'cpp_material'), 
     'material._Material', ('pyne', 'material'), 
     ('{pytype}({var})', '{proxy_name} = {pytype}({var})'),
     ('{proxy_name} = {pytype}({var}, not isinstance({var}, {cytype}))',
      '{proxy_name}.mat_pointer[0]')),
    ]


def main():
    """Entry point for Bright API generation."""
    parser = argparse.ArgumentParser("Generates Bright API")
    parser.add_argument('--debug', action='store_true', default=False, 
                        help='build with debugging flags')
    parser.add_argument('--no-cython', action='store_false', dest='cython', 
                        default=True, help="don't make cython bindings")
    parser.add_argument('--no-cyclus', action='store_false', dest='cyclus', 
                        default=True, help="don't make cyclus bindings")
    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose', 
                        default=False, help="print more output")
    ns = parser.parse_args()

    ns.cyclus = False  # FIXME cyclus bindings don't exist yet!

    # compute all descriptions first 
    env = {}
    for classname, fname, mkcython, mkcyclus in CLASSES:
        print("parsing " + classname)
        desc = env[classname] = describe_class(classname, 
                                               os.path.join('cpp', fname), 
                                               verbose=ns.verbose)
        if ns.verbose:
            pprint(env[classname])

        print("registering " + classname)
        pxd_base = desc['pxd_filename'].rsplit('.', 1)[0]         # eg, fccomp
        cpppxd_base = desc['cpppxd_filename'].rsplit('.', 1)[0]   # eg, cpp_fccomp
        ts.register_class(classname,                              # FCComp
            cython_c_type=cpppxd_base + '.' + classname,          # cpp_fccomp.FCComp
            cython_cimport=('bright.' + cpppxd_base, classname),  # from bright.cpp_fccomp import FCComp
            cython_cy_type=pxd_base + '.' + classname,            # fccomp.FCComp   
            cython_cyimport=pxd_base,                             # fccomp
            )
    cache.dump()

    # now preregister types with the type system
    for prc in PREREGISTER_CLASSES:
        ts.register_class(**dict(zip(PREREGISTER_KEYS, prc)))

    # next, make cython bindings
    for classname, fname, mkcython, mkcyclus in CLASSES:
        if not mkcython or not ns.cython:
            continue
        print("making cython bindings for " + classname)
        # generate first, then write out to ensure this is atomic per-class
        desc = env[classname]
        cpppxd = gencpppxd(desc)
        pxd = genpxd(desc)
        pyx = genpyx(desc)
        with open(os.path.join('bright', desc['cpppxd_filename']), 'w') as f:
            f.write(cpppxd)
        with open(os.path.join('bright', desc['pxd_filename']), 'w') as f:
            f.write(pxd)
        with open(os.path.join('bright', desc['pyx_filename']), 'w') as f:
            f.write(pyx)

    # next, make cyclus bindings
    for classname, fname, mkcython, mkcyclus in CLASSES:
        if not mkcyclus or not ns.cyclus:
            continue
        print("making cyclus bindings for " + classname)


if __name__ == '__main__':
    main()
