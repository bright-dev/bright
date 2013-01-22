"""Top-level automatic API generators for Bright."""
import os
import argparse
from pprint import pprint

from bright.apigen import typesystem as ts
from bright.apigen.cythongen import gencpppxd, genpxd, genpyx
from bright.apigen.autodescribe import describe, merge_descriptions

CLASSES = [
    # classname, base filename, make cython bindings, make cyclus bindings
    ('FCComp', 'fccomp', False, False),
    ('Reprocess', 'reprocess', True, True),
    ]

def describe_class(classname, filename, env=None, verbose=False):
    if env is None:
        env = {}
    descs = [describe(filename + '.cpp', classname=classname, verbose=verbose)]
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
        descs.append(pydesc)
    desc = merge_descriptions(descs)
    basefilename = os.path.split(filename)[-1]
    desc['cpp_filename'] = '{0}.cpp'.format(basefilename)
    desc['header_filename'] = '{0}.h'.format(basefilename)
    desc['metadata_filename'] = '{0}.py'.format(basefilename)
    desc['pxd_filename'] = '{0}.pxd'.format(basefilename)
    desc['pyx_filename'] = '{0}.pyx'.format(basefilename)
    desc['cpppxd_filename'] = 'cpp_{0}.pxd'.format(basefilename)
    return desc


def main():
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
        env[classname] = describe_class(classname, os.path.join('cpp', fname), 
                                        env=env, verbose=ns.verbose)
        if ns.verbose:
            pprint(env[classname])

    # next, make cython bindings
    for classname, fname, mkcython, mkcyclus in CLASSES:
        if not mkcython or not ns.cython:
            continue
        print("making cython bindings for " + classname)

    # next, make cyclus bindings
    for classname, fname, mkcython, mkcyclus in CLASSES:
        if not mkcyclus or not ns.cyclus:
            continue
        print("making cyclus bindings for " + classname)


if __name__ == '__main__':
    main()
