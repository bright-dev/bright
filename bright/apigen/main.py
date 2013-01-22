"""Top-level automatic API generators for Bright."""
import os
import argparse
from pprint import pprint

from autodescribe import describe, merge_descriptions

CLASSES = [
    # classname, filename, make cython bindings, make cyclus bindings
    ('FCComp', 'fccomp.h', False, False),
    ('Reprocess', 'reprocess.h', True, True),
    ]

def describe_class(classname, filename, env=None, description_filename=None, 
                   verbose=False):
    if env is None:
        env = {}
    if description_filename is None:
        description_filename = filename.rsplit('.', 1)[0] + '.py'
    descs = [describe(filename, classname=classname, verbose=verbose)]
    if os.path.isfile(description_filename):
        glbs = globals()
        locs = {}
        execfile(description_filename, glbs, locs)
        if 'desc' not in locs:
            pydesc = {}
        elif callable(locs['desc']):
            pydesc = eval('desc()', glbs, locs)
        else:
            pydesc = locs['desc']
        descs.append(pydesc)
    desc = merge_descriptions(descs)
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
