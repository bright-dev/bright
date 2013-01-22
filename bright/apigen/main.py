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
    parser.add_argument('--debug', action='store_true', default=False)
    parser.add_argument('-v', '--verbose', action='store_true', 
                                           dest='verbose', default=False)
    ns = parser.parse_args()

    env = {}
    for classname, fname, mkcython, mkcyclus in CLASSES:
        print("parsing " + classname)
        env[classname] = describe_class(classname, os.path.join('cpp', fname), 
                                        env=env, verbose=ns.verbose)
        if ns.verbose:
            pprint(env[classname])
        if mkcython:
            pass
        if mkcyclus:
            pass


if __name__ == '__main__':
    main()
