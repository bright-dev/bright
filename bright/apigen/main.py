"""Top-level automatic API generators for Bright."""
import os
import argparse
from pprint import pprint

from .autodescribe import describe, merge_descriptions

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
    descs = [describe(filename, classname=classname, verbose=verbose), {}]
    if os.path.isfile(description_filename):
        execfile(description_filename, globals(), descs[-1])
    desc = merge_descriptions(descs)
    return desc


def main():
    env = {}
    for classname, fname, mkcython, mkcyclus in CLASSES:
        print("parsing " + classname)
        env[classname] = describe_class(classname, os.path.join('cpp', fname), env=env)
        pprint(env[classname])
        if mkcython:
            pass
        if mkcyclus:
            pass


if __name__ == '__main__':
    main()
