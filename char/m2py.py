import numpy as np
import re


if_idx_str = ("""if (exist("idx", "var"));\n"""
              """  idx = idx + 1;\n"""
              """else;\n"""
              """  idx = 1;\n"""
              """end;"""
              )

num_pattern = "([0-9]+[.]?[0-9]*[Ee]?[+-]?[0-9]*)"

array_pattern = r"\[[0-9\sEe+-.]*\]"

def num_replace(matchobj):
    return matchobj.group(0) + ','

def convert(filename):
    """Convert a matlab *.m file to a python file."""
    with open(filename, 'r') as mfile:
        f = mfile.read()

    # Keep comments around
    f = f.replace('%', '#')

    # Grab the number of 'if' statements
    IDX = f.count(if_idx_str)

    # Replace if statements with something more meaningful
    fpart = f.partition(if_idx_str)
    f = fpart[0] + "idx = 0" + fpart[2]
    f = f.replace(if_idx_str, 'idx += 1')

    arrays = re.findall(array_pattern, f)
    for a in arrays:
        new_a = re.sub(num_pattern, num_replace, a)
        f = f.replace(a, new_a)

    print f
