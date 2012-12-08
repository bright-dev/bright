#!/usr/bin/env python
 
import os
import glob
from copy import deepcopy
    
from distutils.core import setup
from distutils.util import get_platform
from distutils.extension import Extension
from distutils.file_util import copy_file, move_file
from distutils.dir_util import mkpath, remove_tree
from distutils.sysconfig import get_python_version
from Cython.Distutils import build_ext

#from setuptools.command.develop import develop
        
import numpy as np

INFO = {'version': '0.6'}



###########################################
### Set compiler options for extensions ###
###########################################
brt_dir = os.path.join('bright')
cpp_dir = os.path.join('cpp')
dat_dir = os.path.join('data')

# Get numpy include dir
numpy_include = np.get_include()

# Get PyNE include and lib dirs
from pyne import pyne_config
pyne_lib = pyne_config.lib
pyne_include = pyne_config.includes

# Path to user's home directory
#user_home = os.path.expanduser('~')

# HDF5 stuff
posix_hdf5_libs = ["z", "m", "hdf5", "hdf5_hl", "hdf5_cpp", "hdf5_hl_cpp",]
nt_hdf5_libs = ["/DEFAULTLIB:szip.lib", "/DEFAULTLIB:zlib1.lib", "/DEFAULTLIB:hdf5dll.lib",
                "/DEFAULTLIB:hdf5_hldll.lib", "/DEFAULTLIB:hdf5_cppdll.lib", "/DEFAULTLIB:hdf5_hl_cppdll.lib", ]
nt_hdf5_extra_compile_args = ["/EHsc"]
nt_hdf5_macros = [("_WIN32", None), ("_HDF5USEDLL_", None), ("HDF5CPP_USEDLL", None), ]



def cpp_ext(name, sources, libs=None, use_hdf5=True):
    """Helper function for setting up extension dictionary.

    Parameters
    ----------
    name : str
        Module name
    sources : list of str
        Files to compile
    libs : list of str
        Additional files to link against
    use_hdf5 : bool
        Link against hdf5?
    """
    ext = {'name': name}

    ext['sources'] = [os.path.join(cpp_dir, s) for s in sources if s.endswith('cpp')] + \
                     [os.path.join(brt_dir, s) for s in sources if s.endswith('pyx')] + \
                     [s for s in sources if not any([s.endswith(suf) for suf in ['cpp', 'pyx']])]

    ext["libraries"] = []
    ext['include_dirs'] = [brt_dir, cpp_dir, pyne_include, numpy_include]
    ext['language'] = "c++"

    # may need to be more general
    ext['library_dirs'] = ['build/lib/bright/lib',
                           'build/lib.{0}-{1}/bright/lib'.format(get_platform(), get_python_version()),
                           pyne_lib,
                           ]
    # perfectly general, thanks to dynamic runtime linking of $ORIGIN
    #ext['runtime_library_dirs'] = ['${ORIGIN}/lib', '${ORIGIN}']
    ext['runtime_library_dirs'] = ['${ORIGIN}/lib', '${ORIGIN}', '${ORIGIN}/.', pyne_lib]
    #ext['runtime_library_dirs'] = ['${ORIGIN}/lib', '${ORIGIN}', '${ORIGIN}/.'] + \
    #                              [os.path.abspath(p) for p in ext['library_dirs']] + \
    #                              [os.path.abspath(p + '/bright/lib') for p in sys.path] + \
    #                              [os.path.abspath(p + '/bright') for p in sys.path] + \
    #                              [os.path.abspath(p) for p in sys.path]


    if os.name == 'posix':
        #ext["extra_compile_args"] = ["-Wno-strict-prototypes"]
        ext["undef_macros"] = ["NDEBUG"]
        if use_hdf5:
            ext["libraries"] += posix_hdf5_libs
        if libs is not None:
            ext["libraries"] += libs
    elif os.name == 'nt':
        ext["extra_compile_args"] = ["/EHsc"]
        ext["define_macros"] = [("_WIN32", None)]

        if use_hdf5:
            ext["libraries"] += nt_hdf5_libs
            ext["extra_compile_args"] += nt_hdf5_extra_compile_args
            ext["define_macros"] += nt_hdf5_macros

        if libs is not None:
            ext["libraries"] += libs

    return ext



#
# For extensions
# 
exts = []
pyne_libs = ['pyne', 'pyne_nucname', 'pyne_data', 'pyne_material']

# Pure C/C++ share libraries
# bright lib
exts.append(cpp_ext("bright.lib.libbright", ['bright.cpp'], pyne_libs))

# fccomp lib
exts.append(cpp_ext("bright.lib.libbright_fccomp", ['fccomp.cpp'], ['bright'] + pyne_libs))

# enrichment lib
exts.append(cpp_ext("bright.lib.libbright_enrichment", ['enrichment.cpp'], ['bright', 'bright_fccomp'] + pyne_libs))

# reprocess lib
exts.append(cpp_ext("bright.lib.libbright_reprocess", ['reprocess.cpp'], ['bright', 'bright_fccomp'] + pyne_libs))

# storage lib
exts.append(cpp_ext("bright.lib.libbright_storage", ['storage.cpp'], ['bright', 'bright_fccomp'] + pyne_libs))

# reactor parameters lib
exts.append(cpp_ext("bright.lib.libbright_reactor_parameters", ['reactor_parameters.cpp'], None, False))

# fluence point lib
exts.append(cpp_ext("bright.lib.libbright_fluence_point", ['fluence_point.cpp'], None, False))

# reactor1g lib
exts.append(cpp_ext("bright.lib.libbright_reactor1g", ['reactor1g.cpp'], 
            ['bright', 'bright_fccomp', 'bright_reactor_parameters', 'bright_fluence_point'] + pyne_libs))

# light_water_reactor1g lib
exts.append(cpp_ext("bright.lib.libbright_light_water_reactor1g", ['light_water_reactor1g.cpp'], 
            ['bright', 'bright_fccomp', 'bright_reactor_parameters', 'bright_fluence_point', 
            'bright_reactor1g'] + pyne_libs))

# fast_reactor1g lib
exts.append(cpp_ext("bright.lib.libbright_fast_reactor1g", ['fast_reactor1g.cpp'], 
            ['bright', 'bright_fccomp', 'bright_reactor_parameters', 'bright_fluence_point', 
            'bright_reactor1g'] + pyne_libs))

# fuel fabrication lib
exts.append(cpp_ext("bright.lib.libbright_fuel_fabrication", ['fuel_fabrication.cpp'], 
            ['bright', 'bright_fccomp', 'bright_reactor_parameters', 'bright_fluence_point', 
            'bright_reactor1g'] + pyne_libs))

# reactormg lib
exts.append(cpp_ext("bright.lib.libbright_reactormg", ['reactormg.cpp'], 
            ['bright', 'bright_fccomp', 'bright_reactor_parameters', 'bright_fluence_point'] + pyne_libs))



# Bright extensions

# bright_config
exts.append(cpp_ext("bright.bright_config", ['bright_config.pyx'], ['bright'] + pyne_libs))

# fccomp
exts.append(cpp_ext("bright.fccomp", ['fccomp.pyx'], ['bright', 'bright_fccomp'] + pyne_libs))

# enrichment
exts.append(cpp_ext("bright.enrichment", ['enrichment.pyx'], 
            ['bright', 'bright_fccomp', 'bright_enrichment'] + pyne_libs))

# reprocess
exts.append(cpp_ext("bright.reprocess", ['reprocess.pyx'], 
            ['bright', 'bright_fccomp', 'bright_reprocess'] + pyne_libs))

# storage
exts.append(cpp_ext("bright.storage", ['storage.pyx'], 
            ['bright', 'bright_fccomp', 'bright_storage'] + pyne_libs))

# reactor parameters
exts.append(cpp_ext("bright.reactor_parameters", ['reactor_parameters.pyx'], 
            ['bright', 'bright_reactor_parameters'], False))

# fluence point
exts.append(cpp_ext("bright.fluence_point", ['fluence_point.pyx'], ['bright', 'bright_fluence_point'], False))

# reactor1g
exts.append(cpp_ext("bright.reactor1g", ['reactor1g.pyx'], 
            ['bright', 'bright_fccomp', 'bright_reactor_parameters', 'bright_fluence_point', 
             'bright_reactor1g'] + pyne_libs))

# light water reactor1g
exts.append(cpp_ext("bright.light_water_reactor1g", ['light_water_reactor1g.pyx'], 
            ['bright', 'bright_fccomp', 'bright_reactor_parameters', 'bright_fluence_point', 
             'bright_reactor1g', 'bright_light_water_reactor1g'] + pyne_libs))

# fast reactor1g
exts.append(cpp_ext("bright.fast_reactor1g", ['fast_reactor1g.pyx'], 
            ['bright', 'bright_fccomp', 'bright_reactor_parameters', 'bright_fluence_point', 
             'bright_reactor1g', 'bright_fast_reactor1g'] + pyne_libs))

# fuel fabrication
exts.append(cpp_ext("bright.fuel_fabrication", ['fuel_fabrication.pyx'], 
            ['bright', 'bright_fccomp', 'bright_reactor_parameters', 'bright_fluence_point', 
             'bright_reactor1g', 'bright_fuel_fabrication'] + pyne_libs))


# reactormg
exts.append(cpp_ext("bright.reactormg", ['reactormg.pyx'], 
            ['bright', 'bright_fccomp', 'bright_reactor_parameters', 'bright_fluence_point', 
             'bright_reactormg'] + pyne_libs))



##########################
### Setup Package Data ###
##########################
packages = ['bright', 'bright.lib', 'bright.gui', 'bright.gui.models', 
            'bright.gui.models.class_models', 'bright.gui.views', 'bright.gui.d3',
            'bright.gui.views.component_views','bright.gui.views.custom_graph_canvas', 
            'bright.gui.views.component_views.views', 'bright.data']

pack_dir = {'bright': 'bright',
            'bright.data': 'data',
            }

pack_data = {'bright': ['includes/*.h', 'includes/bright/*.pxd'],
             'bright.data': ['*.h5'],
            }

ext_modules=[Extension(**ext) for ext in exts]

# Utility scripts
scripts=['scripts/bright_gui']

###################
### Call setup! ###
###################
def main():
    """Perform the Bright setup."""

    # clean includes dir and recopy files over
    if os.path.exists('bright/includes'):
        remove_tree('bright/includes')

    mkpath('bright/includes')
    for header in glob.glob('cpp/*.h') + glob.glob('bright/*.pxd'):
        copy_file(header, 'bright/includes')

    mkpath('bright/includes/bright')
    for header in glob.glob('bright/*.pxd'):
        copy_file(header, 'bright/includes/bright')

    # call setup
    setup(name="bright",
        version = INFO['version'],
        description = 'Bright Nuclear Fuel Cycle Components',
        author = 'Anthony Scopatz',
        author_email = 'scopatz@gmail.com',
        url = 'http://bright-dev.github.com/bright/',
        packages = packages,
        package_dir = pack_dir,
        package_data = pack_data,
        cmdclass = {'build_ext': build_ext},
        ext_modules=ext_modules,
        scripts=scripts,
        )

    # Clean includes after setup has run
    if os.path.exists('bright/includes'):
        remove_tree('bright/includes')

if __name__ == "__main__":
    main()



#############################

#!/usr/bin/env python
 
import os
import sys
import subprocess

import configure

# Thanks to http://patorjk.com/software/taag/  
# and http://www.chris.com/ascii/index.php?art=creatures/dragons
# for ASCII art inspiriation

pyne_logo = """\

                                  /   \       
 _                        )      ((   ))     (                          
(@)                      /|\      ))_((     /|\                          
|-|                     / | \    (/\|/\)   / | \                      (@) 
| | -------------------/--|-voV---\`|'/--Vov-|--\---------------------|-|
|-|                         '^`   (o o)  '^`                          | |
| |                               `\Y/'                               |-|
|-|                                                                   | |
| |        /\             ___           __  __             /\         |-|
|-|       /^~\           / _ \_   _  /\ \ \/__\           /^~\        | |  
| |       /^~\          / /_)/ | | |/  \/ /_\             /^~\        |-|
|-|       /^~\         / ___/| |_| / /\  //__             /^~\        | | 
| |       ^||`         \/     \__, \_\ \/\__/             ^||`        |-|  
|-|        ||                |____/                        ||         | | 
| |       ====                                            ====        |-|
|-|                                                                   | |
| |                                                                   |-|
|-|___________________________________________________________________| |
(@)              l   /\ /         ( (       \ /\   l                `\|-|
                 l /   V           \ \       V   \ l                  (@)
                 l/                _) )_          \I                   
                                   `\ /'
                                     `  
"""

def parse_args():
    distutils = []
    cmake = []
    make = []
    argsets = [distutils, cmake, make]
    i = 0
    for arg in sys.argv:
        if arg == '--':
            i += 1
        else:
            argsets[i].append(arg)
    hdf5opt = [o.split('=')[1] for o in distutils if o.startswith('--hdf5=')]
    if 0 < len(hdf5opt):
        os.environ['HDF5_ROOT'] = hdf5opt[0]  # Expose to CMake
        distutils = [o for o in distutils if not o.startswith('--hdf5=')]
    return distutils, cmake, make


def main_body():
    if not os.path.exists('build'):
        os.mkdir('build')
    sys.argv, cmake_args, make_args = parse_args()
    makefile = os.path.join('build', 'Makefile')
    if not os.path.exists(makefile):
        cmake_cmd = ['cmake', '..'] + cmake_args
        if os.name == 'nt':
            files_on_path = set()
            for p in os.environ['PATH'].split(';')[::-1]:
                if os.path.exists(p):
                    files_on_path.update(os.listdir(p))
            if 'cl.exe' in files_on_path:
                pass
            elif 'sh.exe' in files_on_path:
                cmake_cmd += ['-G "MSYS Makefiles"']
            elif 'gcc.exe' in files_on_path:
                cmake_cmd += ['-G "MinGW Makefiles"']
            cmake_cmd = ' '.join(cmake_cmd)
        rtn = subprocess.check_call(cmake_cmd, cwd='build', shell=(os.name=='nt'))
    rtn = subprocess.check_call(['make'] + make_args, cwd='build')
    cwd = os.getcwd()
    os.chdir('build')
    configure.setup()
    os.chdir(cwd)

def main():
    success = False
    try:
        main_body()
        success = True
    finally:
        configure.final_message(success)

if __name__ == "__main__":
    main()

