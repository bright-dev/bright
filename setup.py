#!/usr/bin/env python
 
import os
from copy import deepcopy
    
from distutils.core import setup
from distutils.extension import Extension
from distutils.file_util import copy_file, move_file
from distutils.dir_util import mkpath, remove_tree
from Cython.Distutils import build_ext

#from setuptools.command.develop import develop
        
import numpy as np

INFO = {'version': '0.5'}



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



def cpp_ext(name, sources, libs=None, use_hdf5=False):
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
                           ]
    # perfectly general, thanks to dynamic runtime linking of $ORIGIN
    #ext['runtime_library_dirs'] = ['${ORIGIN}/lib', '${ORIGIN}']
    ext['runtime_library_dirs'] = ['${ORIGIN}/lib', '${ORIGIN}', '${ORIGIN}/.']
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

# Pure C/C++ share libraries
# pyne lib
exts.append(cpp_ext("pyne.lib.libbright", ['bright.cpp']))




"""\
#
# For bright
#
bright_ext = {'name': "bright.bright"}

bright_ext['sources'] = [
    'bright.cpp', 
    'h5wrap.cpp',
    'isoname.cpp', 
    'MassStream.cpp', 
    'FCComps.cpp', 
    'FCComp.cpp', 
    'Enrichment.cpp', 
    'FuelFabrication.cpp', 
    'Reprocess.cpp', 
    'Storage.cpp', 
    'FluencePoint.cpp', 
    'ReactorParameters.cpp', 
    'Reactor1G.cpp', 
    'LightWaterReactor1G.cpp', 
    'FastReactor1G.cpp', 
    'ReactorMG.cpp', 
    ]
bright_ext['sources'] = [os.path.join(cpp_dir, s) for s in bright_ext['sources']] + \
                        [os.path.join(src_dir, 'bright', 'bright.pyx')]
                  
bright_ext['include_dirs'] = [os.path.join(src_dir, 'bright', ''),
                              os.path.join(src_dir, 'mass_stream', ''),
                              os.path.join(src_dir, 'isoname', ''), 
                              src_dir, cpp_dir, numpy_include]
bright_ext['language'] = "c++"


if os.name == 'posix':
    #bright_ext["extra_compile_args"] = ["-Wno-strict-prototypes"]
    bright_ext["undef_macros"] = ["NDEBUG"]
    bright_ext["libraries"] = [
        "z", 
        "m", 
        "hdf5", 
        "hdf5_hl", 
        "hdf5_cpp", 
        "hdf5_hl_cpp", 
        ] 

elif os.name == 'nt':
    bright_ext["libraries"] = [
        "/DEFAULTLIB:szip.lib",
        "/DEFAULTLIB:zlib1.lib",

        # For Dynamic Libs (dll)
        "/DEFAULTLIB:hdf5dll.lib",
        "/DEFAULTLIB:hdf5_hldll.lib",
        "/DEFAULTLIB:hdf5_cppdll.lib",
        "/DEFAULTLIB:hdf5_hl_cppdll.lib"
        ] 

    bright_ext["extra_compile_args"] = ["/EHsc"]

    bright_ext["define_macros"] = [
        ("_WIN32", None),
        ("_HDF5USEDLL_", None),
        ("HDF5CPP_USEDLL", None),        
        ]


##########################
### Setup Package Data ###
##########################
bright_data_files = [
    'nuc_data.h5', 
    'decay.h5', 
    'KaeriData.h5', 
    'FR.h5', 
    'LWR.h5',
    'lwr_mg.h5',
    ]
        
pack_dlls_hdf5  = [
    "szip.dll",
    "zlib1.dll",
    "hdf5dll.dll",
    "hdf5_cppdll.dll",
    "hdf5_hldll.dll",
    "hdf5_hl_cppdll.dll",
    ]


pack_dir = {
    'isoname': os.path.join('src', 'isoname'),
    'mass_stream': os.path.join('src', 'mass_stream'), 
    'bright': os.path.join('src', 'bright'), 
#    'bright_data': os.path.join('src', 'bright_data'),
    'bright_data': os.path.join('..', 'data'),
    }
    
pack_data = {'bright': [], 'mass_stream': []}


if os.name == 'posix':
    pack_data['bright'].extend(bright_data_files)
elif os.name == "nt":
    pack_data['mass_stream'].extend(pack_dlls_boost)
    pack_data['mass_stream'].extend(pack_dlls_hdf5)

    pack_data['bright'].extend(bright_data_files)
    pack_data['bright'].extend(pack_dlls_boost)
    pack_data['bright'].extend(pack_dlls_hdf5)




###################
### Call setup! ###
###################
setup(name="bright",
    version = INFO['version'],
    description = 'Bright/Python',
    author = 'Anthony Scopatz',
    author_email = 'scopatz@gmail.com',
    url = 'http://www.scopatz.com/',
    packages = ['bright', 'bright_data', 'bright.gui', 'bright.gui.views', 'bright.gui.views.component_views', 'bright.gui.views.component_views.views', 'bright.gui.models', 'bright.gui.models.class_models'],
    package_dir = pack_dir,
    package_data = {'bright_data': bright_data_files},
    cmdclass = {'build_ext': build_ext}, 
    ext_modules=[
        Extension(**stlconv_ext), 
        Extension(**isoname_ext), 
        Extension(**mass_stream_ext),
        Extension(**bright_ext), 
       ],
    scripts = ['src/scripts/bright_gui']
    )
"""

##########################
### Setup Package Data ###
##########################
packages = ['bright', 'bright.lib', 'bright.gui', 'bright.gui.models', 'bright.gui.models.class_models', 
            'bright.gui.views', 'bright.gui.views.component_views', 'bright.gui.views.component_views.views',]

pack_dir = {'bright': 'bright',}

pack_data = {'bright': ['includes/*.h', 'includes/*.pxd'],
            }

ext_modules=[Extension(**ext) for ext in exts]

# Utility scripts
scripts=['scripts/bright_gui']

###################
### Call setup! ###
###################
if __name__ == "__main__":
    # clean includes dir and recopy files over
    if os.path.exists('bright/includes'):
        remove_tree('bright/includes')
    mkpath('bright/includes')
    for header in (glob.glob('cpp/*.h') + glob.glob('bright/*.pxd')):
        copy_file(header, 'bright/includes')

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

