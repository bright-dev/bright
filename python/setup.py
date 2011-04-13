#!/usr/bin/env python
 
import os
from copy import deepcopy

from distutils.core import setup
from distutils.extension import Extension
from distutils.file_util import copy_file, move_file
from distutils.dir_util import mkpath, remove_tree
from Cython.Distutils import build_ext

import numpy as np

from setup_data import INFO

# Clean up old files
cython_cpp_files = ['src/stlconverters.cpp', 
                    'src/isoname/isoname.cpp', 
                    'src/mass_stream/mass_stream.cpp', 
                    'src/bright/bright.cpp', 
                    ]

for f in cython_cpp_files:
    if os.path.exists(f):
        print "Removing {0}...".format(f)
        os.remove(f)



###########################################
### Set compiler options for extensions ###
###########################################
src_dir = os.path.abspath(os.path.join('src', ''))
cpp_dir = os.path.abspath(os.path.join('..', 'cpp'))
dat_dir = os.path.abspath(os.path.join('..', 'data'))

# Get numpy include dir
numpy_include = np.get_include()

# Path to user's home directory
user_home = os.path.expanduser('~')




#
# For stlconverters
# 
stlconv_ext = {'name': "stlconverters"}

stlconv_ext['sources'] = [os.path.join(src_dir, 'stlconverters.pyx')]
stlconv_ext['include_dirs'] = [src_dir, cpp_dir, numpy_include]
stlconv_ext['language'] = "c++"

if os.name == 'posix':
    stlconv_ext["extra_compile_args"] = ["-Wno-strict-prototypes"]
elif os.name == 'nt':
    stlconv_ext["extra_compile_args"] = ["/EHsc"]
    stlconv_ext["define_macros"] = [("_WIN32", None)]



#
# For isoname
# 
isoname_ext = {'name': 'isoname'}
                 
isoname_ext['sources'] = [
    'bright.cpp', 
    'isoname.cpp', 
    ]
isoname_ext['sources'] = [os.path.join(cpp_dir, s) for s in isoname_ext['sources']] + \
                         [os.path.join(src_dir, 'isoname', 'isoname.pyx')]

isoname_ext['include_dirs'] = [os.path.join(src_dir, 'isoname', ''), 
                               src_dir, cpp_dir, numpy_include]
isoname_ext['language'] = "c++"


if os.name == 'posix':
    isoname_ext["extra_compile_args"] = ["-Wno-strict-prototypes"]
elif os.name == 'nt':
    isoname_ext["extra_compile_args"] = ["/EHsc"]
    isoname_ext["define_macros"] = [("_WIN32", None)]



#
# For MassStream
#
mass_stream_ext = {'name': "mass_stream"}

mass_stream_ext['sources'] = [
    'bright.cpp', 
    'isoname.cpp', 
    'MassStream.cpp',
    ]
mass_stream_ext['sources'] = [os.path.join(cpp_dir, s) for s in mass_stream_ext['sources']] + \
                             [os.path.join(src_dir, 'mass_stream', 'mass_stream.pyx')]
                  
mass_stream_ext['include_dirs'] = [os.path.join(src_dir, 'mass_stream', ''),
                                   os.path.join(src_dir, 'isoname', ''), 
                                   src_dir, cpp_dir, numpy_include]
mass_stream_ext['language'] = "c++"


if os.name == 'posix':
    mass_stream_ext["extra_compile_args"] = ["-Wno-strict-prototypes"]
    mass_stream_ext["libraries"] = [
        "z", 
        "m", 
        "hdf5", 
        "hdf5_hl", 
        "hdf5_cpp", 
        "hdf5_hl_cpp",
        ] 
elif os.name == 'nt':
    mass_stream_ext["libraries"] = [
        "/DEFAULTLIB:szip.lib",
        "/DEFAULTLIB:zlib1.lib",

        # For Dynamic Libs (dll)
        "/DEFAULTLIB:hdf5dll.lib",
        "/DEFAULTLIB:hdf5_hldll.lib",
        "/DEFAULTLIB:hdf5_cppdll.lib",
        "/DEFAULTLIB:hdf5_hl_cppdll.lib",
        ] 

    mass_stream_ext["extra_compile_args"] = ["/EHsc"]

    mass_stream_ext["define_macros"] = [
        ("_WIN32", None),
        ("_HDF5USEDLL_", None),
        ("HDF5CPP_USEDLL", None),
        ]


#
# For bright
#
bright_ext = {'name': "bright"}

bright_ext['sources'] = [
    'bright.cpp', 
    'isoname.cpp', 
    'MassStream.cpp', 
    'FCComp.cpp', 
    'Enrichment.cpp', 
    'Reprocess.cpp', 
    'Storage.cpp', 
    'Reactor1G.cpp', 
    'LightWaterReactor1G.cpp', 
    'FastReactor1G.cpp', 
    'FuelFabrication.cpp', 
    ]
bright_ext['sources'] = [os.path.join(cpp_dir, s) for s in bright_ext['sources']] + \
                        [os.path.join(src_dir, 'bright', 'bright.pyx')]
                  
bright_ext['include_dirs'] = [os.path.join(src_dir, 'bright', ''),
                              os.path.join(src_dir, 'mass_stream', ''),
                              os.path.join(src_dir, 'isoname', ''), 
                              src_dir, cpp_dir, numpy_include]
bright_ext['language'] = "c++"


if os.name == 'posix':
    bright_ext["extra_compile_args"] = ["-Wno-strict-prototypes"]
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
    'decay.h5', 
    'KaeriData.h5', 
    'FR.h5', 
    'LWR.h5',
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
    'bright_data': os.path.join('src', 'bright_data'),
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
    packages = ['bright_data'],
    package_dir = pack_dir,
    package_data = {'bright_data': bright_data_files},
    cmdclass = {'build_ext': build_ext}, 
    ext_modules=[
        Extension(**stlconv_ext), 
        Extension(**isoname_ext), 
        Extension(**mass_stream_ext),
        Extension(**bright_ext), 
        ],
    )

