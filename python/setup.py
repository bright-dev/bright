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

# Get numpy include dir
numpy_include = np.get_include()

###########################################
### Set compiler options for extensions ###
###########################################
ext_kwargs = {}

src_dir = os.path.join('..', '/cpp/')
dat_dir = os.path.join('..', '/data/')

if os.name == 'posix':
#    src_dir = 'src/'
#    dat_dir = 'src/bright/'
    ext_kwargs["libraries"] = ["boost_python"]
elif os.name == 'nt':
    ext_kwargs["extra_compile_args"] = ["/EHsc"]
    ext_kwargs["define_macros"] = [("_WIN32", None)]

# Path to user's home directory
user_home = os.path.expanduser('~')

# For MassStream
mass_stream_ext_kwargs = deepcopy(ext_kwargs)

if os.name == 'posix':
    mass_stream_ext_kwargs["libraries"].extend( [
        "z", 
        "m", 
        "hdf5", 
        "hdf5_hl", 
        "hdf5_cpp", 
        "hdf5_hl_cpp",
        ] )
elif os.name == 'nt':
    mass_stream_ext_kwargs["extra_link_args"] = [
        "/DEFAULTLIB:szip.lib",
        "/DEFAULTLIB:zlib1.lib",

        # For Dynamic Libs (dll)
#        "/DEFAULTLIB:szlibdll.lib",
        "/DEFAULTLIB:hdf5dll.lib",
        "/DEFAULTLIB:hdf5_hldll.lib",
        "/DEFAULTLIB:hdf5_cppdll.lib",
        "/DEFAULTLIB:hdf5_hl_cppdll.lib",
        ] 

    mass_stream_ext_kwargs["define_macros"].extend([
        ("_HDF5USEDLL_", None),
        ("HDF5CPP_USEDLL", None),
        ])

# For bright
bright_ext_kwargs = deepcopy(ext_kwargs)

if os.name == 'posix':
    bright_ext_kwargs["libraries"].extend( [
        "hdf5", 
        "z", 
        "m", 
        "hdf5_hl", 
        "hdf5_cpp", 
        "hdf5_hl_cpp", 
        ] )
elif os.name == 'nt':
    bright_ext_kwargs["extra_link_args"] = [
        "/DEFAULTLIB:szip.lib",
        "/DEFAULTLIB:zlib1.lib",

        # For Static Libs (lib)
        #"/DEFAULTLIB:szlib.lib",
        #"/DEFAULTLIB:hdf5.lib",
        #"/DEFAULTLIB:hdf5_hl.lib",
        #"/DEFAULTLIB:hdf5_cpp.lib",
        #"/DEFAULTLIB:hdf5_hl_cpp.lib",

        # For Dynamic Libs (dll)
        #"/DEFAULTLIB:szlibdll.lib",
        "/DEFAULTLIB:hdf5dll.lib",
        "/DEFAULTLIB:hdf5_hldll.lib",
        "/DEFAULTLIB:hdf5_cppdll.lib",
        "/DEFAULTLIB:hdf5_hl_cppdll.lib"
        ] 
    bright_ext_kwargs["define_macros"].extend([
        ("_HDF5USEDLL_", None),
        ("HDF5CPP_USEDLL", None),        
        ])

##########################
### Setup Package Data ###
##########################
pack_dir = {
    'isoname': 'src/isoname',
    'mass_stream': 'src/mass_stream', 
    'bright': 'src/bright', 
    'bright_data': 'src/bright_data',
    }
    
pack_data = {'bright': []}

bright_data_files = [
    os.path.join(dat_dir, 'decay.h5'), 
    os.path.join(dat_dir, 'KaeriData.h5'), 
    os.path.join(dat_dir, 'FR.h5'), 
    os.path.join(dat_dir, 'LWR.h5'),
    ]
        
pack_dlls_boost= ["boost_python-vc90-mt-1_44.dll"]

pack_dlls_hdf5  = [
    "szip.dll",
    "zlib1.dll",
    "hdf5dll.dll",
    "hdf5_cppdll.dll",
    "hdf5_hldll.dll",
    "hdf5_hl_cppdll.dll",
    ]

if os.name == 'posix':
    pack_data['bright'].extend(bright_data_files)
elif os.name == "nt":
    pack_data['isoname'] = []
    pack_data['isoname'].extend(pack_dlls_boost)

    pack_data['mass_stream'] = []
    pack_data['mass_stream'].extend(pack_dlls_boost)
    pack_data['mass_stream'].extend(pack_dlls_hdf5)

    pack_data['bright'].extend(bright_data_files)
    pack_data['bright'].extend(pack_dlls_boost)
    pack_data['bright'].extend(pack_dlls_hdf5)

    # Copy over actual data files, instead of symlinks
    cp_symlinks = True
    if 'build' in os.listdir('.'):
        if 'temp' in os.listdir('build/'):
            cp_symlinks = False

    if cp_symlinks:
        mkpath('build/temp/')
        for f in bright_data_files:
            copy_file(pack_dir['bright'] + '/' + f, 'build/temp/' + f, verbose=True)

        for f in bright_data_files:
            copy_file(dat_dir + f, pack_dir['bright'] + '/' + f, verbose=True)

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
        Extension("stlconverters", 
                  ['src/stlconverters.pyx',
                   ],
                  include_dirs=['src/', numpy_include],
                  libraries=[],
                  language="c++",
                  ), 
        Extension("isoname", 
                  ['src/bright.cpp', 
                   'src/isoname.cpp', 
                   'src/isoname/isoname.pyx',
                   ],
                  include_dirs=['src/isoname/', 'src/', numpy_include],
                  libraries=[],
                  language="c++",
                  ), 
        Extension("mass_stream", 
                  ['src/bright.cpp', 
                   'src/isoname.cpp', 
                   'src/MassStream.cpp',
                   'src/mass_stream/mass_stream.pyx',
                   ],
                  include_dirs=['src/mass_stream/', 'src/isoname/', 'src/', numpy_include],
                  libraries=mass_stream_ext_kwargs["libraries"], 
                  language="c++",
                  ), 
        Extension("bright", 
                  ['src/bright.cpp', 
                   'src/isoname.cpp', 
                   'src/MassStream.cpp', 
                   'src/FCComp.cpp', 
                   'src/Enrichment.cpp', 
                   'src/Reprocess.cpp', 
                   'src/Storage.cpp', 
                   'src/Reactor1G.cpp', 
                   'src/LightWaterReactor1G.cpp', 
                   'src/FastReactor1G.cpp', 
                   'src/FuelFabrication.cpp', 
                   'src/bright/bright.pyx',
                   ],
                  include_dirs=['src/bright/', 'src/mass_stream/', 'src/', numpy_include],
                  libraries=bright_ext_kwargs["libraries"], 
                  language="c++",
                  ), 
        ],
    )

if os.name == 'posix':
    pass
elif os.name == "nt":
    print "Cleaning Windows specific files."

    # Copy symlinks over data files. 
    # Hopefully, leaving the repository in the previous state.
    for f in bright_data_files:
        copy_file('build/temp/' + f, pack_dir['bright'] + '/' + f, verbose=True)

    remove_tree('build/temp/', verbose=True)
