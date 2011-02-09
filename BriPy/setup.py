#!/usr/bin/env python
 
import os
from copy import deepcopy

from distutils.core import setup
from distutils.extension import Extension
from distutils.file_util import copy_file, move_file
from distutils.dir_util import mkpath, remove_tree
from Cython.Distutils import build_ext

###########################################
### Set compiler options for extensions ###
###########################################
ext_kwargs = {}
if os.name == 'posix':
    src_dir = 'src/'
    dat_dir = 'src/BriPy/'
    ext_kwargs["libraries"] = ["boost_python"]
elif os.name == 'nt':
    src_dir = '../FCComps/'
    dat_dir = '../data/'
    ext_kwargs["extra_compile_args"] = ["/EHsc"]
    ext_kwargs["define_macros"] = [("_WIN32", None)]

# Path to user's home directory
user_home = os.path.expanduser('~')

# For MassStream
MassStream_ext_kwargs = deepcopy(ext_kwargs)
if os.name == 'posix':
    MassStream_ext_kwargs["libraries"].extend( [
        "z", 
        "m", 
        "hdf5", 
        "hdf5_hl", 
        "hdf5_cpp", 
        "hdf5_hl_cpp",
        ] )
elif os.name == 'nt':
    MassStream_ext_kwargs["extra_link_args"] = [
        "/DEFAULTLIB:szip.lib",
        "/DEFAULTLIB:zlib1.lib",

        # For Dynamic Libs (dll)
#        "/DEFAULTLIB:szlibdll.lib",
        "/DEFAULTLIB:hdf5dll.lib",
        "/DEFAULTLIB:hdf5_hldll.lib",
        "/DEFAULTLIB:hdf5_cppdll.lib",
        "/DEFAULTLIB:hdf5_hl_cppdll.lib",
        ] 
    MassStream_ext_kwargs["define_macros"].extend([
        ("_HDF5USEDLL_", None),
        ("HDF5CPP_USEDLL", None),
        ])

#For FCComps
FCComps_ext_kwargs = deepcopy(ext_kwargs)
if os.name == 'posix':
    FCComps_ext_kwargs["libraries"].extend( [
        "hdf5", 
        "z", 
        "m", 
        "hdf5_hl", 
        "hdf5_cpp", 
        "hdf5_hl_cpp", 
        ] )
elif os.name == 'nt':
    FCComps_ext_kwargs["extra_link_args"] = [
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
    FCComps_ext_kwargs["define_macros"].extend([
        ("_HDF5USEDLL_", None),
        ("HDF5CPP_USEDLL", None),        
        ])

##########################
### Setup Package Data ###
##########################
pack_dir = {
    'BriPy': 'src/BriPy', 
    'MassStream': 'src/MassStream', 
    'isoname': 'src/isoname'
    }
    
pack_data = {'BriPy': []}

BriPy_data_files = [
    'decay.h5', 
    'KaeriData.h5', 
    'FR.h5', 
    'LWR.h5',
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
    pack_data['BriPy'].extend(BriPy_data_files)
elif os.name == "nt":
    pack_data['isoname'] = []
    pack_data['isoname'].extend(pack_dlls_boost)

    pack_data['MassStream'] = []
    pack_data['MassStream'].extend(pack_dlls_boost)
    pack_data['MassStream'].extend(pack_dlls_hdf5)

    pack_data['BriPy'].extend(BriPy_data_files)
    pack_data['BriPy'].extend(pack_dlls_boost)
    pack_data['BriPy'].extend(pack_dlls_hdf5)

    # Copy over actual data files, instead of symlinks
    cp_symlinks = True
    if 'build' in os.listdir('.'):
        if 'temp' in os.listdir('build/'):
            cp_symlinks = False

    if cp_symlinks:
        mkpath('build/temp/')
        for f in BriPy_data_files:
            copy_file(pack_dir['BriPy'] + '/' + f, 'build/temp/' + f, verbose=True)

        for f in BriPy_data_files:
            copy_file(dat_dir + f, pack_dir['BriPy'] + '/' + f, verbose=True)

###################
### Call setup! ###
###################
setup(name="BriPy",
    version = '0.23',
    description = 'Bright/Python',
    author = 'Anthony Scopatz',
    author_email = 'scopatz@gmail.com',
    url = 'http://www.scopatz.com/',
#    packages = ['BriPy', 'MassStream', 'isoname'],
#    package_dir = pack_dir,
#    package_data = pack_data,
    cmdclass = {'build_ext': build_ext}, 
    ext_modules=[
        Extension("isoname", 
                  ['src/bright.cpp', 
                   'src/isoname.cpp', 
                   'src/isoname/isoname_wrapper.pyx',
                   ],
                  include_dirs=['src/isoname/'],
                  libraries=[],
                  language="c++",
                  ), 
        Extension("MassStream", 
                  ['src/bright.cpp', 
                   'src/isoname.cpp', 
                   'src/MassStream.cpp', 
                   'src/MassStream/mass_stream_wrapper.pyx',
                   ],
                  include_dirs=['src/MassStream/'],
                  libraries=MassStream_ext_kwargs["libraries"], 
                  language="c++",
                  ), 
        Extension("BriPy", 
                  ['src/bright.cpp', 
                   'src/isoname.cpp', 
                   'src/MassStream.cpp', 
                   'src/FCComp.cpp', 
                   'src/BriPy/bright.pyx',
                   ],
                  include_dirs=['src/BriPy/', 'src/MassStream/'],
                  libraries=FCComps_ext_kwargs["libraries"], 
                  language="c++",
                  ), 
#        Extension("BriPy.FCComps", [
#            src_dir + "bright.cpp", 
#            src_dir + "isoname.cpp", 
#            src_dir + "MassStream.cpp", 
#            src_dir + "FCComp.cpp", 
#            src_dir + "Reprocess.cpp", 
#            src_dir + "Storage.cpp", 
#            src_dir + "Enrichment.cpp", 
#            src_dir + "Reactor1G.cpp", 
#            src_dir + "LightWaterReactor1G.cpp", 
#            src_dir + "FastReactor1G.cpp", 
#            src_dir + "FuelFabrication.cpp", 
#            "BriPy_FCComps.cpp"], **FCComps_ext_kwargs),
        ],
    )

if os.name == 'posix':
    pass
elif os.name == "nt":
    print "Cleaning Windows specific files."

    # Copy symlinks over data files. 
    # Hopefully, leaving the repository in the previous state.
    for f in BriPy_data_files:
        copy_file('build/temp/' + f, pack_dir['BriPy'] + '/' + f, verbose=True)

    remove_tree('build/temp/', verbose=True)
