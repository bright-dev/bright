#!/usr/bin/env python
 
from distutils.core import setup
from distutils.extension import Extension
import subprocess
import os
from copy import deepcopy

###########################################
### Set compiler options for extensions ###
###########################################
ext_kwargs = {}
if os.name == 'posix':
    src_dir = 'src/'
    ext_kwargs["libraries"] = ["boost_python"]
elif os.name == 'nt':
    src_dir = '../FCComps/'
    ext_kwargs["extra_compile_args"] = ["/EHsc"]
    ext_kwargs["define_macros"] = [("_WIN32", None)]

# For MassStream
MassStream_ext_kwargs = deepcopy(ext_kwargs)
if os.name == 'posix':
    MassStream_ext_kwargs["libraries"].extend( [
        "z", 
        "m", 
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
pack_data = {'BriPy': ['decay.h5', 'KaeriData.h5', 'FR.h5', 'LWR.h5']}

pack_dlls_boost= ["boost_python-vc90-mt-1_44.dll"]
pack_dlls_hdf5  = [
    "szip.dll",
    "zlib1.dll",
    "hdf5dll.dll",
    "hdf5_cppdll.dll",
    "hdf5_hldll.dll",
    "hdf5_hl_cppdll.dll",
    ]

if os.name == "nt":
    pack_data['isoname'] = []
    pack_data['isoname'].extend(pack_dlls_boost)

    pack_data['MassStream'] = []
    pack_data['MassStream'].extend(pack_dlls_boost)
    pack_data['MassStream'].extend(pack_dlls_hdf5)

    pack_data['BriPy'].extend(pack_dlls_boost)
    pack_data['BriPy'].extend(pack_dlls_hdf5)

###################
### Call setup! ###
###################
setup(name="BriPy",
    version = '0.23',
    description = 'Bright/Python',
    author = 'Anthony Scopatz',
    author_email = 'scopatz@gmail.com',
    url = 'http://www.scopatz.com/',
    packages = ['BriPy', 'MassStream', 'isoname'],
    package_dir = {'BriPy': 'src/BriPy', 'MassStream': 'src/MassStream', 'isoname': 'src/isoname'},
    package_data = pack_data,
    #py_modules=["BriPy.__init__"],
    ext_modules=[
        Extension("isoname.isoname", [
            src_dir + "isoname.cpp", 
            src_dir + "bright.cpp", 
            "BriPy_isoname.cpp"], **ext_kwargs),
        Extension("MassStream.MassStream", [ 
            src_dir + "MassStream.cpp", 
            src_dir + "isoname.cpp", 
            src_dir + "bright.cpp", 
            "BriPy_MassStream.cpp"], **MassStream_ext_kwargs),
        Extension("BriPy.FCComps", [
            src_dir + "bright.cpp", 
            src_dir + "isoname.cpp", 
            src_dir + "MassStream.cpp", 
            src_dir + "FCComp.cpp", 
            src_dir + "Reprocess.cpp", 
            src_dir + "Storage.cpp", 
            src_dir + "Enrichment.cpp", 
            src_dir + "Reactor1G.cpp", 
            src_dir + "LightWaterReactor1G.cpp", 
            src_dir + "FastReactor1G.cpp", 
            "BriPy_FCComps.cpp"], **FCComps_ext_kwargs),
        ],
    )
