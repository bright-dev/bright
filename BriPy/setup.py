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
	ext_kwargs["libraries"] = ["boost_python"]
elif os.name == 'nt':
	ext_kwargs["extra_compile_args"] = ["/EHsc"]
        ext_kwargs["define_macros"] = [("_WIN32", None)]

#For MassStream
MassStream_ext_kwargs = deepcopy(ext_kwargs)
if os.name == 'posix':
	MassStream_ext_kwargs["libraries"].extend( ["z", "m", "hdf5_cpp", "hdf5_hl_cpp"] )
elif os.name == 'nt':
        FCComps_ext_kwargs["extra_link_args"] = [
                "/DEFAULTLIB:zlib1.lib",

                #For Dynamic Libs (dll)
                "/DEFAULTLIB:szlibdll.lib",
                "/DEFAULTLIB:hdf5dll.lib",
                "/DEFAULTLIB:hdf5_hldll.lib",
                "/DEFAULTLIB:hdf5_cppdll.lib",
                "/DEFAULTLIB:hdf5_hl_cppdll.lib"
                ] 
        FCComps_ext_kwargs["define_macros"].extend( [("_HDF5USEDLL_", None)] )

#For FCComps
FCComps_ext_kwargs = deepcopy(ext_kwargs)
if os.name == 'posix':
	FCComps_ext_kwargs["libraries"].extend( ["hdf5", "z", "m", "hdf5_hl", "hdf5_cpp", "hdf5_hl_cpp"] )
elif os.name == 'nt':
        FCComps_ext_kwargs["extra_link_args"] = [
                "/DEFAULTLIB:zlib1.lib",

                #For Static Libs (lib)
                #"/DEFAULTLIB:szlib.lib",
                #"/DEFAULTLIB:hdf5.lib",
                #"/DEFAULTLIB:hdf5_hl.lib",
                #"/DEFAULTLIB:hdf5_cpp.lib",
                #"/DEFAULTLIB:hdf5_hl_cpp.lib",

                #For Dynamic Libs (dll)
                "/DEFAULTLIB:szlibdll.lib",
                "/DEFAULTLIB:hdf5dll.lib",
                "/DEFAULTLIB:hdf5_hldll.lib",
                "/DEFAULTLIB:hdf5_cppdll.lib",
                "/DEFAULTLIB:hdf5_hl_cppdll.lib"
                ] 
        FCComps_ext_kwargs["define_macros"].extend( [("_HDF5USEDLL_", None)] )

##########################
### Setup Package Data ###
##########################
pack_data = {'BriPy': ['decay.h5', 'KaeriData.h5', 'FR.h5', 'LWR.h5']}

if os.name == "nt":
        pack_data['isoname']    = ["boost_python-vc90-mt-1_42.dll"]
        pack_data['MassStream'] = ["boost_python-vc90-mt-1_42.dll"]
        pack_data['BriPy'].extend( [
                "boost_python-vc90-mt-1_42.dll",
                "hdf5dll.dll",
                "hdf5_cppdll.dll",
                "hdf5_hldll.dll",
                "hdf5_hl_cppdll.dll",
                "szlibdll.dll",
                "szlib1.dll",
                ])

###################
### Call setup! ###
###################
setup(name="BriPy",
	version = '0.21',
	description = 'Bright/Python',
	author = 'Anthony Scopatz',
	author_email = 'scopatz@gmail.com',
	url = 'http://www.scopatz.com/',
	packages = ['BriPy', 'MassStream', 'isoname'],
	package_dir = {'BriPy': 'src/BriPy', 'MassStream': 'src/MassStream', 'isoname': 'src/isoname'},
	package_data = pack_data,
	#py_modules=["BriPy.__init__"],
	ext_modules=[
		Extension("isoname.isoname",    
			["src/isoname.cpp", "src/bright.cpp", "BriPy_isoname.cpp"], 
			**ext_kwargs),
		Extension("MassStream.MassStream", 
			["src/MassStream.cpp", "src/isoname.cpp", "src/bright.cpp", 
				"BriPy_MassStream.cpp"], 
			**MassStream_ext_kwargs),
		Extension("BriPy.FCComps",    
			["src/bright.cpp", 
				"src/isoname.cpp", 
				"src/MassStream.cpp", 
				"src/FCComp.cpp", 
				"src/Reprocess.cpp", 
				"src/Storage.cpp", 
				"src/Reactor1G.cpp", 
				"src/FastReactor1G.cpp", 
				"src/LightWaterReactor1G.cpp", 
				"BriPy_FCComps.cpp"], 
			**FCComps_ext_kwargs),
		],
	)

