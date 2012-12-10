#!/usr/bin/env python
 

"""\


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
"""



#############################

#!/usr/bin/env python
 
import os
import sys
import subprocess

import configure

# Thanks to http://patorjk.com/software/taag/  
# and http://www.chris.com/ascii/index.php?art=creatures/dragons
# for ASCII art inspiriation


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

