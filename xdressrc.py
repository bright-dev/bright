import pyne
import xdress.typesystem as ts

package = 'bright'
sourcedir = 'cpp'
packagedir = 'bright'

includes = [pyne.includes]

extra_types = 'pyne.extra_types'  # non-default value
make_extra_types = False

stlcontainers = []
make_stlcontainers = False
stlcontainers_module = 'pyne.stlcontainers'

classes = [
    # classname, source filename[, bindings filename]
    ('FCComp', 'fccomp'),
    ('EnrichmentParameters', 'enrichment_parameters'),
    ('Enrichment', 'bright_enrichment', 'enrichment'),
    ('Reprocess', 'reprocess'),
    ('decay_nuc', 'storage'),
    ('from_nuc_struct', 'storage', None),
    ('Storage', 'storage'),
    ('FluencePoint', 'fluence_point'),
    ('ReactorParameters', 'reactor_parameters'),
    ('Reactor1G', 'reactor1g'),
    ('LightWaterReactor1G', 'light_water_reactor1g'),
    ('FastReactor1G', 'fast_reactor1g'),
    ('FuelFabrication', 'fuel_fabrication'),
    ('ReactorMG', 'reactormg'),
    ]

functions = []


# hack in some material registrations
ts.register_class('Material', 
    cython_c_type='cpp_material.Material', cython_cimport=('pyne', 'cpp_material'),
    cython_cy_type='material._Material', cython_py_type='material.Material', 
    cython_template_class_name='Material', cython_cyimport=('pyne', 'material'),
    cython_pyimport=('pyne', 'material'), 
    cython_c2py=('{pytype}({var})',
                 ('{proxy_name} = {pytype}()\n'
                  '{proxy_name}.mat_pointer[0] = {var}'),
                 ('if {cache_name} is None:\n'
                  '    {proxy_name} = {pytype}(free_mat=False)\n'
                  '    {proxy_name}.mat_pointer = &{var}\n'
                  '    {cache_name} = {proxy_name}\n')),
    cython_py2c=(
      '{proxy_name} = {pytype}({var}, free_mat=not isinstance({var}, {cytype}))',
      '{proxy_name}.mat_pointer[0]'),
    )

ts.register_specialization(('map', 'str', ('Material', '*'), 0),
    cython_c_type='material._MapStrMaterial',
    cython_cy_type='material._MapStrMaterial',
    cython_py_type='material.MapStrMaterial',
    cython_cimport=(('pyne', 'material'),),
    cython_cyimport=(('pyne', 'material'),),
    cython_pyimport=(('pyne', 'material'),),
    )
