import os

# prepare temporary filename
if not os.path.isdir('build'):
    os.mkdir('build')
tempfile = os.path.join('build', 'xdress_temp')

package = 'bright'
sourcedir = 'cpp'
packagedir = 'bright'

extra_types = 'pyne.extra_types'  # non-default value

stlcontainers = []

stlcontainers_module = 'pyne.stlcontainers'

classes = [
    # classname, source filename[, bindings filename]
    ('FCComp', 'fccomp'),
    ('EnrichmentParameters', 'enrichment_parameters'),
    ('Enrichment', 'bright_enrichment', 'enrichment'),
    ('Reprocess', 'reprocess'),
    ('decay_nuc', 'storage', tempfile),
    ('from_nuc_struct', 'storage', tempfile),
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
