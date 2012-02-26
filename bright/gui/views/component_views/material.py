from __future__ import absolute_import

import numpy as np
from enthought.traits.api import HasTraits, Instance, Str, Float, File, List, on_trait_change, Range
from enthought.traits.ui.api import View, Item, Group, VGroup, HGroup, TableEditor, InstanceEditor, Tabbed

from enthought.traits.ui.table_column import ObjectColumn
from enthought.traits.ui.table_filter import EvalFilterTemplate, MenuFilterTemplate, RuleFilterTemplate, EvalTableFilter, TableFilter

import bright
from pyne import nucname
from pyne.material import Material

###############################################################################
# General Material View Helpers
###############################################################################


class _IsoEntry(HasTraits):
    """Isotope entry (row) for the mass stream table."""

    isotope = Str

    #mass_weight = Range(0.0, np.inf, 0.0)
    mass_weight = Float(0.0)

# Mass Stream isotopic filters
UFilter   = TableFilter(allowed=lambda x: isoname.mixed_2_zzaaam(str(x.isotope))/10000 == 92,        name='Uranium')
PUFilter  = TableFilter(allowed=lambda x: isoname.mixed_2_zzaaam(str(x.isotope))/10000 == 94,        name='Plutonium')
LANFilter = TableFilter(allowed=lambda x: isoname.mixed_2_zzaaam(str(x.isotope))/10000 in isoname.lan, name='Lanthanides')
ACTFilter = TableFilter(allowed=lambda x: isoname.mixed_2_zzaaam(str(x.isotope))/10000 in isoname.act, name='Actinides')
TRUFilter = TableFilter(allowed=lambda x: isoname.mixed_2_zzaaam(str(x.isotope))/10000 in isoname.tru, name='Transuranics')
MAFilter  = TableFilter(allowed=lambda x: isoname.mixed_2_zzaaam(str(x.isotope))/10000 in isoname.ma,  name='Minor Actinides')
FPFilter  = TableFilter(allowed=lambda x: isoname.mixed_2_zzaaam(str(x.isotope))/10000 in isoname.fp,  name='Fission Products')


material_editor = TableEditor(
    columns = [ ObjectColumn(name='isotope',
                    label='Isotope',
                    #tooltip='Isotope name in natrual (LLAAAM) form.',
#                    width=0.20,
                    ),
                ObjectColumn(name='mass_weight',
                    label='Mass Weight', 
                    #tooltip='Mass weight of the corresponding isotope.',
#                    width=0.20, 
                    ),
                ],
    deletable   = True,   
    sort_model  = True,   
    auto_size   = False,
    orientation = 'vertical',
    filters     = [ UFilter, 
                    PUFilter, 
                    LANFilter,
                    ACTFilter,
                    TRUFilter,
                    MAFilter,
                    FPFilter,
                    EvalFilterTemplate, 
                    MenuFilterTemplate, 
                    RuleFilterTemplate 
                    ],
    search      = EvalTableFilter(expression="isotope == 'U235'"),
    row_factory = _IsoEntry, 
    row_factory_kw = {'isotope': "New Isotope", "mass_weight": 0.0}, 
#    auto_add = True,
    sortable = True, 
    configurable = False,
    editable = True,
    show_toolbar = True,
    reorderable = True,
    selected = 'table_selection',
    edit_on_first_click = False,
    )



class _Material_View(HasTraits):
    """At last! A Material view."""

    material = Instance(Material)

    name = Str('')

    mass = Float

    file = File(filter=["Text file (*.txt)|*.txt", "HDF5 (*.h5)|*.h5", "All files|*"], auto_set=False)

    iso_entries = List(_IsoEntry)

    traits_view = View(
        VGroup(
            Item('file', label='Load From File'),
            HGroup(
                Item('name', label='Mass Stream Name', width=0.5),
                Item('mass', label='Total Mass', width=0.5),
                #show_border=True,
                #padding=10,
                ),
            Item('iso_entries',
                show_label  = False,
                editor      = material_editor,
                ),
            ),
        width     = 0.4, 
        height    = 0.4,
        resizable = True,
        )


    def get_iso_entries_from_comp(self, comp=None):
        if comp == None:
            comp = self.material.mult_by_mass()

        iso_entries = [ _IsoEntry(isotope=isoname.zzaaam_2_LLAAAM(iso).capitalize(), mass_weight=comp[iso]) 
                        for iso in sorted(comp.keys()) ]
        return iso_entries

    def get_comp_from_iso_entries(self, mass=None):
        comp = {}

        if mass == None:
            mass = self.mass

        if mass == 0.0:
            return comp

        for iso in self.iso_entries:
            try:
                isozz = isoname.mixed_2_zzaaam(str(iso.isotope))
            except RuntimeError:
                continue

            comp[isozz] = float(iso.mass_weight / mass)

        return comp 

    def get_mass_from_iso_entries(self):
        mass = 0.0
        for iso in self.iso_entries:
            mass += iso.mass_weight
        return mass

    @on_trait_change('iso_entries.isotope')
    def update_material(self):
        comp = self.get_comp_from_iso_entries()
        ms = Material(comp, self.mass, str(self.name))
        self.material = ms

    def update_from_material(self, ms=None):
        if ms == None:
            ms = self.material

        self.name = ms.name
        self.mass = ms.mass
        iso_entries = self.get_iso_entries_from_comp(ms.multByMass())
        self.iso_entries = iso_entries

    #
    # Default trait values 
    #
    
    def _iso_entries_default(self):
        iso_entries = self.get_iso_entries_from_comp()
        return iso_entries

    def _name_default(self):
        return self.material.name

    def _mass_default(self):
        return self.material.mass

    #
    # Changed trait methods
    #

    def _name_changed(self):
        self.update_material()

    def _mass_changed(self, old, new):
        # Don't zero out isotopics 
        if (new == 0.0):
            return 

        # Recover mass from isotopics from the entries
        old = self.get_mass_from_iso_entries()
        if (old == 0.0):
            return 

        # Update entries if they weren't edited
        mass_ratio = new / old
        if mass_ratio != 1.0:
            for iso in self.iso_entries:
                iso.mass_weight = iso.mass_weight * mass_ratio

        # Update mass stream
        self.update_material()

    @on_trait_change('iso_entries.mass_weight')
    def update_from_mass_weights(self):
        mass = self.get_mass_from_iso_entries()
        self.mass = mass


    def _file_changed(self):
        ms = bright.Material()

        rpart = self.file.rpartition('.')

        # HDF5 does not work!
        # Should launch a view to select group and row here
        row = None
        groupname = ''
        if (rpart[2] in ['h5', 'hdf5', 'H5', 'HDF5']) and (groupname != ''):
            if isinstance(row, int):
                ms.load_from_hdf5(str(self.file), groupname, row)
            else:
                ms_out.load_from_hdf5(str(self.file), groupname)

        elif (rpart[2] in ['txt', 'TXT']):
            ms.load_from_text(str(self.file))

        self.update_from_material(ms)
        self.material = ms

    def _iso_entries_items_changed(self):
        self.update_material()


###############################################################################
# Mass Stream Block
###############################################################################

def Material(ms_in):
    """Mass Stream
    ***********
    This provides a simple, standard way to represent fuel cycle mass flows.  The ``Material`` objects
    represent the flows (arrows) between fuel cycle component objects.  A mass stream has three main attributes:

    * ``Material.mass``: The current mass of the flow in the units of the problem.
    * ``Material.comp``: A dictionary (map) of isotopes to their normalized fractional weights.  Isotopic keys are given in `zzaaam` (integer) form.
    * ``Material.name``: (Optional) The label or name  of the stream.

    """

    if isinstance(ms_in, bright.Material):
        ms_out = ms_in

    return ms_out


class MaterialView(HasTraits):

    ms_out = Instance(Material)

    ms_out_view = Instance(_Material_View)

    traits_view = View(
        Item('ms_out_view',
            editor=InstanceEditor(view='traits_view'),
            label = "Mass Stream",
            show_label=False,
            style='custom',
            resizable=True,
            width = 700,
            height = 600,
            ),
        )

    def _ms_out_default(self):
        ms_out = pyne.Material()
        return ms_out

    def _ms_out_view_default(self):
        _msov = MaterialView(material=self.ms_out)
        return _msov

    @on_trait_change('ms_out_view.material')
    def _ms_out_view_changed(self):
        self.ms_out = self.ms_out_view.material


###############################################################################
# Natural Uranium Block
###############################################################################

def NaturalUranium(nu_in):
    """Natural Uranium
    ***************
    A natural uranium mass stream.

    """
    if isinstance(nu_in, Material):
        nu_out = nu_in
    else:
        nu_out = Material({922340: 0.000055, 922350: 0.00720, 922380: 0.992745}, 1.0, "Natural Uranium")
    return nu_out


class _NaturalUranium_view(HasTraits):

    nu_in = nu_out = Instance(Material)

    nu_view = Instance(MaterialView)

    traits_view = View(
        Item('nu_view',
            editor=InstanceEditor(view='traits_view'),
            label = "Natural Uranium",
            show_label=False,
            style='custom',
            resizable=True,
            width = 700,
            height = 600,
            ),
        )

    def _nu_in_default(self):
        nu = NaturalUranium(None)
        return nu

    def _nu_view_default(self):
        _msov = MaterialView(material=self.nu_in)
        return _msov

    @on_trait_change('nu_view.material')
    def _nu_view_changed(self):
        self.nu_in = self.nu_out = self.nu_view.material


###############################################################################
# Depleted Uranium Block
###############################################################################

def DepletedUranium(du_in):
    """Depleted Uranium
    ****************
    A depleted uranium mass stream.

    """
    if isinstance(du_in, Material):
        du_out = du_in
    else:
        du_out = Material({922350: 0.0025, 922380: 0.9975}, 1.0, "Depleted Uranium")
    return du_out


class _DepletedUranium_view(HasTraits):

    du_in = du_out = Instance(Material)

    du_view = Instance(MaterialView)

    traits_view = View(
        Item('du_view',
            editor=InstanceEditor(view='traits_view'),
            label = "Depleted Uranium",
            show_label=False,
            style='custom',
            resizable=True,
            width = 700,
            height = 600,
            ),
        )

    def _du_in_default(self):
        du = DepletedUranium(None)
        return du

    def _du_view_default(self):
        _msov = MaterialView(material=self.du_in)
        return _msov

    @on_trait_change('du_view.material')
    def _du_view_changed(self):
        self.du_in = self.du_out = self.du_view.material


###############################################################################
# Low Enriched Uranium Block
###############################################################################

def LowEnrichedUranium(leu_in):
    """Low Enriched Uranium
    ********************
    A low enriched uranium mass stream.

    """
    if isinstance(leu_in, Material):
        leu_out = leu_in
    else:
        leu_out = Material({922350: 0.05, 922380: 0.95}, 1.0, "Low Enriched Uranium")
    return leu_out


class _LowEnrichedUranium_view(HasTraits):

    leu_in = leu_out = Instance(Material)

    leu_view = Instance(MaterialView)

    traits_view = View(
        Item('leu_view',
            editor=InstanceEditor(view='traits_view'),
            label = "Low Enriched Uranium",
            show_label=False,
            style='custom',
            resizable=True,
            width = 700,
            height = 600,
            ),
        )

    def _leu_in_default(self):
        leu = LowEnrichedUranium(None)
        return leu

    def _leu_view_default(self):
        _msov = MaterialView(material=self.leu_in)
        return _msov

    @on_trait_change('leu_view.material')
    def _du_view_changed(self):
        self.leu_in = self.leu_out = self.leu_view.material



# A sample of how the view is suppossed to work
if __name__ == "__main__":
    nu = Material({922340: 0.000055, 922350: 0.00720, 922380: 0.992745}, 42.0, "Natural Uranium")

    _msview = MaterialView(material=nu)
    _msview.configure_traits()
#    _msview.edit_traits()
