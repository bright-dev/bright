from traits.api import HasTraits, Float, Str, Range, Instance, on_trait_change

from traitsui.api import View, Item, Tabbed, Group, InstanceEditor

import material

import bright

from utils import _init_comp, _get_comp_pass_number
from pyne.material import Material
from bright.storage import Storage

class StorageView(HasTraits):

    IsosIn  = Instance(material.Material)
    IsosOut = Instance(material.Material)

    IsosIn_view  = Instance(Material)
    IsosOut_view = Instance(Material)

    name = Str("Storage")

    decay_time = Range(0.0, 10.0**4)

    traits_view = View(
        Tabbed(
            Group(
                Item('name', label="Component Name"), 
                Item('decay_time', label="Storage Time"), 
                label = "Storage"
                ),
             Item('IsosIn_view',
                editor=InstanceEditor(view='traits_view'),
                label = "Isotopic Input",
                show_label=False,
                style='custom',
                resizable=True,
                width = 700,
                height = 600,
                ),
             Item('IsosOut_view',
                editor=InstanceEditor(view='traits_view'),
                label = "Isotopic Output",
                show_label=False,
                style='custom',
                resizable=True,
                width = 700,
                height = 600,
                ),
            ),
        )

    # IsosIn functions

    def _IsosIn_default(self):
        IsosIn = Material({922340: 0.000055, 922350: 0.00720, 922380: 0.992745})
        return IsosIn

    def _IsosIn_view_default(self):
        #iiv = Material(material=self.IsosIn)
        iiv = Material({})
        return iiv

    @on_trait_change('IsosIn_view.material')
    def _IsosIn_view_changed(self):
        self.IsosIn = self.IsosIn_view.material

    # IsosIn functions

    def _IsosOut_default(self):
        IsosOut = Material({922340: 0.000055, 922350: 0.00720, 922380: 0.992745})
        return IsosOut

    def _IsosOut_view_default(self):
        #iov = Material(material=self.IsosOut)
        iov = Material({})
        return iov

    @on_trait_change('IsosOut_view.material')
    def _IsosOut_view_changed(self):
        self.IsosOut = self.IsosOut_view.material


if __name__ == "__main__":
    nu =  Material({922340: 0.000055, 922350: 0.00720, 922380: 0.992745}, 42.0, "Natural Uranium")  

    _storageview = StorageView(IsosIn = nu, IsosOut = nu )
    _storageview.configure_traits()

#    Storage(nu)
