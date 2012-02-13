from enthought.traits.api import HasTraits, Float, Str, Range, Instance, on_trait_change

from enthought.traits.ui.api import View, Item, Tabbed, Group, InstanceEditor

import material

import bright

from utils import _init_comp, _get_comp_pass_number
from material import MaterialView

class StorageView(HasTraits):

    IsosIn  = Instance(material.Material)
    IsosOut = Instance(material.Material)

    IsosIn_view  = Instance(MaterialView)
    IsosOut_view = Instance(MaterialView)

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
        IsosIn = material.Material({922340: 0.000055, 922350: 0.00720, 922380: 0.992745})
        return IsosIn

    def _IsosIn_view_default(self):
        iiv = MaterialView(material=self.IsosIn)
        return iiv

    @on_trait_change('IsosIn_view.material')
    def _IsosIn_view_changed(self):
        self.IsosIn = self.IsosIn_view.material

    # IsosIn functions

    def _IsosOut_default(self):
        IsosOut = material.Material({922340: 0.000055, 922350: 0.00720, 922380: 0.992745})
        return IsosOut

    def _IsosOut_view_default(self):
        iov = MaterialView(material=self.IsosOut)
        return iov

    @on_trait_change('IsosOut_view.material')
    def _IsosOut_view_changed(self):
        self.IsosOut = self.IsosOut_view.material


if __name__ == "__main__":
    nu = material.Material({922340: 0.000055, 922350: 0.00720, 922380: 0.992745}, 42.0, "Natural Uranium")  

#    _storageview = _StorageView(IsosIn = nu, IsosOut = nu )
#    _storageview.configure_traits()

    Storage(nu)
