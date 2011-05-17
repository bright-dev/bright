from enthought.traits.api import HasTraits, Float, Str, Range, Instance, on_trait_change

from enthought.traits.ui.api import View, Item, Tabbed, Group, InstanceEditor

import mass_stream

import bright

from utils import _init_comp, _get_comp_pass_number
from mass_stream_views import _MassStreamView

def Storage(IsosIn, decay_time, name):
    """Storage:  
    --------
    This is a storage object that runs a decay calculation on a mass stream.
    """
    _init_comp()
    stor = bright.Storage(name)
    stor.PassNum = _get_comp_pass_number(stor.natural_name)
    stor.doCalc(IsosIn, decay_time)
    stor.writeout()
    IsosOut = stor.IsosOut
    return IsosOut


class _StorageView(HasTraits):

    IsosIn  = Instance(mass_stream.MassStream)
    IsosOut = Instance(mass_stream.MassStream)

    IsosIn_view  = Instance(_MassStreamView)
    IsosOut_view = Instance(_MassStreamView)

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
        IsosIn = bright.MassStream()
        return IsosIn

    def _IsosIn_view_default(self):
        iiv = _MassStreamView(mass_stream=self.IsosIn)
        return iiv

    @on_trait_change('IsosIn_view.mass_stream')
    def _IsosIn_view_changed(self):
        self.IsosIn = self.IsosIn_view.mass_stream

    # IsosIn functions

    def _IsosOut_default(self):
        IsosOut = bright.MassStream()
        return IsosOut

    def _IsosOut_view_default(self):
        iov = _MassStreamView(mass_stream=self.IsosOut)
        return iov

    @on_trait_change('IsosOut_view.mass_stream')
    def _IsosOut_view_changed(self):
        self.IsosOut = self.IsosOut_view.mass_stream


if __name__ == "__main__":
    nu = mass_stream.MassStream({922340: 0.000055, 922350: 0.00720, 922380: 0.992745}, 42.0, "Natural Uranium")  

    _storageview = _StorageView(IsosIn = nu, IsosOut = nu )
    _storageview.configure_traits()
