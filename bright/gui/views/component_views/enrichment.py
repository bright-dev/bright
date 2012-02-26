from enthought.traits.api import HasTraits, Float, Str, Range, Instance, on_trait_change

from enthought.traits.ui.api import View, Item, Tabbed, Group, InstanceEditor

#import mass_stream

import bright

from utils import _init_comp, _get_comp_pass_number
from material import MaterialView

class EnrichmentView(HasTraits):

    IsosIn  = Instance(mass_stream.MassStream)
    IsosOut = Instance(mass_stream.MassStream)

    IsosIn_view  = Instance(MassStreamView)
    IsosOut_view = Instance(MassStreamView)

    name = Str("Enrichment")

    enrichment = Range(0.0, 1.0)

    traits_view = View(
        Tabbed(
            Group(
                Item('name', label="Component Name"), 
                Item('enrichment', label="Target Enrichment"), 
                label = "Enrichment"
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
        IsosIn = mass_stream.MassStream()
        return IsosIn

    def _IsosIn_view_default(self):
        iiv = MassStreamView(mass_stream=self.IsosIn)
        return iiv

    def _IsosIn_view_changed(self):
        self.IsosIn = self.IsosIn_view.mass_stream

    # IsosIn functions

    def _IsosOut_default(self):
        IsosOut = mass_stream.MassStream()
        return IsosOut

    def _IsosOut_view_default(self):
        iov = MassStreamView(mass_stream=self.IsosOut)
        return iov

    def _IsosOut_view_changed(self):
        self.IsosOut = self.IsosOut_view.mass_stream


if __name__ == "__main__":
    nu = mass_stream.MassStream({922340: 0.000055, 922350: 0.00720, 922380: 0.992745}, 42.0, "Natural Uranium")
    
    _enrichmentview = EnrichmentView(IsosIn = nu, IsoOut = nu)
    _enrichmentview.configure_traits() 
