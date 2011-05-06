from enthought.traits.api import HasTraits, Float, Str, Range, Instance, on_trait_change

from enthought.traits.ui.api import View, Item, Tabbed, Group, InstanceEditor

import bright

from utils import _init_comp, _get_comp_pass_number
from mass_stream import _MassStreamView

def Enrichment(IsosIn, enrichment, name):
    """Enrichment:  
    -----------
    This is an enrichment object that enrichs a mass stream.
    """
    _init_comp()
    params = bright.UraniumEnrichmentDefaults()
    params.xP_j = enrichment
    enr = bright.Enrichment(params, name)
    enr.PassNum = _get_comp_pass_number(enr.natural_name)
    enr.doCalc(IsosIn)
    enr.writeout()
    IsosOut = enr.IsosOut
    return IsosOut


class _Enrichment_view(HasTraits):

    IsosIn  = Instance(bright.MassStream)
    IsosOut = Instance(bright.MassStream)

    IsosIn_view  = Instance(_MassStreamView)
    IsosOut_view = Instance(_MassStreamView)

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
        IsosIn = bright.MassStream()
        return IsosIn

    def _IsosIn_view_default(self):
        iiv = _MassStreamView(mass_stream=self.IsosIn)
        return iiv

    def _IsosIn_view_changed(self):
        self.IsosIn = self.IsosIn_view.mass_stream

    # IsosIn functions

    def _IsosOut_default(self):
        IsosOut = bright.MassStream()
        return IsosOut

    def _IsosOut_view_default(self):
        iov = _MassStreamView(mass_stream=self.IsosOut)
        return iov

    def _IsosOut_view_changed(self):
        self.IsosOut = self.IsosOut_view.mass_stream
