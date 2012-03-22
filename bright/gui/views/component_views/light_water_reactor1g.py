from traits.api import HasTraits, Float, Str, Range, Instance, on_trait_change, Int, Array, Bool, List, Enum

from traitsui.api import View, Item, Tabbed, Group, InstanceEditor, HGroup, Label, spring, RangeEditor

from chaco.api import OverlayPlotContainer, Plot, ArrayPlotData
from enable.component_editor import ComponentEditor
from chaco.tools.cursor_tool import CursorTool, BaseCursorTool
from chaco.base import n_gon
from numpy import linspace, sin, sqrt, transpose

from pyne.material import Material
import bright
import numpy as np

from utils import _init_comp, _get_comp_pass_number, _get_lwr_data_path
from material import MaterialView

from bright.gui.views.component_views.views.burnup_criticality_estimator import _burnup_criticality_estimator_plot, _k_BU

def _make_pincell_plot(radius, length):
    npoints1 = n_gon(center=(0,0), r=length/(2**0.5), nsides=4, rot_degrees=45)
    npoints2 = n_gon(center=(0,0), r=radius, nsides=100, rot_degrees=45)
    plotdata = ArrayPlotData()
    pincell_plot = Plot(plotdata)
    nxarray1, nyarray1 = transpose(npoints1)
    nxarray2, nyarray2 = transpose(npoints2)
    plotdata.set_data("x" + str(4), nxarray1)
    plotdata.set_data("y" + str(4), nyarray1)
    pincell_plot.plot(("x"+str(4), "y"+str(4)),
                        type="polygon",
                        face_color="lightblue")
                        
    plotdata.set_data("x" + str(100), nxarray2)
    plotdata.set_data("y" + str(100), nyarray2)
    pincell_plot.plot(("x"+str(100), "y"+str(100)),
                        type="polygon",
                        face_color="red")
    pincell_plot.title = "Fuel Pin Cell"
    return pincell_plot

class LightWaterReactor1gView(HasTraits):

    IsosIn  = Instance(Material)
    IsosOut = Instance(Material)

    IsosIn_view  = Instance(MaterialView)
    IsosOut_view = Instance(MaterialView)

    name = Str("Light Water Reactor")

    burnup = Range(10.0, 150.0)
    radius = Range(0.1, 50.0)
    length = Range(0.1, 100.0)
    
    flux = Range(0.0, 10**17)
    fuel_density = Range(0.0, 20.0)
    coolant_density = Range(0.0, 20.0)
    pnl = Range(0.0, 1.1)
    use_disadvantage = Bool()
    hydrogen_rescale = Bool()
    lattice_type = Enum(["Cylindrical", "Planar", "Spherical"])
    open_slots = Float()
    total_slots = Float() 
    
    
    batches = Int(3)

    burnup_plot = Instance(OverlayPlotContainer)
    pincell_plot = Instance(Plot)
    top_cursor = Instance(BaseCursorTool)

    BU_F_ = Array
    k_F_  = Array

    traits_view = View(
        Tabbed(
            Group(
                Group( Item('name', label="Component Name"), 
                       ), 
                #Item('burnup_plot',
                #    editor=ComponentEditor(),
                #    resizable=True, 
                #    springy=True,
                #    show_label=False, 
                #    ),
                HGroup( Item('burnup', label="Burnup", width=0.7),
                       Item('batches', label='Batches', width=0.3),  
                       ), 
                label = "Light Water Reactor"
                ),
            Group(
                Group( Item('radius', label="Fuel Region Radius"),
                       Item('length', label="Pitch"),
                       ),
                Group( Item('pincell_plot',
                            editor=ComponentEditor(),
                            resizable=True,
                            springy=True,
                            show_label=False
                            ),
                        ),
                label = "Fuel Pin Cell"
                ),
            Group(
                Item('open_slots', label="Open Slots"),
                Item('total_slots', label="Total Slots"),
                Item('lattice_type', label="Lattice Type"),
                label = "Lattice"
                ),
            Group(
                Item('flux', label="Flux"),
                Item('fuel_density', label="Fuel Density"),
                Item('coolant_density', label="Coolant Density"),
                Item('pnl', label="Non-leakage Probability"),
                Item('use_disadvantage', label = "Use Disadvantage"),
                Item('hydrogen_rescale', label = "Hydrogen Rescale"),
                label = "Reactor Parameters"
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
        iiv = MaterialView(mass_stream=self.IsosIn)
        return iiv

    @on_trait_change('IsosIn_view.mass_stream')
    def _IsosIn_view_changed(self):
        self.IsosIn = self.IsosIn_view.mass_stream

    # IsosIn functions

    def _IsosOut_default(self):
        IsosOut = Material({922340: 0.000055, 922350: 0.00720, 922380: 0.992745})
        return IsosOut

    def _IsosOut_view_default(self):
        iov = MaterialView(mass_stream=self.IsosOut)
        return iov

    @on_trait_change('IsosOut_view.mass_stream')
    def _IsosOut_view_changed(self):
        self.IsosOut = self.IsosOut_view.mass_stream

    def load_bu_k(self):
        _init_comp()
        lwr_data = _get_lwr_data_path()

        params = bright.LWRDefaults()
        params.BUt = self.burnup
        params.B = self.batches

        lwr = bright.LightWaterReactor1G(lwr_data, params, self.name)
        lwr.IsosIn = self.IsosIn
        lwr.foldMassWeights()

        BU_F_ = np.array(lwr.BU_F_)
        k_F_  =  np.array(lwr.k_F_)

        self.BU_F_ = BU_F_
        self.k_F_  = k_F_

        return 

    #
    # Traits Defaults
    #

    def _burnup_plot_default(self):
        self.load_bu_k()
        bcep = _burnup_criticality_estimator_plot(self.BU_F_, self.k_F_, self.burnup, self.batches)

        line = bcep.plot_components[0]
        top_cursor = line.overlays[self.batches - 1]
        self.top_cursor = top_cursor

        return bcep
    
    def _pincell_plot_default(self):
        pincell_plot = _make_pincell_plot(self.radius, self.length)
        return pincell_plot
    
    def _batches_default(self):
        return 3

    def _burnup_default(self):
        return 40.0
    
    def _radius_default(self):
        return 0.412
    
    def _length_default(self):
        return 1.33
    
    def _flux_default(self):
        return 4.0 * 10.0**14
    
    def _fuel_density_default(self):
        return 10.7
    
    def _coolant_density_default(self):
        return 0.73
    
    def _pnl_default(self):
        return 0.98
    
    def _use_disadvantage_default(self):
        return True
    
    def _lattice_type_default(self):
        return "Cylindrical"
    
    def _hydrogen_rescale_default(self):
        return True
    
    def _open_slots_default(self):
        return 25.0
    
    def _total_slots_default(self):
        return 289.0
        
    #
    # Trait changes
    #

    def _burnup_plot_changed(self):
        line = self.burnup_plot.plot_components[0]
        top_cursor = line.overlays[self.batches - 1]
        self.top_cursor = top_cursor

    def _batches_changed(self):
        self.load_bu_k()
        bcep = _burnup_criticality_estimator_plot(self.BU_F_, self.k_F_, self.burnup, self.batches)
        self.burnup_plot = bcep

    def _burnup_changed(self):
        if self.top_cursor.current_position[0] == self.burnup:
            return

        k = _k_BU(self.BU_F_, self.k_F_, self.burnup)
        self.top_cursor.current_position = self.burnup, k
        
    def _radius_changed(self):
        pincell_plot = _make_pincell_plot(self.radius, self.length)
        self.pincell_plot = pincell_plot
        
    def _length_changed(self):
        pincell_plot = _make_pincell_plot(self.radius, self.length)
        self.pincell_plot = pincell_plot
        
        

    @on_trait_change('top_cursor.current_position')
    def cursor_changed(self):
        if self.top_cursor.current_position[0] == self.burnup:
            return

        self.burnup = self.top_cursor.current_position[0]

if __name__ == "__main__":
    nu = Material({922340: 0.000055, 922350: 0.00720, 922380: 0.992745}, 42.0, "Natural Uranium")

    _lightwaterreactorview = LightWaterReactorView(IsosIn = nu, IsosOut = nu)
    _lightwaterreactorview.configure_traits()
