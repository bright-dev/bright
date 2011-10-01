# Major library imports
import numpy as np
from numpy import arange, sort, compress, arange
from numpy.random import random
from enthought.enable.example_support import DemoFrame, demo_main
# Enthought library imports
from enthought.enable.api import Component, ComponentEditor, Window
from enthought.traits.api import HasTraits, Instance, Array, Str
from enthought.traits.ui.api import Item, Group, View
# Chaco imports
from enthought.chaco.api import AbstractDataSource, ArrayPlotData, Plot, \
                                 HPlotContainer, LassoOverlay
from enthought.chaco.tools.api import LassoSelection, ScatterInspector

from enthought.chaco.api import ArrayPlotData
from enthought.chaco.tools.api import PanTool, ZoomTool


### Display Scatter Plot ######################################################

def _create_display_scatter_plot_component(ind, z):
    # Create some data
    x = ind
    y = z
    # Create a plot data obect and give it this data
    pd = ArrayPlotData()
    pd.set_data("index", x)
    pd.set_data("value", y)
    # Create the plot
    plot = Plot(pd)
    plot.plot(("index", "value"),
              type="scatter",
              marker="circle",
              index_sort="ascending",
              color="orange",
              marker_size=3,
              bgcolor="white")
    # Tweak some of the plot properties
    plot.title = "Scatter Plot"
    plot.line_width = 0.5
    plot.padding = 50
    # Attach some tools to the plot
    plot.tools.append(PanTool(plot, constrain_key="shift"))
    zoom = ZoomTool(component=plot, tool_mode="box", always_on=False)
    plot.overlays.append(zoom)
   
    return plot

class ScatterPlotView(HasTraits):
    plot = Instance(Component)

    ind = Array
    z = Array
   
    traits_view = View(
                    Group(
                        Item('plot', editor=ComponentEditor(size=(500, 500),
                                                            bgcolor='lightgray'), 
                             show_label=False),
                        orientation = "vertical"),
                    resizable=True,
                    )

    def _z_changed(self, old, new):
        print "z changed!"
        self.ind = np.arange(len(self.z))
        self.plot.plots['plot0'][0].index.set_data(self.ind)        
        self.plot.plots['plot0'][0].value.set_data(self.z) 
   
    def _plot_default(self):
        self.ind = np.arange(len(self.z))    
        return _create_display_scatter_plot_component(self.ind, self.z)
   
def scatter_plot(ind, z):
    return 


if __name__ == "__main__":
    
    _spview = ScatterPlotView()
    _spview.configure_traits()
