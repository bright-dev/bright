# Major library imports
import numpy as np
from numpy import arange, sort, compress, arange
from numpy.random import random
from enable.example_support import DemoFrame, demo_main
# Enthought library imports
from enable.api import Component, ComponentEditor, Window

from traits.api import HasTraits, Instance, Array, Str, File, Any, on_trait_change
from traitsui.api import Item, Group, View, Tabbed, VGroup, HGroup

# Chaco imports
from chaco.api import AbstractDataSource, ArrayPlotData, Plot, \
                                 HPlotContainer, LassoOverlay
from chaco.tools.api import LassoSelection, ScatterInspector

from chaco.api import ArrayPlotData
from chaco.tools.api import PanTool, ZoomTool

import tables as tb

from views.hdf5_tree import _hdf5_tree, _hdf5_tree_editor, Hdf5FileNode

### Display Scatter Plot ######################################################

def _fuel_cycle_plot_component(x, y, x_name, y_name):
    # Create some data

    # Create a plot data obect and give it this data
    pd = ArrayPlotData()
    pd.set_data("index", x)
    pd.set_data("value", y)
    # Create the plot
    plot = Plot(pd)
    plot.plot(("index", "value"),
              type="line",
              marker="circle",
              index_sort="ascending",
              color="red",
              marker_size=3,
              bgcolor="white")
    # Tweak some of the plot properties
    plot.title = "Fuel Cycle Plot"
    plot.line_width = 0.5
    plot.padding = 100

    plot.x_axis.title = x_name
    plot.x_axis.title_font = "Roman 16"
    plot.x_axis.tick_label_font = "Roman 12"
        
    plot.y_axis.title = y_name
    plot.y_axis.title_font = "Roman 16"
    plot.y_axis.tick_label_font = "Roman 12"

    # Attach some tools to the plot
    plot.tools.append(PanTool(plot))
    zoom = ZoomTool(component=plot, tool_mode="box", always_on=False)
    plot.overlays.append(zoom)
   
    return plot

class FuelCyclePlotView(HasTraits):
    plot = Instance(Component)

    file = File('fuel_cycle.h5')

    h5file = Any

    x_path = Str
    y_path = Str

    x_name = Str
    y_name = Str

    x = Array
    y = Array

    x_tree = Instance(Hdf5FileNode)
    y_tree = Instance(Hdf5FileNode)


    @on_trait_change('x, y, x_name, y_name')    
    def update_plot(self):
        try:
            len(self.x)
            len(self.y)
        except:
            return 

        if len(self.x) == len(self.y):
            x = self.x
            y = self.y

        elif len(self.x) < len(self.y):
            x = np.arange(len(self.y))
            y = self.y

        elif len(self.y) < len(self.x):
            x = self.x
            y = np.arange(len(self.x))
            
        fcpc = _fuel_cycle_plot_component(x, y, self.x_name, self.y_name)
        self.plot = fcpc

#    def update_x_path(self, node):
#        print node.path
#        self.x_path = node.path[1:].replace('/', '.')

    x_node = Any
    y_node = Any

    traits_view = View(
                        HGroup(
                            Item('x_tree', editor = _hdf5_tree_editor(selected='x_node'), resizable = True, show_label=False, width=0.15), 
                            Item('y_tree', editor = _hdf5_tree_editor(selected='y_node'), resizable = True, show_label=False, width=0.15), 
                            Item('plot', editor=ComponentEditor(size=(700, 700), bgcolor='lightgray'), show_label=False),
                            ),
                    resizable=True,
                    )
   
    #
    # Trait Defaults 
    #

    def _plot_default(self):
        fcpc = _fuel_cycle_plot_component(self.x, self.y, 'No Data', 'No Data')
        return fcpc

    def _h5file_default(self):
        h5file = tb.openFile(self.file, 'r')
        return h5file

    def _x_default(self):
        return np.array([])

    def _y_default(self):
        return np.array([])

    def _x_tree_default(self):
        x_tree = _hdf5_tree(self.file)
        return x_tree

    def _y_tree_default(self):
        y_tree = _hdf5_tree(self.file)
        return y_tree

    #
    # Trait Changed
    #

    def _x_node_changed(self):
        x_path = self.x_node.path[1:].replace('/', '.')
        self.x_path = x_path

    def _y_node_changed(self):
        y_path = self.y_node.path[1:].replace('/', '.')
        self.y_path = y_path

    def _x_path_changed(self):
        node_str = '.'.join(['self.h5file.root', self.x_path])
        try:
            x = eval(node_str)
            self.x = np.array(x)
            self.x_name = node_str.rpartition('.')[2]
        except:
            return

    def _y_path_changed(self):
        node_str = '.'.join(['self.h5file.root', self.y_path])
        try:
            y = eval(node_str)
            self.y = np.array(y)
            self.y_name = node_str.rpartition('.')[2]
        except:
            return 
     

#    def __del__(self):
#        self.h5file.close()
#        del super(_fuel_cycle_plot_view, self)
   
def fuel_cycle_plot():
    return 

if __name__ == "__main__":
       
    
    _fcpview = FuelCyclePlotView ()
    _fcpview.configure_traits()
