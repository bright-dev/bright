# Chaco imports
from enthought.chaco.api import AbstractDataSource, ArrayPlotData, Plot, \
                                 HPlotContainer, LassoOverlay
from enthought.chaco.tools.api import LassoSelection, ScatterInspector

from enthought.chaco.api import ArrayPlotData
from enthought.chaco.tools.api import PanTool, ZoomTool

import metasci

def stairstep_plot(energy, data, data_name):
    # Munge the data into a plotable form
    x, y = metasci.stair_step(energy, data)
    x = np.array(x)
    y = np.array(y)

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
    plot.title = data_name
    plot.line_width = 0.5
    plot.padding = 100

    plot.x_axis.title = "Energy [MeV]"
    plot.x_axis.title_font = "Roman 16"
    plot.x_axis.tick_label_font = "Roman 12"
        
    plot.y_axis.title = "Data"
    plot.y_axis.title_font = "Roman 16"
    plot.y_axis.tick_label_font = "Roman 12"

    # Attach some tools to the plot
    plot.tools.append(PanTool(plot))
    zoom = ZoomTool(component=plot, tool_mode="box", always_on=False)
    plot.overlays.append(zoom)
   
    return plot
