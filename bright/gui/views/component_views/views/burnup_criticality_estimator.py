"""A Burnup criticality estimator"""
# Major library imports
import numpy as np

# Enthought library imports
from chaco.api import create_line_plot, OverlayPlotContainer, Plot, PlotLabel
from chaco.tools.api import PanTool, ZoomTool
from chaco.tools.cursor_tool import CursorTool, BaseCursorTool
from enable.component_editor import ComponentEditor
from traits.api import HasTraits, Instance, on_trait_change
from traitsui.api import View, Item, HGroup, VGroup

from metasci import SolveLine

def _k_BU(BU_F_, k_F_, BU):
    """Finds the k from a set of data from and a BU value."""

    BU_ind = np.where(BU <= BU_F_)[0]
    if len(BU_ind) < 2:
        return -1.0 

    BU_val = BU_F_[BU_ind]
    k_val  = k_F_[BU_ind]

    k = SolveLine(BU, BU_val[1], k_val[1], BU_val[0], k_val[0])

    return k

def _burnup_criticality_estimator_plot(BU_F_, k_F_, BUt, batches):

    plot = OverlayPlotContainer(padding=40)
    
    # Assign the data
    index = np.array(BU_F_)
    value = np.array(k_F_)
    ones = np.ones(len(BU_F_), dtype=float)
        
    # Create a LinePlot instance and add it to the subcontainer
    line = create_line_plot([index, value], add_grid=True, add_axis=True, index_sort='ascending', orientation = 'h', color='black')

    line.x_axis.title = 'Burnup [MWd/kgIHM]'
    line.x_axis.title_font = "Roman 16"
    line.x_axis.tick_label_font = "Roman 12"

    line.y_axis.title = 'Multiplication Factor'
    line.y_axis.title_font = "Roman 16"
    line.y_axis.tick_label_font = "Roman 12"

    plot.padding_left = 50

    #k_eq_1 = create_line_plot([index, ones])

    plot.add(line)
#    plot.add(k_eq_1)

    def update_positions(new):
        # Find the index, n, of the line that was changed
        for ol in line.overlays:
            if ol.current_position == new:
                n = line.overlays.index(ol)
                break

        new_BUt = new[0] * float(batches)/float(n+1)

        partial_bu = np.arange(0.0, new_BUt * (batches+1)/batches, new_BUt/batches)[1:]
        partial_k = np.array( [_k_BU(BU_F_, k_F_, BU) for BU in partial_bu] )

        current_bu = np.array([line.overlays[b].current_position[0] for b in range(batches)])
        current_k  = np.array([line.overlays[b].current_position[1] for b in range(batches)])

        if (partial_bu == current_bu).all() and (partial_k == current_k).all():
                return 

        if n == 0:
            line.overlays[batches-1].current_position = (partial_bu[batches-1], partial_k[batches-1])
        else:
            line.overlays[n-1].current_position = (partial_bu[n-1], partial_k[n-1])

    # Do the same for the sub components
    partial_bu = np.arange(0.0, BUt, BUt/batches)[1:]
    for b in range(len(partial_bu)):

        csr = CursorTool(line, 
            drag_button="left", 
            color='red', 
            line_style='dash',
            show_value_line=False, 
            )
        csr.current_position = partial_bu[b], _k_BU(BU_F_, k_F_, partial_bu[b])
        line.overlays.append(csr)
        csr.on_trait_change(update_positions, 'current_position')
        
    # Give me a cursor
    cursor = CursorTool(line, drag_button="left", color='red', show_value_line=False)

    # and set it's initial position (in data-space units)
    cursor.current_position = BUt, _k_BU(BU_F_, k_F_, BUt)
        
    # This is a rendered component so it goes in the overlays list
    line.overlays.append(cursor)

    # Link it up
    cursor.on_trait_change(update_positions, 'current_position')
        
    # Some other standard tools
    line.tools.append(PanTool(line, drag_button="right"))
    line.overlays.append(ZoomTool(line))

    plot.overlays.append(PlotLabel('Burnup-Criticality Estimator',
        component=plot,
        font = "Roman 20",
        overlay_position="top")
        )
        
    return plot
    

if __name__=='__main__':
    from traits.api import Array

    class BurnupView(HasTraits):
        burnup_plot = Instance(OverlayPlotContainer)

        traits_view = View(Item('burnup_plot',
                    editor=ComponentEditor(),
                    resizable=True,
                    springy=True,
                    show_label=False
                    ),
                    width = 600, 
                    height = 600, 
                    )

        def _burnup_plot_default(self):
            burnup = np.arange(20.0, 80.0, 60.0/100.0)
            k = np.arange(0.8, 1.2, 0.4/100.0)[::-1]

            bp = _burnup_criticality_estimator_plot(burnup, k, 60.0, 3)
            return bp

    bv = BurnupView()
    bv.configure_traits()
#    bv.edit_traits()
