"""Provides a UI application on top of char.  Currently, this is provided as a way of visualizing the
cross-section database output.  It the future, it may also drive the char system and be able to spwan 
and monitor runs."""

from enthought.traits.api import HasTraits, Float, File, Str, Int, Array, Instance, Any, on_trait_change
from enthought.traits.ui.api import View, Item, HGroup, VGroup, InstanceEditor, spring, TreeEditor

from enthought.enable.api import Component, ComponentEditor

import numpy as np
import tables as tb

from char.ui.hdf5_table import Hdf5Table
from char.ui.hdf5_tree import TablesFile, TablesGroup, TablesArray, TablesTable
from char.ui.stairstep_plot import stairstep_plot

import os

class Application(HasTraits):
    """Front-facing char application."""

    # Pytables Traits
    rx_h5_path = File(filter=["H5 (*.h5)|*.h5", "HDF5 (*.hdf5)|*.hdf5", "All files|*"], auto_set=False)
    rx_h5 = Instance(tb.File)

    # Object selected from the tree view
    tree_selected = Any 

    # Perturbation table
    perturbations_table = Instance(Hdf5Table)

    # Plot traits 
    plot = Instance(Component)

    traits_view = View(
                    VGroup(
                        HGroup( Item('rx_h5_path', label="Path to data:", width=1.0) ),

                        Item("_"),

                        HGroup(
                           Item('rx_h5',
                                editor = TreeEditor(editable=False,
                                                    hide_root=True, 
                                                    lines_mode='on',
                                                    selected='tree_selected',
                                                    ),
                                show_label=False, 
                                width = 0.15, 
                                ),

                            Item('plot', 
                                editor = ComponentEditor(size=(500, 500), bgcolor='lightgray'), 
                                show_label=False, 
                                resizable=True,
                                ),
                            ),

                        Item("_"),

                        Item('perturbations_table', 
                            editor=InstanceEditor(view='traits_view'), 
                            style='custom', 
                            show_label=False,
                            resizable=True, 
                            ),
                        ),

                    width=500,
                    height=500,
                    resizable=True,
                  )

    #
    # Trait default handlers
    #

    def _plot_default(self):
        E_g = np.array([0.1, 1.0])
        data = np.array([1.0])
        return stairstep_plot(E_g, data, "No Data Selected")


    def _rx_h5_default(self):
        f = tb.File('none.h5', 'w')
        os.remove('none.h5')
        return f

    #
    # Traits changed functions
    # 

    def _rx_h5_path_changed(self, new):
        # Close old file
        if self.rx_h5 is not None:
            self.rx_h5.close()

        # Open new file
        self.rx_h5 = tb.openFile(new, 'r')

        self.perturbations_table = Hdf5Table(h5=self.rx_h5, path_to_table="/perturbations")


    @on_trait_change('tree_selected, perturbations_table.selection_index')
    def change_plot(self):
        # Selected values
        node = self.tree_selected
        n = self.perturbations_table.selection_index

        # Don't change the plot under certain conditions
        if not isinstance(node, tb.Array):
            return 

        if len(node) != len(self.perturbations_table.table_data):
            return 

        # Read in the energy array
        if 'hi_res' in node._v_pathname:
            E_g = self.rx_h5.root.hi_res.energy.read()
        else:
            E_g = np.array(self.rx_h5.root.energy[n])

        # Read in the data value
        data = np.array(node[n])

        # Confirm that the data is graphable
        if data.shape == ():
            return 

        if len(E_g) != len(data) + 1:
            return 

        # Plot the data
        self.plot = stairstep_plot(E_g, data, node._v_pathname)

