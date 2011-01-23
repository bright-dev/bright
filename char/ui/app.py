"""Provides a UI application on top of char.  Currently, this is provided as a way of visualizing the
cross-section database output.  It the future, it may also drive the char system and be able to spwan 
and monitor runs."""

from enthought.traits.api import HasTraits, Float, File, Str, Int, Array, Instance
from enthought.traits.ui.api import View, Item, HGroup, VGroup, InstanceEditor, spring

import tables as tb

from char.ui.hdf5_table import Hdf5Table
from char.ui.hdf5_tree import Hdf5Viewer

class Application(HasTraits):
    """Front-facing char application."""

    # Pytables Traits
    rx_h5_path = File(filter=["H5 (*.h5)|*.h5", "HDF5 (*.hdf5)|*.hdf5", "All files|*"], auto_set=False)
    rx_h5 = Instance(tb.File)

    perturbations_table = Instance(Hdf5Table)

    traits_view = View(
                    VGroup(
                        HGroup( Item('rx_h5_path', label="Path to data:", width=1.0) ),
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
    # Traits changed functions
    # 

    def _rx_h5_path_changed(self, new):
        # Close old file
        if self.rx_h5 is not None:
            self.rx_h5.close()

        # Open new file
        self.rx_h5 = tb.openFile(new, 'r')

        self.perturbations_table = Hdf5Table(h5=self.rx_h5, path_to_table="/perturbations")
