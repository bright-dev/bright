"""Provides a UI application on top of char.  Currently, this is provided as a way of visualizing the
cross-section database output.  It the future, it may also drive the char system and be able to spwan 
and monitor runs."""

from enthought.traits.api import HasTraits, Float, File, Str, Int, Array, Instance
from enthought.traits.ui.api import View, Item, HGroup, VGroup

import tables as tb


class Application(HasTraits):
    """Front-facing char application."""

    # Pytables Traits
    rx_h5_path = File
    rx_h5 = Instance(tb.File)


    traits_view = View(
                    Item('rx_h5_path', label="Path to data:")
                  )


    #
    # Traits changed functions
    # 

    def _rx_h5_file_changed(self, new):
        # Close old file
        if self.rx_h5 is not None:
            self.rx_h5.close()

        # Open new file
        self.rx_h5 = tb.openFile(new, 'r')
