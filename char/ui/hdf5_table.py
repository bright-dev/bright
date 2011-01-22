from enthought.traits.api import HasTraits, Str, Instance, Array, NO_COMPARE
from enthought.traits.ui.ui_editors.array_view_editor import ArrayViewEditor

import tables as tb

class Hdf5Table(HasTraits):
    """An general HDF5 table object."""

    h5 = Instance(tb.File)

    path_to_table = Str

    table_data = Array(comparison_mode=NO_COMPARE)

    editor = ArrayViewEditor(font='Mono 10')


    def _table_data_default(self):
        tab = self.h5.getNode(self.path_to_table)
        tab_data = np.array(tab)
        return tab_data

#    def _editor_default(self):
        

#class 
