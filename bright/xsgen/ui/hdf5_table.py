import re

from enthought.traits.api import HasTraits, Instance, Array, NO_COMPARE, List, Class, \
    Bool, Int, Float, Complex, Str, Unicode, Any

from enthought.traits.ui.api import View, Item, TableEditor, InstanceEditor, Group, HGroup, VGroup, spring

from enthought.traits.ui.table_column import ObjectColumn

import tables as tb
import numpy as np

dtype_to_traits = {
    'b': Bool, 
    'i': Int, 
    'u': Int, 
    'f': Float,
    'c': Complex, 
    'S': Str, 
    'a': Str, 
    'U': Unicode, 
    'V': Any,
    }

dtype_pattern = ".*?([biufcSaUV]).*?"

class Hdf5Table(HasTraits):
    """An general HDF5 table object."""

    h5 = Instance(tb.File)

    path_to_table = Str

    table_data = Array(comparison_mode=NO_COMPARE)

    fields = List(Str)
    columns = List(ObjectColumn)

    row_class = Any
    rows = List

    # 
    # Default trait handlers
    #

    def _table_data_default(self):
        tab = self.h5.getNode(self.path_to_table)
        tab_data = tab.read()
        return tab_data

    def _fields_default(self):
        return list(self.table_data.dtype.names)


    def _columns_default(self):
        c = [ObjectColumn(name=f, label=f.replace('_', ' ').capitalize()) for f in self.fields]
        return c


    def _row_class_default(self):
        """Create a custom HasTriats class for this table's rows."""
        class RowClass(HasTraits):
            pass

        for field, type in self.table_data.dtype.descr:
            m = re.match(dtype_pattern, type)
            dtype_key = m.group(1)

            ttype = dtype_to_traits[dtype_key]

            RowClass.add_class_trait(field, ttype)

        return RowClass


    def _rows_default(self):
        r = [self.row_class(**dict([(name, self.table_data[n][name]) for name in self.table_data.dtype.names]))
                for n in range(len(self.table_data))]
        return r


    # 
    # Some table selection exposure
    #

    table_selection = Any
    selection_index = Int

    def _table_selection_changed(self, new):
        # Find the index of the selection, whenever the selection changes.
        for n in range(len(self.rows)):
            if new is self.rows[n]:
                self.selection_index = n
                break


    # The view here must be a function so the columns can be figured out dynamically.
    def traits_view(self):
        return View(
            VGroup(
            HGroup( Item('path_to_table', style='readonly'), ), 
            HGroup(
                Item('rows', 
                    editor=TableEditor(
                                    columns = self.columns,
                                    deletable   = False,
                                    sort_model  = False,
                                    auto_size   = True,
                                    orientation = 'vertical',
                                    filters     = [],
                                    row_factory = 'row_class',
                                    configurable = False,
                                    editable = False,
                                    show_toolbar = False,
                                    reorderable = False,
                                    selected = 'table_selection',
                                    edit_on_first_click = False,
                                    ), 
                    show_label=False, 
                    resizable=True,
                    ),
                ),
                ),
            resizable=True,
            )


# Demonstrate the table view
if __name__ == '__main__':
    h5 = tb.openFile("/home/scopatz/MultiGroupPaper/DataXS/lwr/lwr.h5", 'r')
    p = "/perturbations"

    class HasTable(HasTraits):
        h5table = Instance(Hdf5Table)

        def _h5table_default(self):
            return Hdf5Table(h5=h5, path_to_table=p)

        traits_view = View(Item('h5table', editor=InstanceEditor(view='traits_view'), style='custom', show_label=False), 
            width=500,
            height=500,
            resizable=True,
            )
        

    ht = HasTable()
    ht.configure_traits()

    h5.close()
