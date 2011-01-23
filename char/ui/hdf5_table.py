import re

from enthought.traits.api import HasTraits, Instance, Array, NO_COMPARE, List, Class, \
    Bool, Int, Float, Complex, Str, Unicode, Any
from enthought.traits.ui.api import View, Item, TableEditor, InstanceEditor

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

    row_class = Class
    rows = List

    table_editor = TableEditor


    def _table_data_default(self):
        tab = self.h5.getNode(self.path_to_table)
        tab_data = tab.read()

        #self.fields = list(tab_data.dtype.names)

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


    def _table_editor_default(self):
        return  TableEditor(
            columns = self.columns,
            deletable   = False,
            sort_model  = False,
            auto_size   = False,
            orientation = 'vertical',
            filters     = [],
            row_factory = self.row_class,
            configurable = False,
            editable = False,
            show_toolbar = False,
            reorderable = False,
            selected = 'table_selection',
            edit_on_first_click = False,
            )

    traits_view = View(Item('rows', editor=_table_editor_default, show_label=False))


if __name__ == '__main__':

    from enthought.traits.ui.api import InstanceEditor

    h5 = tb.openFile("/home/scopatz/MultiGroupPaper/DataXS/lwr/lwr.h5", 'r')
    p = "/perturbations"

    class HasTable(HasTraits):
        h5table = Hdf5Table(h5=h5, path_to_table=p)
        traits_view = View(Item('h5table', editor=InstanceEditor(view='traits_view')), 
            width=500,
            height=500,
            resizable=True,
            )
        

    ht = HasTable()
    ht.configure_traits()

    h5.close()
