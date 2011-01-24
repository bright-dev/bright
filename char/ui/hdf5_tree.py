from enthought.traits.api import HasTraits, Instance, adapts, File,  Property
from enthought.traits.ui.api import View, Item
from enthought.traits.ui.api import TreeEditor, ITreeNode, ITreeNodeAdapter

import tables as tb

class TablesFile(ITreeNodeAdapter):
    adapts(tb.File, ITreeNode)
    
    #-- ITreeNodeAdapter Method Overrides --------------------------------------

    def get_name(self):
        return self.adaptee.filename
    
    def get_children(self):
        return self.adaptee.listNodes(where='/')
    
    def allows_children ( self ):
        return True

    def has_children ( self ):
        """ Returns whether the object has children.
        """
        return True

    def get_label ( self ):
        """ Gets the label to display for a specified object.
        """
        return self.adaptee.filename 
        
    def get_icon ( self, is_expanded ):
        """ Returns the icon for a specified object.
        """
        if is_expanded:
            return '<open>'
            
        return '<open>'

    def can_auto_close ( self ):
        """ Returns whether the object's children should be automatically 
            closed.
        """
        return True

    def select(self):
        return self.adaptee._v_pathname


class TablesGroup(ITreeNodeAdapter):
    adapts(tb.Group, ITreeNode)

    def get_name(self):
        return self.adaptee._v_name
    
    def get_children(self):
        groups =  self.adaptee._v_groups.values()
        leaves  = self.adaptee._v_leaves.values()
        return groups + leaves
    
    def allows_children ( self ):
        return True    
    
    def get_label ( self ):
        """ Gets the label to display for a specified object.
        """
        return self.adaptee._v_name 
    
    def has_children ( self ):
        """ Returns whether the object has children.
        """
        children_count = len(self.adaptee._v_groups)  + \
                       len(self.adaptee._v_leaves) 
        if children_count == 0:
            has_children = False
        else:
            has_children = True
            
        return has_children
    
    def get_icon ( self, is_expanded ):
        """ Returns the icon for a specified object.
        """
        if is_expanded:
            return '<open>'
            
        return '<open>'
    

class TablesArray(ITreeNodeAdapter):
    adapts(tb.Array, ITreeNode)
    
    def get_name(self):
        return self.adaptee._v_name
    
    def has_children(self):
        return False
    
    def get_label ( self ):
        """ Gets the label to display for a specified object.
        """
        return self.get_name() 
    
    def get_icon ( self, is_expanded ):
        """ Returns the icon for a specified object.
        """
        return '<item>'
   
    def select(self):
        return self.adaptee._v_pathname


class TablesTable(ITreeNodeAdapter):
    adapts(tb.Table, ITreeNode)
    
    def get_name(self):
        return self.adaptee._v_name
    
    def has_children(self):
        return False
    
    def get_label ( self ):
        """ Gets the label to display for a specified object.
        """
        return self.get_name() 
    
    def get_icon ( self, is_expanded ):
        """ Returns the icon for a specified object.
        """
        return '<item>'
   
    def select(self):
        return self.adaptee._v_pathname


class Hdf5Viewer(HasTraits):
    
    tableFile = Instance(tb.File)
    
    traits_view = View(
                       Item('tableFile', 
                            editor = TreeEditor(editable=False, 
                                                auto_open = 1, 
                                                ),
                            show_label=False, resizable=True),
                       width = 0.33,
                       height = 0.55,
                       resizable  =True
                       )
    
