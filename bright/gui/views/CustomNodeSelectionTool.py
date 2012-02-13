from enable.api import BaseTool
from traits.api import HasTraits, Dict

class CustomNodeSelectionTool(BaseTool):
    """ Listens for double-clicks and tries to open a traits editor on the
        graph node under the mouse.
    """
    class_views = Dict
    variables_available = Dict
    classes_available = Dict


    def normal_left_dclick(self, event):
        for node in self.component.components:
            if node.is_in(event.x, event.y):
                    for key, value in self.classes_available.items():
                        if self.variables_available[str(node.value)].__class__ == value:
                            #x =  self.class_views[str(key)]()
                            #print self.class_views
                            x =  self.class_views[str(key + "View")]()
                    x.edit_traits(kind='livemodal')
                    break

