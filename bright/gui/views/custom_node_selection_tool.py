from enable.api import BaseTool
from traits.api import HasTraits, Dict, Float

class CustomNodeSelectionTool(BaseTool):
    """ Listens for double-clicks and tries to open a traits editor on the
        graph node under the mouse.
    """
    class_views = Dict
    variables_available = Dict
    classes_available = Dict

    offset_x = Float
    offset_y = Float

    def normal_left_dclick(self, event):
        for node in self.component.components:
            if node.is_in(event.x, event.y):
                    for key, value in self.classes_available.items():
                        if self.variables_available[str(node.value)].__class__ \
                        == value:
                            #x =  self.class_views[str(key)]()
                            #import pdb; pdb.set_trace()
                            #print self.class_views
                            x =  self.class_views[str(key + "View")]()
                    x.edit_traits(kind='livemodal')
                    break
    
    
    """def normal_left_down(self, event):
        for node in self.component.components:
            if node.is_in(event.x, event.y):
                print node.value
                node.position = [event.x + 100, event.y + 100]
        #import pdb; pdb.set_trace()
        self.event_state = "moving"
        #event.window.set_pointer(self.moving_pointer)
        #event.window.set_mouse_owner(self, event.net_transform())
        #self.offset_x = event.x - self.x
        #self.offset_y = event.y - self.y
        event.handled = True
        return

    def moving_mouse_move(self, event):
        self.position = [event.x-self.offset_x, event.y-self.offset_y]
        event.handled = True
        self.request_redraw()
        return

    def moving_left_up(self, event):
        self.event_state = "normal"
        #event.window.set_pointer(self.normal_pointer)
        #event.window.set_mouse_owner(None)
        event.handled = True
        self.request_redraw()
        return"""

