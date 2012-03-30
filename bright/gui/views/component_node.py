import numpy
from math import sqrt

from enable.api import Component
from traits.api import Tuple, List, Int, Float, Any, Str, cached_property, Property, Set

class ComponentNode(Component):
    #import pdb; pdb.set_trace()
    x = Float()
    y = Float()
    height = Float()
    length = Float()
    value = Any
    fill_color = Tuple()

    def draw(self, gc, view_bounds=None, mode="default"):
        gc.save_state()
        gc.begin_path()
        #import pdb; pdb.set_trace()        
        gc.arc(self.x, self.y+(self.height*self.length), 3, numpy.pi/2, -numpy.pi/2)
        gc.arc(self.x, self.y+(self.height*self.length), -3, numpy.pi/2, -numpy.pi/2)      
        #gc.set_fill_color((red,green,blue,shade)) range between 0.0 to 1.0
        #import pdb; pdb.set_trace()
        gc.set_fill_color(self.fill_color)
        #gc.set_fill_color((0.8,0.0,0.1,1.0))
        gc.draw_path()  
        gc.fill_path()
        gc.restore_state()
        
    def normal_left_down(self, event):
        event.handled = True
        print "Node clicked."
        event.handled = True
        return

