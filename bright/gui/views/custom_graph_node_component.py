import numpy
from math import sqrt

from enable.api import Component
from kiva.constants import MODERN
from kiva.fonttools import Font
from traits.api import List, Int, Any, Str, cached_property, Property
from traitsui.api import View, Item, spring, HGroup


class CustomGraphNodeComponent(Component):
    """ An Enable Component which represents a graph node.
    """

    # The level from the root. This is used for layout and may not be
    # meaningful in graphs with no root level.
    level = Int(0)

    # The object contained in the graph node
    value = Any

    # The label which will be shown on the graph node
    label = Property(Str, depends_on='value')

    # The key on the graph for this node. This should not be
    # changed
    _key = Any

    """padding_left = 5
    padding_right = 5
    padding_top = 5
    padding_bottom = 5"""

    padding_left = 25
    padding_right = 25
    padding_top = 25
    padding_bottom = 25

    

    traits_view = View(HGroup(
                           spring,
                           Item('value', style='readonly', show_label=False),
                           spring),
                        width=200, resizable=True)

    def draw(self, gc, view_bounds=None, mode="default"):
        """ Draws the graph node
       """
       
        font = Font(family=MODERN)
        gc.set_font(font)        

        # update the size to match the text extent.
        x, y, width, height = gc.get_text_extent(self.label)
        #import pdb; pdb.set_trace()
        
        self.width = width + self.padding_left + self.padding_right
        self.height = height + self.padding_bottom + self.padding_top 


        self._draw_border(gc, view_bounds, mode)
        self._draw_text(gc, view_bounds, mode)
       
    def _draw_text(self, gc, view_bounds, mode):
        pos = (self.x + self.padding_left,
               self.y2 - self.padding_bottom)

        gc.show_text(self.label, pos)


    def _draw_border(self, gc, view_bounds, mode):
        """ Draws a nicely shaded border around the graph node
        """
        end_radius = 4
        """if 'MassStream' in str(self.value.__class__):
            starting_color = numpy.array([0.0, 0.0, 1.0, 1.0, 0.0])
            ending_color = numpy.array([1.0, 0.0, 0.0, 0.0, 1.0])
        elif 'Storage' in str(self.value.__class__):
            starting_color = numpy.array([0.0, 0.5, 0.5, 1.0, 1.0])
            ending_color = numpy.array([1.0, 0.2, 0.2, 0.2, 1.0])
        elif 'Reactor' in str(self.value.__class__):
            starting_color = numpy.array([0.6, 0.0, 1.0, 1.0, 1.0])
            ending_color = numpy.array([1.0, 0.2, 0.0, 0.0, 1.0])
        elif 'Reprocess' in str(self.value.__class__):
            starting_color = numpy.array([0.0, 0.5, 1.0, 1.0, 1.0])
            ending_color = numpy.array([0.8, 0.0, 0.0, 0.0, 1.0])
        elif 'Enrichment' in str(self.value.__class__):
            starting_color = numpy.array([0.2, 0.0, 0.8, 1.0, 1.0])
            ending_color = numpy.array([1.0, 0.2, 1.0, 0.2, 1.0])

        else:    
            starting_color = numpy.array([0.0, 1.0, 1.0, 1.0, 1.0])
            ending_color = numpy.array([1.0, 0.0, 0.0, 0.0, 1.0])"""
        starting_color = numpy.array([0.0, 1.0, 1.0, 1.0, 1.0])
        ending_color = numpy.array([1.0, 0.0, 0.0, 0.0, 1.0])


        x = self.x
        y = self.y

        gc.save_state()
        gc.begin_path()
        
	#gc.move_to(x,y)
	#gc.move_to(x + end_radius, y)
	#import pdb; pdb.set_trace()
	if y == 0:
            gc.arc(x, (y + self.height)/2, 2, numpy.pi/2, -numpy.pi/2)
            gc.arc(x, (y + self.height)/2, -2, numpy.pi/2, -numpy.pi/2)
	if y != 0:
            gc.arc(x, (self.height/2), 2, numpy.pi/2, -numpy.pi/2)
            gc.arc(x, (self.height/2), -2, numpy.pi/2, -numpy.pi/2)

	#import pdb; pdb.set_trace()

        gc.move_to(x + end_radius, y)
	#gc.rect(x,y,x+self.width,y+self.height)       
	
        gc.arc_to(x + self.width, y,
                x + self.width, y + end_radius,
                end_radius)
        gc.arc_to(x + self.width, y + self.height,
                x + self.width - end_radius, y + self.height,
                end_radius)
        
        gc.arc_to(x, y + self.height,
                x, y + self.height - end_radius,
                end_radius)
        
        gc.arc_to(x, y,
                x + end_radius, y,
                end_radius)
	
        gc.linear_gradient(x, y, x, y+100,
                numpy.array([starting_color, ending_color]),
                "pad")
	
	#import pdb; pdb.set_trace()

        #gc.set_fill_color((0.8,0.0,0.1,1.0))
        #gc.set_fill_color(color)
        gc.draw_path()
        gc.restore_state()
        #gc.fill_path()
    def __key_default(self):
        return self.value

    @cached_property
    def _get_label(self):
        if hasattr(self.value, 'label'):
            text = self.value.label
        else:
            text = str(self.value)
        #if len(text) > 20:
         #   text =  text[0:17] + "..."
        #if len(text) > 20:
         #   text =  text[0:17] + "\n    " + text[18:len(text)]
        return text

    def _value_changed(self):
        self.request_redraw()

