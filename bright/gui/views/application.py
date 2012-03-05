<<<<<<< HEAD
from traits.api import HasTraits, Any, Instance, on_trait_change, DelegatesTo, Dict, List, Set, Str, File, Button, Enum, Bool, Event, Int
=======
from traits.api import HasTraits, Any, Instance, on_trait_change, DelegatesTo, Dict, List, Set, Str, File, Button, Enum, Bool, Event
>>>>>>> f79612095f5f0bea9975226fcd825a16cfc6241f
from traitsui.api import View, InstanceEditor, Item, Group, HGroup, VGroup, Tabbed, TreeEditor, TreeNode, CodeEditor, ShellEditor, FileEditor, TitleEditor, TableEditor, ListEditor, ListStrEditor, Handler, ToolBar, Action, MenuBar, Menu
from traitsui.file_dialog import open_file, save_file
from enable.api import ComponentEditor
from bright.gui.models.fuel_cycle_model import FuelCycleModel
from graph_view import GraphView
#from graphcanvas.api import GraphVie
import os
import re
from graphcanvas.graph_node_hover_tool import GraphNodeHoverTool
from CustomNodeSelectionTool import CustomNodeSelectionTool
from CustomGraphNodeComponent import CustomGraphNodeComponent
from traits.trait_handlers import BaseTraitHandler, TraitHandler
from graph_container import GraphContainer
from CustomDagContainer import CustomDAGContainer
import random

class E_handler(Handler):
    file_name = File

    def save_file(self, info):
        
        if info.object.save == True:
            info.object.save = False
        else:
            info.object.save = True
        return
    def open_file(self, info):
        if info.object.open == True:
            info.object.open = False
        else:
            info.object.open = True
        return

    def preconfigured_lwru(self, info):
        info.object.model.add_instance("LWRU Nu","Material")
        info.object.model.add_instance("LWRU Enrich","Enrichment")
        info.object.model.add_instance("LWRU Storage","Storage")
        info.object.model.add_instance("LWRU Reactor","Reactor")
        info.object.model.calc_comp("LWRU Storage","LWRU Nu")
        info.object.model.calc_comp("LWRU Enrich","LWRU Nu")
        info.object.model.calc_comp("LWRU Reactor","LWRU Enrich")
        return

    def preconfigured_lwrmox(self, info):
        info.object.model.add_instance("LWRMOX Reprocess","Reprocess")
        info.object.model.add_instance("LWRMOX Reactor","Reactor")
        info.object.model.add_instance("LWRMOX Storage","Storage")
        info.object.model.calc_comp("LWRMOX Storage","LWRMOX Reprocess")
        info.object.model.calc_comp("LWRMOX Reactor","LWRMOX Reprocess")
        return

    def preconfigured_candu(self, info):
        info.object.model.add_instance("LWRCANDU Nu","Material")
        info.object.model.add_instance("LWRCANDU Reactor","Reactor")
        info.object.model.calc_comp("LWRCANDU Reactor","LWRCANDU Nu")
        return

    save = Action(name = "Save", action = "save_file")
    open = Action(name = "Open", action = "open_file")
    preset1 = Action(name = "LWRU Fuel Cycle", action = "preconfigured_lwru")
    preset2 = Action(name = "LWRMOX Fuel Cycle", action = "preconfigured_lwrmox")
    preset3 = Action(name = "CANDU Fuel Cycle", action = "preconfigured_candu")


    



class Application(HasTraits):
    model = Instance(FuelCycleModel)
    graph_view = Instance(GraphView)
    graph_container = Instance(CustomDAGContainer)
    _container = DelegatesTo('graph_view')
    script = DelegatesTo('model')
    model_context = Dict
    #graph_changed_event = DelegatesTo('model')
    #classes_available = Dict
    classes_list = List(Str)
    variables_list = List
    file_name = File
    file_name2 = File
    open = Bool
    save = Bool
    class_title = Enum('Classes Available')
    component_views = Dict
    loaded_script = Str
    handle = E_handler()
    activated_formation = Any
    instancekey = Dict
    count = Int(1)

    def register_views(self):
        localdict = {}
        comp_list = os.listdir('component_views/')
        comp_list.remove('views')
        comp_list.remove('__init__.py')
        
        #temporarily remove these classes (some are missing in the bright directory or code needs to be rewritten to compensate for change 
        comp_list.remove('scatter_plot.py')
        comp_list.remove('fuel_cycle_plot.py')
        comp_list.remove('light_water_reactor1g.py')
        comp_list.remove('enrichment.py')
        
        
        for i in comp_list:
            if 'init' not in i and 'util' not in i and 'lwr' not in i:
                match = re.search('(.+).py',i)
                vname_list = match.group(1).split("_")
                for n in vname_list:
                    vname_list[vname_list.index(n)] = n.capitalize()
                vname = ''.join(vname_list)
                if match.group(1) == 'material':
      
                    exec('from bright.gui.views.component_views.{name} import {view_name}View'.format(name=match.group(1), view_name=vname), {}, localdict)
                #exec('from bright.gui.views.component_views.{name} import {view_name}View'.format(name=match.group(1), view_name=vname), {}, localdict)
         #       else:
          #          exec('from bright.gui.views.component_views.{name} import {view_name}View'.format(name=match.group(1), view_name=vname), {}, localdict)
        for key, value in localdict.items():
            self.component_views[key] = value
        #import pdb; pdb.set_trace()

        
    traits_view =View(
                     VGroup(
                        HGroup(
                            Item('classes_list', editor = ListStrEditor(activated = 'activated_formation', title = 'Classes Available', editable = False, operations = []), show_label = False, width =.10),
                            #Item('classes_list', editor = tree_editor, resizable = False),
                            #Item('classes_list', editor = ListEditor(trait_handler=instance_handler), style = 'readonly', show_label = False, resizable = True, width =.25),
                            Item('_container', editor = ComponentEditor(), show_label = False, resizable = True, width =.52),
                            Item('script', editor = CodeEditor(), show_label = False, resizable = True, width = .38)
                            ), 
                        HGroup(
                            Item('variables_list', editor = ListStrEditor(title = 'Variables In Use', editable = False, operations = []), show_label = False, resizable = True, width =.10),
                            #Item('variables_list', editor = ListEditor(), style = 'readonly', show_label = False, resizable = True, width =.17),
                            Item('model_context', editor = ShellEditor(share = True), show_label = False, width = .80)
                              )
                          ),
                  resizable = True,
                  width = 1,
                  height = 1,
                  handler = handle, 
                  #handler = TreeHandler(),
                  title = "Fuel Cycle Model",
                  menubar = MenuBar(Menu(handle.open, handle.save, name = "File"),Menu(handle.preset1, handle.preset2, handle.preset3, name = "PreSets")),
                       
                    
                    )
    def _activated_formation_changed(self):	
        self.model.add_instance(self.activated_formation.strip() + " " + str(self.instancekey[self.activated_formation][1]), self.instancekey[self.activated_formation][0])
        self.instancekey[self.activated_formation][1] += 1
        #self.model.add_instance(self.instancekey[self.activated_formation] + str(random.randint(0,9)), self.activated_formation) 
	
        self.model.add_instance(self.instancekey[self.activated_formation] + str(random.randint(0,9)), self.activated_formation) 
        


    def _model_default(self):
        fcm = FuelCycleModel()
        fcm.add_instance("nu", "Material", {922380:0.992745, 922350:0.0072, 922340:0.000055})
        fcm.add_instance("sr1", "Storage")
        fcm.calc_comp("sr1","nu")
        self.register_views()
        return fcm

    #@on_trait_event('model.graph_changed_event')
    def update_graph_view(self):
        #print "yo dudes i'm workin"
        self.graph_view.graph = self.model.graph

        #Either one can be used; which one is better?#
        self.graph_view._graph_changed(self.graph_view, self.model.graph)
 #       self.graph_view._graph_changed(self.model.graph)
        #self.graph_view._GraphView__canvas_default(self.graph_view)
        
        self.graph_view._canvas.tools.pop(1)            
        self.graph_view._canvas.tools.append(CustomNodeSelectionTool(classes_available = self.model.classes_available, variables_available = self.model.variables, class_views = self.component_views, component=self.graph_view._canvas))

    def _graph_view_default(self):
        self.on_trait_event(self.update_graph_view, 'model.graph_changed_event')
        gv = GraphView(graph = self.model.graph)
        
                  
        #import pdb; pdb.set_trace()


        gv._graph_changed = _graph_changed
        gv._graph_changed(gv, self.model.graph)


    



        #Either this or the three lines above can be used; which one is better though?#
        #GraphView._graph_changed = _graph_changed
        #gv = GraphView(graph = self.model.graph)
        #gv._graph_changed(self.model.graph)




        gv._canvas.tools.pop(0)
        gv._canvas.tools.append(CustomNodeSelectionTool(classes_available = self.model.classes_available, variables_available = self.model.variables, class_views = self.component_views, component=gv._canvas))

        return gv
    def _graph_container_default(self):
        #self.graph_container.draw = draw
        return CustomDAGContainer()
    def _model_context_default(self):
        return {'fc': self.model}

    def _classes_list_default(self):

        fcm = FuelCycleModel()
        list_temp = []
        x = ["    Uranium Mine","    Thorium Mine", "    Pressurized Water Reactor", "    Sodium Fast Reactor", "    CANDU",
                     "    Aqueous Reprocess Plant", "    Electrochemical Reprocessing Plant","    Interim Storage Facility",
        x = ["    Uranium Mine","    Thorium Mine", "    Pressurized Water Reactor", "    Sodium Fast", "    CANDU",
                     "    Aqueous Reprocess Plant", "    Electrochemical Reprocessing Plant","    Interlm Storage Facility",
                      "    Geologic Repository"]

        for key, value in fcm.classes_available.items():
            list_temp.append(key)
            if(key == "Reactor"):
               list_temp.append(x[2])
               list_temp.append(x[3])
	       list_temp.append(x[4])
            elif(key =="Material"):
	       list_temp.append(x[0])
               list_temp.append(x[1])
            elif(key == "Reprocess"):
               list_temp.append(x[5])
               list_temp.append(x[6])
            elif(key == "Storage"):
               list_temp.append(x[7])
               list_temp.append(x[8])
            
        #import pdb; pdb.set_trace()
       # import pdb; pdb.set_trace()
        

        return list_temp
    
    def _variables_list_default(self):
        temp_list = []    
        for key, value in self.model.variables.items():
            temp_list.append(key)
        return temp_list

    def _script_changed(self):
        temp_list = []
        for key, value in self.model.variables.items():
            temp_list.append(key)
        self.variables_list  = temp_list

    def _open_changed(self):
        file_name = open_file()
        if file_name != '':
            f = open(file_name)
            for line in f:
                self.loaded_script += line
            self.model.script = self.loaded_script
            self.model.script_to_graph(self.model.script)
           
            #self.graph_view.graph = self.model.graph_from_script
            #self.file_name = file_name
        #self.open = False

    def _save_changed(self):
        file_name = save_file()
        if file_name != '':
            with open(file_name, 'w') as f:
                f.write(self.script)
        #self.save = False
        #if file_name != '':
         #   self.file_name2 = save_file
    
    def _instancekey_default(self):
        tempdict = {}
        """for i in self.classes_list:
	    if (i[0] == ' '):
	        tempdict[i] = i[4] + i[5] + i[6] + i[7] 
            else:
		tempdict[i] = i[0] + i[1] + i[2]"""
	#import pdb; pdb.set_trace()
        for i in self.classes_list:
<<<<<<< HEAD
            if ("Reactor" in i) or ("CANDU" in i):
                tempdict[i] = ["Reactor",0]
            elif "Mine" in i:
                tempdict[i] = ["Material",0]
            elif "Storage" in i or "Repository" in i:
                tempdict[i] = ["Storage",0]
            elif "Reprocess" in i:
                tempdict[i] = ["Reprocess",0]
            else:
                tempdict[i] = [i,0]
=======
	    if (i[0] == ' '):
	        tempdict[i] = i[4] + i[5] + i[6] + i[7] 
            else:
		tempdict[i] = i[0] + i[1] + i[2]
	#import pdb; pdb.set_trace()
>>>>>>> f79612095f5f0bea9975226fcd825a16cfc6241f
        return tempdict


def _graph_changed(self, new):
    for component in self._canvas.components:
        component.container = None
    self._canvas._components = []

    for node in new.nodes():
     # creating a component will automatically add it to the canvas
        CustomGraphNodeComponent(container=self._canvas, value=node)

    self._canvas.graph = new
    self._canvas._graph_layout_needed = True
    self._canvas.request_redraw()

def _GraphView__canvas_default(self):
    """ default setter for _canvas
    """
    if self.graph.is_directed():
        container = CustomDAGContainer(style=self.layout)
    else:
        container = GraphContainer(style=self.layout)

    container.tools.append(CustomNodeSelectionTool(component=container))
    container.tools.append(GraphNodeHoverTool(component=container,
                                                  callback=self._on_hover))
    return container



def draw(self, gc, view_bounds=None, mode="default"):
    if self._layout_needed:
        self.do_layout()
    # draw each component first to ensure their position and size
    # are more or less finalized
    component_dict = {}
    for component in self.components:
        component.draw(gc, view_bounds, mode)
        component_dict[component.value] = component
    # draw the connectors
    # connectors will always originate on a side
    # and terminate on the top or bottom
    line_starts = []
    line_ends = []
    for edge in self.graph.edges():
        orig = component_dict[edge[0]]
        dest = component_dict[edge[1]]
        if orig.y < dest.y:
            # up
            orig_y = orig.y + dest.height/2
            dest_y = dest.y
        else:
            # down
            orig_y = orig.y + dest.height/2
            dest_y = dest.y + dest.height

        if orig.x < dest.x:
            # right
            orig_x = orig.x + orig.width
            dest_x = dest.x + dest.width/2
        else:
            # left
            orig_x = orig.x
            dest_x = dest.x + dest.width/2

        line_starts.append([orig_x, orig_y])
        line_ends.append([dest_x, dest_y])

        with gc:
            gc.set_stroke_color((.5,.5,.5))
            gc.set_fill_color((1,1,1,0))

        # TODO: expose weighed parameters
        attributes = self.graph.get_edge_data(*edge)
        if 'weight' in attributes:
            weight = attributes['weight']
            if weight < 0.5:
                phase = 3 * 2.5;
                pattern = 3 * numpy.array((5,5))
                gc.set_line_dash(pattern,phase)
                gc.set_line_cap(CAP_BUTT)

            if self.graph.is_directed():
                gc.set_fill_color((.5,.5,.5,1))
                if orig.x < dest.x:
                    gc.arc(orig_x, orig_y, 3, -numpy.pi/2, numpy.pi/2)
                else:
                    gc.arc(orig_x, orig_y, -3, -numpy.pi/2, numpy.pi/2)

            gc.move_to(orig_x, orig_y)
            gc.line_to(dest_x, dest_y)
            gc.draw_path()

    line_starts = numpy.array(line_starts)
    line_ends = numpy.array(line_ends)


    if self.graph.is_directed():
        a = 0.707106781   # sqrt(2)/2
        vec = line_ends - line_starts
        unit_vec = vec / numpy.sqrt(vec[:,0] ** 2 + vec[:,1] ** 2)[:, numpy.newaxis]
            
        with gc:
            gc.set_fill_color((1,1,1,0))

            # Draw the left arrowhead (for an arrow pointing straight up)
            arrow_ends = line_ends - numpy.array(unit_vec*numpy.matrix([[a, a], [-a, a]])) * 10
            gc.begin_path()
            gc.line_set(line_ends, arrow_ends)
            gc.stroke_path()

            # Draw the right arrowhead (for an arrow pointing straight up)
            arrow_ends = line_ends - numpy.array(unit_vec*numpy.matrix([[a, -a], [a, a]])) * 10
            gc.begin_path()
            gc.line_set(line_ends, arrow_ends)
            gc.stroke_path()



if __name__ == '__main__':
    app = Application()
    #app.register_views()
    app.configure_traits()
    



