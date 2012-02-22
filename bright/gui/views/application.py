from traits.api import HasTraits, Any, Instance, on_trait_change, DelegatesTo, Dict, List, Set, Str, File, Button, Enum, Bool, Event
from traitsui.api import View, InstanceEditor, Item, HGroup, VGroup, Tabbed, CodeEditor, ShellEditor, FileEditor, TitleEditor, TableEditor, ListEditor, ListStrEditor, Handler, ToolBar, Action, MenuBar, Menu
from traitsui.file_dialog import open_file, save_file
from enable.api import ComponentEditor
from bright.gui.models.fuel_cycle_model import FuelCycleModel
from graphcanvas.api import GraphView
import os
import re
from CustomNodeSelectionTool import CustomNodeSelectionTool
from CustomGraphNodeComponent import CustomGraphNodeComponent
from traits.trait_handlers import BaseTraitHandler, TraitHandler
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
    
    save = Action(name = "Save", action = "save_file")
    open = Action(name = "Open", action = "open_file")


class Application(HasTraits):
    model = Instance(FuelCycleModel)
    graph_view = Instance(GraphView)
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
    

    def register_views(self):
        localdict = {}
        comp_list = os.listdir('component_views/')
        comp_list.remove('views')
        comp_list.remove('__init__.py')
        
        #temporarily remove these classes (some are missing in the bright directory or code needs to be rewritten to compensate for change 
        comp_list.remove('scatter_plot.py')
        comp_list.remove('fuel_cycle_plot.py')
        comp_list.remove('light_water_reactor1g.py')
        
      
        for i in comp_list:
            if 'init' not in i and 'util' not in i and 'lwr' not in i:
                match = re.search('(.+).py',i)
                vname_list = match.group(1).split("_")
                for n in vname_list:
                    vname_list[vname_list.index(n)] = n.capitalize()
                vname = ''.join(vname_list)
                if match.group(1) == 'material':
      
                    exec('from pyne.{name} import {view_name}'.format(name=match.group(1), view_name=vname), {}, localdict)
                #exec('from bright.gui.views.component_views.{name} import {view_name}View'.format(name=match.group(1), view_name=vname), {}, localdict)
                else:
                    exec('from bright.{name} import {view_name}'.format(name=match.group(1), view_name=vname), {}, localdict)
        for key, value in localdict.items():
            self.component_views[key] = value
        
    traits_view = View(
                     VGroup(
                        HGroup(
                            Item('classes_list', editor = ListStrEditor(activated = 'activated_formation', title = 'Classes Available', editable = False, operations = []), show_label = False, width =.10),
                            #Item('classes_list', editor = ListEditor(trait_handler=instance_handler), style = 'readonly', show_label = False, resizable = True, width =.25),
                            Item('_container', editor = ComponentEditor(), show_label = False, resizable = True, width =.42),
                            Item('script', editor = CodeEditor(), show_label = False, resizable = True, width = .38)
                            ),
                        HGroup(
                            Item('variables_list', editor = ListStrEditor(title = 'Variables In Use', editable = False, operations = []), show_label = False, resizable = True, width =.10),
                            #Item('variables_list', editor = ListEditor(), style = 'readonly', show_label = False, resizable = True, width =.17),
                            Item('model_context', editor = ShellEditor(share = True), show_label = False, width = .80)
                              )
                          ),
                  resizable = True,
                  width = .90,
                  height = .90,
                  handler = handle,
                  title = "Fuel Cycle Model",
                  menubar = MenuBar(Menu(handle.open, handle.save, name = "File"))
                    )
    
    def _activated_formation_changed(self):
        self.model.add_instance(self.instancekey[self.activated_formation] + str(random.randint(0,9)), self.activated_formation) 
        


    def _model_default(self):
        fcm = FuelCycleModel()
        fcm.add_instance("nu", "MassStream", {922380:0.992745, 922350:0.0072, 922340:0.000055})
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
        

        self.graph_view._canvas.tools.pop(1)            
        self.graph_view._canvas.tools.append(CustomNodeSelectionTool(classes_available = self.model.classes_available, variables_available = self.model.variables, class_views = self.component_views, component=self.graph_view._canvas))
        

    def _graph_view_default(self):
        self.on_trait_event(self.update_graph_view, 'model.graph_changed_event')
        
        gv = GraphView(graph = self.model.graph)
        gv._graph_changed = _graph_changed
        gv._graph_changed(gv, self.model.graph)



        



        #Either this or the three lines above can be used; which one is better though?#
        #GraphView._graph_changed = _graph_changed
        #gv = GraphView(graph = self.model.graph)
        #gv._graph_changed(self.model.graph)




        gv._canvas.tools.pop(0)
        gv._canvas.tools.append(CustomNodeSelectionTool(classes_available = self.model.classes_available, variables_available = self.model.variables, class_views = self.component_views, component=gv._canvas))

        return gv
    
    def _model_context_default(self):
        return {'fc': self.model}

    def _classes_list_default(self):
        fcm = FuelCycleModel()
        list_temp = []
        for key, value in fcm.classes_available.items():
            list_temp.append(key)
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
        for i in self.classes_list:
            tempdict[i] = i[0] + i[1] + i[2]
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

if __name__ == '__main__':
    app = Application()
    #app.register_views()
    app.configure_traits()
    


