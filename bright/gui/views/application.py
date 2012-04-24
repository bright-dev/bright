from traits.api import HasTraits, Any, Instance, on_trait_change, \
                       DelegatesTo, Dict, List, Set, Str, File, Button, \
                       Enum, Bool, Event, Int
from traitsui.api import View, InstanceEditor, Item, Group, HGroup, \
                         VGroup, Tabbed, TreeEditor, TreeNode, CodeEditor, \
                         ShellEditor, FileEditor, TitleEditor, TableEditor, \
                         ListEditor, ListStrEditor, Handler, ToolBar, Action, \
                         MenuBar, Menu
from traitsui.file_dialog import open_file, save_file
from enable.api import ComponentEditor
from bright.gui.models.fuel_cycle_model import FuelCycleModel
from graph_view import GraphView
from enthought.plugins import ipython_shell
#from graphcanvas.api import GraphView
import os
import re
from graphcanvas.graph_node_hover_tool import GraphNodeHoverTool
from bright.gui.views.custom_node_selection_tool import CustomNodeSelectionTool
from custom_graph_node_component import CustomGraphNodeComponent
from traits.trait_handlers import BaseTraitHandler, TraitHandler
from graph_container import CustomGraphContainer
from bright.gui.views.custom_dag_container import CustomDAGContainer
from bright.gui.views.custom_graph_canvas.io_coordinate import IOPair
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
        info.object.model.add_instance("lwru_nu","Material")
        info.object.model.add_instance("lwru_enrich","Enrichment")
        info.object.model.add_instance("lwru_storage","Storage")
        info.object.model.add_instance("lwru_reactor","Reactor")
        info.object.model.calc_comp("lwru_storage","lwru_enrich")
        info.object.model.calc_comp("lwru_enrich","lwru_nu")
        info.object.model.calc_comp("lwru_reactor","lwru_enrich")
        return

    def preconfigured_lwrmox(self, info):
        info.object.model.add_instance("lwrmox_reprocess","Reprocess")
        info.object.model.add_instance("lwrmox_reactor","Reactor")
        info.object.model.add_instance("lwrmox_storage","Storage")
        info.object.model.calc_comp("lwrmox_storage","lwrmox_reprocess")
        info.object.model.calc_comp("lwrmox_reactor","lwrmox_reprocess")
        return

    def preconfigured_candu(self, info):
        info.object.model.add_instance("lwrcandu_nu","Material")
        info.object.model.add_instance("lwrcandu_reactor","Reactor")
        info.object.model.calc_comp("lwrcandu_reactor","lwrcandu_nu")
        return

    save = Action(name = "Save", action = "save_file")
    open = Action(name = "Open", action = "open_file")
    preset1 = Action(name = "lwru_fuel_cycle", action = "preconfigured_lwru")
    preset2 = Action(name = "lwrmox_fuel_cycle", action = "preconfigured_lwrmox")
    preset3 = Action(name = "candu_fuel_cycle", action = "preconfigured_candu")


    



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
    selected = Any
    instancekey = Dict
    count = Int(1)

    def register_views(self):
        localdict = {}
        local_dir = os.path.split(__file__)[0]
        comp_list = os.listdir(os.path.join(local_dir, 'component_views/'))
        comp_list.remove('views')
        comp_list.remove('__init__.py')
        
        #temporarily remove these classes (some are missing in the bright 
        #directory or code needs to be rewritten to compensate for change 
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
      
                    exec('from bright.gui.views.component_views.{name} import \
                         {view_name}View'.format(name=match.group(1), \
                         view_name=vname), {}, localdict)
                #exec('from bright.gui.views.component_views.{name} import 
                #{view_name}View'.format(name=match.group(1), 
                #view_name=vname), {}, localdict)
         

                else:
                    exec('from bright.gui.views.component_views.{name} import \
                    {view_name}View'.format(name=match.group(1), \
                    view_name=vname), {}, localdict)
     
        for key, value in localdict.items():
            self.component_views[key] = value
        #import pdb; pdb.set_trace()

        
    traits_view =View(
                     VGroup(
                        HGroup(
                            Item(
                                 'classes_list', 
                                 editor = ListStrEditor(
                                             activated = 'activated_formation', 
                                             selected = 'selected',
                                             title = 'Classes Available', 
                                             editable = False, 
                                             operations = []
                                                       ), 
                                 show_label = False, 
                                 width =.10,
                                 height = .90,
                                ),
                            Item(
                                 '_container', 
                                 editor = ComponentEditor(), 
                                 show_label = False, 
                                 resizable = True, 
                                 width =.52,
                                 height = .90
                                ),
                            Item(
                                 'script', 
                                 editor = CodeEditor(), 
                                 show_label = False, 
                                 resizable = True, 
                                 width = .38,
                                 height = .90
                                )
                            ), 
                        HGroup(
                            Item(
                                 'variables_list', 
                                 editor = ListStrEditor(
                                             title = 'Variables In Use', 
                                             editable = False, 
                                             operations = []
                                                       ), 
                                 show_label = False, 
                                 resizable = True, 
                                 width =.10,
                                 height = .10
                                ),

                            Item(
                                 'model_context', 
                                 editor = ShellEditor(share = True), 
                                 show_label = False, 
                                 width = .80,
                                 height = .10
                                )
                              )
                          ),
                  resizable = True,
                  width = 1,
                  height = 1,
                  handler = handle, 
                  #handler = TreeHandler(),
                  title = "Fuel Cycle Model",
                  menubar = MenuBar(
                                    Menu(
                                         handle.open, handle.save, 
                                         name = "File"
                                        ),
                                    Menu(
                                         handle.preset1, 
                                         handle.preset2, 
                                         handle.preset3, 
                                         name = "PreSets"
                                        )
                                    ),
                       
                    
                    )
    def _activated_formation_changed(self):	
        if self.activated_formation == "":
            return
        variable_name = self.activated_formation.strip().lower()
        #variable_name = self.selected.strip().lower()
        variable_name = '_'.join(variable_name.split(' '))
        if (variable_name + str(self.instancekey[self.activated_formation][1]))\
        in self.variables_list:
            self.instancekey[self.activated_formation][1] += 1 
        self.model.add_instance(variable_name + \
            str(self.instancekey[self.activated_formation][1]), \
            self.instancekey[self.activated_formation][0])
        self.instancekey[self.activated_formation][1] += 1
        self.activated_formation = ""
        #self.model.add_instance(self.instancekey[self.activated_formation] 
        #+ str(random.randint(0,9)), self.activated_formation) 

            


    def _model_default(self):
        fcm = FuelCycleModel()
        fcm.add_instance("uranium_mine0", "Material", \
            {922380:0.992745, 922350:0.0072, 922340:0.000055})
        fcm.add_instance("interim_storage_facility0", "Storage")
        fcm.calc_comp("interim_storage_facility0","uranium_mine0")
        self.register_views()
        return fcm

    #@on_trait_event('model.graph_changed_event')
    def update_graph_view(self):
        self.graph_view.graph = self.model.graph

        #Either one can be used; which one is better?#

        self.graph_view._graph_changed(self.graph_view, self.model.graph)


 #       self.graph_view._graph_changed(self.model.graph)
        #self.graph_view._GraphView__canvas_default(self.graph_view)
        
        self.graph_view._canvas.tools.pop(1)            
        self.graph_view._canvas.tools.append(CustomNodeSelectionTool(
                            classes_available = self.model.classes_available, 
                            variables_available = self.model.variables, 
                            class_views = self.component_views, 
                            component=self.graph_view._canvas
                                                                   )
                                            )

    def _graph_view_default(self):
        self.on_trait_event(self.update_graph_view, 'model.graph_changed_event')
        gv = GraphView(graph = self.model.graph)
        
                  
        #import pdb; pdb.set_trace()


        gv._graph_changed = _graph_changed
        gv._graph_changed(gv, self.model.graph)


    



        #Either this or the three lines above can be used; 
        #which one is better though?#
        #GraphView._graph_changed = _graph_changed
        #gv = GraphView(graph = self.model.graph)
        #gv._graph_changed(self.model.graph)




        gv._canvas.tools.pop(0)
        gv._canvas.tools.append(CustomNodeSelectionTool(
                            classes_available = self.model.classes_available, 
                            variables_available = self.model.variables, 
                            class_views = self.component_views, 
                            component=gv._canvas
                                                        )
                                )

        return gv
    def _graph_container_default(self):
        #self.graph_container.draw = draw
        return CustomDAGContainer()
    def _model_context_default(self):
        return {'fc': self.model}

    def _classes_list_default(self):

        fcm = FuelCycleModel()
        list_temp = []
	#circle_temp = {}
	#pair = IOPair()
        x = ["    Uranium Mine",
             "    Thorium Mine", 
             "    Pressurized Water Reactor", 
             "    Sodium Fast Reactor", 
             "    CANDU",
             "    Aqueous Reprocess Plant", 
             "    Electrochemical Reprocessing Plant",
             "    Interim Storage Facility",
             "    Geologic Repository"]

        for key, value in fcm.classes_available.items():
            list_temp.append(key)
            if(key == "Reactor"):
                list_temp.append(x[2])
                list_temp.append(x[3])
                list_temp.append(x[4])
            elif(key == "Material"):
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

        for i in self.classes_list:
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
        container = CustomGraphContainer(style=self.layout)

    container.tools.append(CustomNodeSelectionTool(component=container))
    container.tools.append(GraphNodeHoverTool(component=container,
                                                  callback=self._on_hover))
    return container



if __name__ == '__main__':
    app = Application()
    #app.register_views()
    app.configure_traits()
    



