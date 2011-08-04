from traits.api import HasTraits, Instance, on_trait_change, DelegatesTo, Dict, List, Set, Str, File, Button, Enum
from traitsui.api import View, InstanceEditor, Item, HGroup, VGroup, Tabbed, CodeEditor, ShellEditor, FileEditor, TitleEditor, TableEditor, ListEditor
from traitsui.file_dialog import open_file, save_file
from enable.api import ComponentEditor
from bright.gui.models.fuel_cycle_model import FuelCycleModel
from graphcanvas.api import GraphView
import os
import re

class Application(HasTraits):
    model = Instance(FuelCycleModel)
    graph_view = Instance(GraphView)
    _container = DelegatesTo('graph_view')
    script = DelegatesTo('model')
    model_context = Dict
    #graph_changed_event = DelegatesTo('model')
    #classes_available = Dict
    classes_list = List    
    variables_list = List
    file_name = File
    file_name2 = File
    open = Button('Open')
    save = Button('Save')
    class_title = Enum('Classes Available')
    component_views = Dict


def register_views(self):
    localdict = {}
    comp_list = os.listdir('component_views')
    for i in dirlist:
        if i != '__init.py':
            exec('from bright.gui.views.component_views.{0} import *'.format(i), {}, localdict)
    for key, value in localdict:
        component_views[key] = value

    traits_view = View(
                    VGroup(
                        HGroup(
                            Item('open', show_label = False, width = .05),
                            Item('save', show_label = False, width = .05),
                            Item('file_name', show_label = False, width = .05)
                            
                              ),
                        HGroup(
                            Item('classes_list', editor = ListEditor(), style = 'readonly', show_label = False, resizable = True, width =.25),
                            Item('_container', editor = ComponentEditor(), show_label = False, resizable = True, width =.25),
                            Item('script', editor = CodeEditor(), show_label = False, resizable = True, width = .50)
                            ),
                        HGroup(
                            Item('variables_list', editor = ListEditor(), style = 'readonly', show_label = False, resizable = True, width =.17),
                            Item('model_context', editor = ShellEditor(share = True), show_label = False)
                              )
                          ),
                  resizable = True
                    )
    
    def _model_default(self):
        fcm = FuelCycleModel()
        fcm.add_instance("nu", "MassStream", {922380:0.992745, 922350:0.0072, 922340:0.000055})
        fcm.add_instance("sr1", "Storage")
        fcm.calc_comp("sr1","nu")
        return fcm

    #@on_trait_event('model.graph_changed_event')
    def update_graph_view(self):
        #print "yo dudes i'm workin"
        self.graph_view.graph = self.model.graph
        self.graph_view._graph_changed(self.model.graph)
        #self.graph_view = GraphView(graph =self.model.graph)
           
    def _graph_view_default(self):
        self.on_trait_event(self.update_graph_view, 'model.graph_changed_event')
        return GraphView(graph = self.model.graph)
    
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
        self.blah = temp_list

    def _open_changed(self):
        file_name = open_file()
        if file_name != '':
            self.file_name = file_name

    def _save_changed(self):
        file_name = save_file()
        if file_name != '':
            self.file_name2 = save_file
    
if __name__ == '__main__':
    app = Application()
    app.configure_traits()
    


