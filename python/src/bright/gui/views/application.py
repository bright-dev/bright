from traits.api import HasTraits, Instance, on_trait_change, DelegatesTo, Dict
from traitsui.api import View, InstanceEditor, Item, HGroup, VGroup, Tabbed, CodeEditor, ShellEditor
from enable.api import ComponentEditor
from bright.gui.models.fuel_cycle_model import FuelCycleModel
from graphcanvas.api import GraphView


class Application(HasTraits):
    model = Instance(FuelCycleModel)
    graph_view = Instance(GraphView)
    _container = DelegatesTo('graph_view')
    script = DelegatesTo('model')
    model_context = Dict
    

    traits_view = View(
                    VGroup(
                        HGroup(
                            Item('_container', editor = ComponentEditor(), show_label = False, resizable = True, width =.25),
                            Item('script', editor = CodeEditor(), show_label = False, resizable = True, width = .50)
                            ),
                    
                        Item('model_context', editor = ShellEditor(share = True), label = 'Shell')
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
        print "yo dudes i'm workin"
        self.graph_view.graph = self.model.graph
    
    def _graph_view_default(self):
        self.on_trait_event(self.update_graph_view, 'model.graph_changed_event')
        return GraphView(graph = self.model.graph)
    
    def _model_context_default(self):
        return {'fc':self.model}

if __name__ == '__main__':
    app = Application()
    app.configure_traits()
    

