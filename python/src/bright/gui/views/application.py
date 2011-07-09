from traits.api import HasTraits, Instance, on_trait_change, DelegatesTo
from traitsui.api import View, InstanceEditor, Item
from enable.api import ComponentEditor
from bright.gui.models.fuel_cycle_model import FuelCycleModel
from graphcanvas.api import GraphView


class Application(HasTraits):
    model = Instance(FuelCycleModel)
    graph_view = Instance(GraphView)
    _container = DelegatesTo('graph_view')
    
    traits_view = View(Item('_container',editor = ComponentEditor(), show_label = False))
    
    def _model_default(self):
        fcm = FuelCycleModel()
        fcm.add_instance("nu", "MassStream", {922380: 0.992745, 922350: 0.0072, 922340: 0.000055})
        fcm.add_instance("sr1", "Storage")
        fcm.calc_comp("sr1","nu")
        return fcm
    
    def _graph_view_default(self):
        return GraphView(graph=self.model.graph, layout='tree')
    
    @on_trait_change('model.graph')
    def update_graph_view(self):
        self.graph_view.graph = self.model.graph

    

