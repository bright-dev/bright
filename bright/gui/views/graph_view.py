import networkx

from enable.api import ComponentEditor, Scrolled,Viewport
from enable.tools.api import ViewportPanTool
from traits.api import HasTraits, Instance, Dict, Any, Enum, \
        on_trait_change, Property, cached_property, List
from traitsui.api import View, Item

from custom_dag_container import CustomDAGContainer
#from bright.gui.views.custom_graph_canvas.graph_container import GraphContainer
#from bright.gui.views.graph_container import GraphContainer
from graph_container import CustomGraphContainer
from bright.gui.views.custom_graph_node_component import CustomGraphNodeComponent
from bright.gui.views.custom_node_selection_tool import CustomNodeSelectionTool
from graphcanvas.graph_node_hover_tool import GraphNodeHoverTool

def graph_from_dict(d):
    """ Creates a NetworkX Graph from a dictionary

    Parameters
    ----------
    d : dict

    Returns
    -------
    Graph: NetworkX Graph

    Examples
    --------
    >>> g = graph_from_dict({'a':['b'], 'b':['c', 'd'], 'c':[], 'd':[], 'e':['d']})
    """

    g = networkx.DiGraph()
    for key, children in d.items():
        for child in children:
            g.add_edge(key, child)
    return g

class GraphView(HasTraits):
    """ View containing visualization of a networkx graph.
    """
    # The graph to be visualized
    graph = Instance(networkx.Graph)
    nodes = Property(List, depends_on='graph')

    # How the graph's visualization should be layed out
    layout = Enum('spring', 'tree', 'shell', 'circular')

    
    # Scrolled contained which holds the canvas in a viewport
    _container = Instance(Scrolled)
    # The canvas which the graph will be drawn on
    _canvas = Instance(CustomDAGContainer)

    traits_view = View(Item('_container', editor=ComponentEditor(),
                            show_label=False),
                        width=400,
                        height=400,
                        resizable=True)

    def __init__(self, *args, **kw):
        super(GraphView, self).__init__(*args, **kw)

        if isinstance(self.graph.nodes()[0], HasTraits):
            self.on_trait_change(self.node_changed, 'nodes.+')

    def __canvas_default(self):
        """ default setter for _canvas
        """
        if self.graph.is_directed():
            #container = CustomDAGContainer(style=self.layout)
            container = CustomDAGContainer(style=self.layout)
        else:
            container = CustomGraphContainer(style=self.layout)

        container.tools.append(CustomNodeSelectionTool(component=container))
        container.tools.append(GraphNodeHoverTool(component=container,
                                                  callback=self._on_hover))

        return container

    def __container_default(self):
        """ default setter for _container
        """

        viewport = Viewport(component=self._canvas, enable_zoom=True)
        viewport.view_position = [0,0]
        viewport.tools.append(ViewportPanTool(viewport))

        return Scrolled(self._canvas,
                        viewport_component = viewport)

    @cached_property
    def _get_nodes(self):
        return self.graph.nodes()

    def _graph_changed(self, new):
        #""" handler for changes to graph attribute
        #"""
        for component in self._canvas.components:
            component.container = None

        self._canvas._components = []
        
        for node in new.nodes():
            # creating a component will automatically add it to the canvas
            CustomGraphNodeComponent(container=self._canvas, value=node)
        
        self._canvas.graph = new
        self._canvas._graph_layout_needed = True
        self._canvas.request_redraw()
        
    def _layout_changed(self, new):
        self._canvas.style = new

    def _on_hover(self, label):
        print "hovering over:", label

#    @on_trait_change('nodes.+')
    def node_changed(self, name, obj, old, new):
        print "node changed"
        self._canvas.request_redraw()
