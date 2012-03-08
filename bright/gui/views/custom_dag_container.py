import networkx

from enable.api import Container
from traits.api import Instance, Enum

from graph_container import CustomGraphContainer

class CustomDAGContainer(CustomGraphContainer):
    """ Enable Container for Directed Acyclic Graphs
    """
    graph = Instance(networkx.DiGraph)
    
