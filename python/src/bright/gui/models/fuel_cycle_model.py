from traits.api import HasTraits, Str, Dict, Set, on_trait_change, Instance, Property, Event
from bright.gui.models.class_models.class_model import ClassModel

import networkx as nx
import os
import re
import ast

class FuelCycleModel(HasTraits):

    def __init__(self, *args, **kwargs):
        super(FuelCycleModel, self).__init__(*args, **kwargs)
        self.register_classes_available()
       
    graph = Instance(nx.DiGraph)
    script = Str
    script_bright_config = Str
    script_imports = Str("from bright import bright_config\n")
    script_variables = Str
    script_execution = Str
    variables = Dict
    classes_available = Dict
    classes_imported = Set
    edges_set = Set
    graph_changed_event = Event
        
    def _graph_default(self):
        dag = nx.DiGraph()
        return dag
    
    def add_node (self, node):
        """Adds a node to an existing graph.
           
           Format:
           add_node([name of node])"""
        
        self.graph.add_node(node)
        self.convert_to_script()
  
    def add_instance(self, varname, class_name, data_dict = {}):
        """Add an instance of a class to the classes_imported dictionary as well as the script_imports string variable.
           
           Format:
           add_instance([variable name], [class name], {optional dictionary for additional data1}) """
           
        # in-memory representation
     
        var = self.classes_available[class_name](var=varname, extra_data_parameter = data_dict) #class definition is stored in var
        self.classes_imported.add(class_name) #add the name of class selected to a set
        self.variables[varname] = var #store in dictionary with varname as key and var as value

        self.add_node(var)
        self.graph_changed_event = True
        
    def register_classes_available(self):
        """Check the class_models directory for all available models and record them into the classes_available dictionary."""
        
        localdict = {}
        dirlist = os.listdir(os.path.split(os.path.abspath(__file__))[0] + '/class_models')
        #Find all files within class_models directory
        for i in dirlist:
            match = re.search('(.+_model).py$', i)  #select models only
            if match:
                exec('from bright.gui.models.class_models.{0} import *'.format(match.group(1)), {}, localdict)  #import all classes and store into a local dictionary
                for key, value in localdict.items():     #check for a subclass of ClassModel
                    if issubclass(value, ClassModel) and value != ClassModel:
                       self.classes_available[key] = value   #append to classes_available dictionary
 
    @on_trait_change ('script_imports, script_bright_config, script_variables, script_execution')
    def update_script(self):
        """Update the script if any of its components (script_imports, script_bright_config, script_variables, script_execution) change."""

        self.script = self.script_imports + self.script_bright_config + self.script_variables + self.script_execution
    
    def add_edge(self, node1, node2, ms_value = None):
        """Add an edge between two specified nodes.  If an additional mass stream name is provided, add_edge() will add
           an attribute called msname in the edge.
           Format:
           add_edge([node1], [node2], [optional mass stream name]"""
               
        #calculate the walk_value for edge
        #temp_list = [self.graph.get_edge_data(j,k)['walk_value'] for j,k in self.graph.edges
        if self.graph.out_degree(node1) > 0:
            temp_list = [self.graph[j][k]['walk_value'] for j,k in self.graph.out_edges(node1)]            
            walk_value = max(temp_list) + 1
        else:
            walk_value = 1
    
        #create edge with calculated walk_value, walked flag, and msname        
        self.graph.add_edge(node1, node2, walk_value=walk_value, walked=False, msname=ms_value)
        self.convert_to_script()
        self.graph_changed_event = True

    def calc_comp(self, varname, msname, ms_value = None):
        """Inserts an execution line into the script.
           
           Format:
           calc_comp([name of variable to be calculated], [name of mass stream variable to be used])
           
           Note: Both variable names MUST be defined previously."""
        
        #self.graph.add_edge(self.variables[varname], self.variables[msname])    
        if ms_value is not None:
            self.add_edge(self.variables[msname], self.variables[varname], ms_value)
        else:
            self.add_edge(self.variables[msname], self.variables[varname])        
        self.graph_changed_event = True
    
    def remove_variable(self, varname):
        """Deletes everything related to a specified variable from memory.
           
           Format:
           remove_variable([name of variable to be deleted])"""
        self.graph.remove_node(self.variables[varname])
        
        #removes the class of varname (should there be no other variables that share that class) 
        #from classes_imported set 
        for key,value in self.classes_available.items():
            if isinstance(self.variables[varname], value):
                self.classes_imported.remove(key)           

        #remove the varname from the variables dictionary                
        del self.variables[varname]
        self.convert_to_script()
        self.graph_changed_event = True

    def configure_bright(self, **bright_options):
        """Adds an additional bright configuration line to the script.
           
           Format:
           configure_bright([bright_config attribute]=[custom value])"""
        
        temp_script = ""
        for key, value in bright_options.items():
            temp_script = temp_script + "bright_config." + key + " = " + repr(value) + "\n"           
        self.script_bright_config = temp_script
    
    @on_trait_change ('graph')
    def convert_to_script (self):
        
            
        #import line
        temp_script = "from bright import bright_config\n"
        import_set = set()
        for key, value in self.variables.items():
            import_set.add(value.add_import())
        for k in import_set:
            temp_script = temp_script + k + '\n'
        self.script_imports = temp_script
        
        #variable instantiation lines
        temp_script2 = ""
        for i in self.graph.nodes():
            temp_script2 = temp_script2 + i.add_instance() + '\n'
        self.script_variables = temp_script2
                
        nodes = self.graph.nodes()
        pred = {node:self.graph.predecessors(node) for node in nodes}       
        pred_len = {node:len(pred[node]) for node in nodes}
        
        success = {node:self.graph.successors(node) for node in nodes}
        success_len = {node:len(success[node]) for node in nodes}
        
        for node in nodes:
            if (pred_len[node] == 0) and (success_len[node] > 0):
                self.follow_path(node)

                
    
    def follow_path(self, node):
        """Given a head node, follow_path will parse through all connecting nodes and translate the path into the
           execution portion of the script. """
        #declare variables
        temp_script3 = ""    
        var = 1
        original_node = node
        cycle_element_used = set()

        while var == 1:
            #check for edges that have been traversed and ignore them
            original_walked = [self.graph.edge[j][k]['walked'] for j,k in self.graph.out_edges(original_node)]
            edges = self.graph.out_edges(node)
            for j,k in edges:
                if self.graph.edge[j][k]['walked'] == True:
                    edges.remove((j,k))
                #if a loop is detected, parse through and translate everything into a for loop for the script
                if len(nx.simple_cycles(self.graph)) > 0:
                    cycles_list = nx.simple_cycles(self.graph)[0]
                    if j in cycles_list and not cycle_element_used:
                        temp_script3 = temp_script3 + "for n in range(10):\n"
                        
                        #cycles_list.pop(), this was commented out due to the fact that the order of the cycle was incorrect
                        cycles_list.remove(cycles_list[0])
                        for i in range(len(cycles_list)):
                            temp_script3 = temp_script3 + "    " + cycles_list[i].add_calc(cycles_list[i-1].var, self.graph.edge[cycles_list[i-1]][cycles_list[i]]['msname']) + '\n'                             
                        cycle_element_used.add(j)
                        cycle_element_used.add(k)
  
            #move to the next node and add any lines of code that have not already been included in the script
            walk_list = {self.graph.edge[j][k]['walk_value']:(j,k) for j,k in edges}
            if len(walk_list) > 0:
                x,y =  walk_list[min(walk_list)]
                self.graph.edge[x][y]['walked'] = True
                if y.add_calc(x.var, self.graph.edge[x][y]['msname']) not in temp_script3:
                    temp_script3 = temp_script3 + y.add_calc(x.var, self.graph.edge[x][y]['msname']) + '\n'

                #if there are more than one edges branching off of a node, take the one with the lowest walk_value
                #else:simply take the edge and follow it to the next node
                if len(self.graph.successors(node)) > 1:
                    node = y
                else:
                    node = self.graph.successors(node)[0]
            
            #if all edges down one path has been traveres, check the head node for any other paths and follow the one
            #with the next least walk_value
            else:
                node = original_node
        
            #if everything has been traversed, reset all walked values to False to ensure that the next runthrough 
            #succeeds
            if False not in original_walked:
                for j,k in self.graph.out_edges():
                    self.graph.edge[j][k]['walked'] = False
                break

        #apply changes to the execution portion of the script
        self.script_execution = temp_script3
        
    def script_to_graph(self, fcm_script):
        """Given a script (in string format), script_to_graph creates an instance of the ScriptToParser class to
           parse through and convert the string into its graphical equivalent."""
        
        #create an instance of the parser class
        stg = ScriptToGraphParser()

        #convert the script to an abstract syntax tree
        astrep = ast.parse(fcm_script)

        #traverse through the nodes within the AST
        stg.visit(astrep)


        #print nx.simple_cycles(stg.graph_from_script)

class ScriptToGraphParser (ast.NodeVisitor):
    """The ScriptToGraphParser class takes any given script (in string format) and converts it into the proper
       graphical representation.  More specifically, the class consists of methods that will break the script
       into an abstract syntax tree and traverse through nodes relevant to the conversion process (i.e. variable
       names, arguments within expressions, etc). """
    
    #create a new instance of a directed graph derived from networkx
    graph_from_script = nx.DiGraph()

    def __init__ (self):
        pass
    
    #find all imports and print out what it is trying to import (possibly not needed if error checking is implemented)
    def visit_ImportFrom (self, stmt_importfrom):
        for alias in stmt_importfrom.names:
            print alias.name
    
    #find all assignment statements and add a node using the variable name as a reference
    def visit_Assign(self, stmt_assign):
        for alias in stmt_assign.targets:
            self.graph_from_script.add_node(alias.id)
        
        #print stmt_assign.value.func.id
    
    
    #find all 'calc' lines and add an edge between them (may be uneeded if the script is always correct with error checking)
    def visit_Expr(self, stmt_expr):
        print stmt_expr.value.func.value.id, stmt_expr.value.func.attr
        for i in stmt_expr.value.args:
            self.graph_from_script.add_edge(i.id, stmt_expr.value.func.value.id)
    
    #find all for loops and extract its statements inside, making edges to the graph (once again, error checking implementation most likely needed)
    def visit_For(self, stmt_for):
        for n in stmt_for.body:
            y = n.value.func.value.id
            
            for p in n.value.args:
                x = p.value.id
            self.graph_from_script.add_edge(x,y)


    
    
if __name__ == "__main__":
    fcm = FuelCycleModel()
    fcm.add_instance("sr2","Storage")
    fcm.add_instance("sr1","Storage")
    fcm.add_instance("sr3","Storage")
    fcm.add_instance("sr4","Storage")
    fcm.add_instance("sr5","Storage")
    fcm.add_instance("sr6","Storage")
    fcm.add_instance("ms1","MassStream",{922350:1.0})
    fcm.calc_comp("sr1","ms1")
    fcm.calc_comp("sr2","sr1", "ms_prod")
    fcm.calc_comp("sr3","sr2", "ms_tail")
    fcm.calc_comp("sr4","sr3", "ms_prod")
    fcm.calc_comp("sr2","sr4", "ms_prod23")
    fcm.calc_comp("sr5","ms1")
    fcm.calc_comp("sr6","sr5")
    #fcm.configure_bright(write_text = False, write_hdf5 = True)  
    #fcm.configure_bright(track_isos = set([10010, 80160, 922380]))
    #fcm.remove_variable("sr1")
    #fcm.remove_variable("sr2")
    #fcm.remove_variable("ms1")
    fcm.script_to_graph("for n in range(10):\n    sr3.calc(sr2.ms_prod)\n    sr4.calc(sr3.ms_tail)\n    sr2.calc(sr4.ms_prod)")
    print fcm.script
