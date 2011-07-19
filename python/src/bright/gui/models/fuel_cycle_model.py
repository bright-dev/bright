from traits.api import HasTraits, Str, Dict, Set, on_trait_change, Instance, Property, Event
from bright.gui.models.class_models.class_model import ClassModel
#rework code to use nodes, then convert nodes to script

import networkx as nx
import os
import re

class FuelCycleModel(HasTraits):

    def __init__(self, *args, **kwargs):
        super(FuelCycleModel, self).__init__(*args, **kwargs)
        self.register_classes_available()
       
    #doc strings, other class models
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
        
        '''\if self.graph.number_of_nodes() == 1:
            first_node = self.graph.nodes()[0]
            print first_node'''
          
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
    
    #l = [n**2 for n in range(10) if n%2 == 0] ==> list comprehension
    #
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
    
        
        '''\if len(temp_list) == 0:
            walk_value = 1
        else:
            walk_value = max(temp_list) + 1'''
    
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

        ##################################################################
        #now find a way to remove the class from the classes_imported set#
        ##################################################################

        #remove the varname from the variables dictionary                
        del self.variables[varname]
        self.convert_to_script()
        self.graph_changed_event = True

    #def configure_bright() put in b/t imports and variables, add from bright import bright_config by default,  
    def configure_bright(self, **bright_options):
        """Adds an additional bright configuration line to the script.
           
           Format:
           configure_bright([bright_config attribute]=[custom value])"""
        
        temp_script = ""
        for key, value in bright_options.items():
            temp_script = temp_script + "bright_config." + key + " = " + repr(value) + "\n"           
        self.script_bright_config = temp_script
    
    #only watch graph, remove all sets
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
                
        #self.script_execution = self.script_execution + self.variables[varname].add_calc(msname) + '\n'
        #execution lines
        '''/temp_script3 = ""
        for key, value in self.graph.edges():
            temp_script3 = temp_script3 + key.add_calc(value.var) + '\n'
        self.script_execution = temp_script3'''

        '''\ temp_list = []
        for j,k in self.graph.edges():
            temp_list.append(self.graph[j][k]['walk_value'])
            temp_list.sort()
            for i in temp_list:
                for j,k in self.graph.edges():
                    if self.graph[j][k]['walk_value'] == i:'''
        nodes = self.graph.nodes()
        pred = {node:self.graph.predecessors(node) for node in nodes}       
        pred_len = {node:len(pred[node]) for node in nodes}
        
        success = {node:self.graph.successors(node) for node in nodes}
        success_len = {node:len(success[node]) for node in nodes}
        
        for node in nodes:
            if (pred_len[node] == 0) and (success_len[node] > 0):
                self.follow_path(node)

                
    
    def follow_path(self, node):
        temp_script3 = ""    
        var = 1
        original_node = node
        cycle_element_used = set()
        while var == 1:
            original_walked = [self.graph.edge[j][k]['walked'] for j,k in self.graph.out_edges(original_node)]
            edges = self.graph.out_edges(node)
            for j,k in edges:
                if self.graph.edge[j][k]['walked'] == True:
                    edges.remove((j,k))
                if len(nx.simple_cycles(self.graph)) > 0:
                   cycles_list = nx.simple_cycles(self.graph)[0]
                   if j in cycles_list and not cycle_element_used:
                        temp_script3 = temp_script3 + "for n in range(10):\n"
                        cycles_list.pop()
                        for i in range(len(cycles_list)):
                            temp_script3 = temp_script3 + "    " + cycles_list[i].add_calc(cycles_list[i-1].var, self.graph.edge[cycles_list[i-1]][cycles_list[i]]['msname']) + '\n'                             
                        cycle_element_used.add(j)
                        cycle_element_used.add(k)
            walk_list = {self.graph.edge[j][k]['walk_value']:(j,k) for j,k in edges}
            if len(walk_list) > 0:
                x,y =  walk_list[min(walk_list)]
                self.graph.edge[x][y]['walked'] = True
                if y.add_calc(x.var, self.graph.edge[x][y]['msname']) not in temp_script3:
                    temp_script3 = temp_script3 + y.add_calc(x.var, self.graph.edge[x][y]['msname']) + '\n'
             
                if len(self.graph.successors(node)) > 1:
                    node = y
                else:
                    node = self.graph.successors(node)[0]
            
            else:
                node = original_node
            if False not in original_walked:
                for j,k in self.graph.out_edges():
                    self.graph.edge[j][k]['walked'] = False
                break
            #else:
             #   for j,k in self.graph.out_edges():
              #      self.graph.edge[j][k]['walked'] = False
               # break
        self.script_execution = temp_script3
        




    
    
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
    print fcm.script
