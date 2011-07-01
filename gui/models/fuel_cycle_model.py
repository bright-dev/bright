from traits.api import HasTraits, Str, Dict, Set, on_trait_change, Instance, Property
from class_models.class_model import ClassModel
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
          
    def register_classes_available(self):
        """Check the class_models directory for all available models and record them into the classes_available dictionary."""
        
        localdict = {}
        
        dirlist = os.listdir('class_models')
        #Find all files within class_models directory
        for i in dirlist:
            match = re.search('(.+_model).py$', i)  #select models only
            if match:
                exec('from class_models.{0} import *'.format(match.group(1)), {}, localdict)  #import all classes and store into a local dictionary
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
        temp_list = [self.graph.get_edge_data(j,k)['walk_value'] for j,k in self.graph.edges()]
        if len(temp_list) == 0:
            walk_value = 1
        else:
            walk_value = max(temp_list) + 1

        #create edge with calculated walk_value, walked flag, and msname        
        self.graph.add_edge(node1, node2, walk_value=walk_value, walked=False, msname=ms_value)
        
        self.convert_to_script()

    def calc_comp(self, varname, msname, ms_value = None):
        """Inserts an execution line into the script.
           
           Format:
           calc_comp([name of variable to be calculated], [name of mass stream variable to be used])
           
           Note: Both variable names MUST be defined previously."""
        
        #self.graph.add_edge(self.variables[varname], self.variables[msname])    
        if ms_value is not None:
            self.add_edge(self.variables[varname], self.variables[msname], ms_value)
        else:
            self.add_edge(self.variables[varname], self.variables[msname])        
    
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
        temp_script3 = ""
        for key, value in self.graph.edges():
            temp_script3 = temp_script3 + key.add_calc(value.var) + '\n'
        self.script_execution = temp_script3   

if __name__ == "__main__":
    fcm = FuelCycleModel()
    fcm.add_instance("sr1","Storage")
    fcm.add_instance("sr2","Storage")
    #fcm.add_instance("enr1", "Enrichment")
    fcm.add_instance("ms1","MassStream",{922350:1.0})
    fcm.calc_comp("sr1","ms1", "ms_prod")
    fcm.calc_comp("sr2","ms1", "ms_prod")
    fcm.calc_comp("sr1","sr2")
    #fcm.configure_bright(write_text = False, write_hdf5 = True)  
    #fcm.configure_bright(track_isos = set([10010, 80160, 922380]))
    #fcm.remove_variable("sr1")
    #fcm.remove_variable("sr2")
    #fcm.remove_variable("ms1")
    print fcm.script
