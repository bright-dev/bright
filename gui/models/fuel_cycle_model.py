from traits.api import HasTraits, Str, Dict, Set, on_trait_change, Instance
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
    import_set = Set
    edges_set = Set
    
    def _graph_default(self):
        return nx.DiGraph()
   
    
    def add_instance(self, varname, class_name, data_dict = {}):
        """Add an instance of a class to the classes_imported dictionary as well as the script_imports string variable.
           
           Format:
           add_instance([variable name], [class name], {optional dictionary for additional data1}) """
           
        # in-memory representation
        var = self.classes_available[class_name](var=varname, extra_data_parameter = data_dict) #class definition is stored in var
        self.classes_imported.add(class_name) #add the name of class selected to a set
        self.variables[varname] = var #store in dictionary with varname as key and var as value
        
        self.import_set.add(var.add_import())
        self.graph.add_node(var)
        
          
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

    
    
    def calc_comp(self, varname, msname):
        """Inserts an execution line into the script.
           
           Format:
           calc_comp([name of variable to be calculated], [name of mass stream variable to be used])
           
           Note: Both variable names MUST be defined previously."""
        
        self.graph.add_edge(self.variables[varname], self.variables[msname])    
        self.edges_set = set(self.graph.edges())
        
 
    def remove_variable(self, varname):
        """Deletes everything related to a specified variable from memory.
           
           Format:
           remove_variable([name of variable to be deleted])"""
        
        #declare local variables
        instance_count = 0
        classname_to_delete = ""
                           
        for key in self.variables:
            if self.variables[key].__class__ == self.variables[varname].__class__:
                instance_count = instance_count + 1
                
        #if there are no other similar instances, remove the import line
        if instance_count == 1:
            for key in self.classes_available:
                if self.classes_available[key] == self.variables[varname].__class__:
                    classname_to_delete = key
            imports_set = self.script_imports.splitlines()
            for i in imports_set:
                if classname_to_delete in i:
                    imports_set.remove(i)
                    self.script_imports = "\n".join(imports_set) +'\n'
        
            #remove the class from the classes_imported set
            self.classes_imported.remove(classname_to_delete)
            
        
        self.graph.remove_node(self.variables[varname])
   
        #remove the variable instantiation line for varname
        variables_set = self.script_variables.splitlines()
        for i in variables_set:
            if varname in i:
                variables_set.remove(i)
                self.script_variables = "\n".join(variables_set) + '\n'
        
        
        #remove the execution line for varname
        execution_set = self.script_execution.splitlines()
        for i in execution_set:
            if varname in i:
                execution_set.remove(i)
                self.script_execution = "\n".join(execution_set) + '\n'

        #remove the varname from the variables dictionary                
        del self.variables[varname]
 
    #def configure_bright() put in b/t imports and variables, add from bright import bright_config by default,  
    def configure_bright(self, **bright_options):
        temp_script = ""
        for key, value in bright_options.items():
            temp_script = temp_script + "bright_config." + key + " = " + repr(value) + "\n"           
        self.script_bright_config = temp_script
    
    @on_trait_change ('import_set, edges_set')
    def convert_to_script (self):
        #script_imports representation (during add_instance, add it onto a set and print out on script)
        temp_script = "from bright import bright_config\n"
        temp_script2 = ""
        temp_script3 = ""
        #import line
        for k in self.import_set:
            temp_script = temp_script + k + '\n'
        self.script_imports = temp_script
        
        #variable instantiation lines
        for i in self.graph.nodes():
            temp_script2 = temp_script2 + i.add_instance() + '\n'
        self.script_variables = temp_script2
        #self.script_execution = self.script_execution + self.variables[varname].add_calc(msname) + '\n'
        #execution lines
        for key, value in self.edges_set:
            temp_script3 = temp_script3 + key.add_calc(value.var) + '\n'
        self.script_execution = temp_script3   

if __name__ == "__main__":
    fcm = FuelCycleModel()
    fcm.add_instance("sr1","Storage")
    fcm.add_instance("sr2","Storage")
    #fcm.add_instance("enr1", "Enrichment")
    fcm.add_instance("ms1","MassStream",{922350:1.0})
    fcm.calc_comp("sr1","ms1")
    fcm.calc_comp("sr2","ms1")
    #fcm.configure_bright(write_text = False, write_hdf5 = True)  
    #fcm.configure_bright(track_isos = set([10010, 80160, 922380]))
    #fcm.remove_variable("sr1")
    #fcm.remove_variable("sr2")
    #fcm.remove_variable("ms1")
    print fcm.script
