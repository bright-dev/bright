from traits.api import HasTraits, Str, Dict, Set, on_trait_change
from class_models.class_model import ClassModel


import os
import re

class FuelCycleModel(HasTraits):

    def __init__(self, *args, **kwargs):
        super(FuelCycleModel, self).__init__(*args, **kwargs)
        self.register_classes_available()
       

    script = Str
    script_imports = Str
    script_variables = Str
    script_execution = Str
    variables = Dict
    classes_available = Dict
    classes_imported = Set

    
        
    
    def add_instance(self, varname, class_name, data_dict = {}):
        # in-memory representation
        var = self.classes_available[class_name](var=varname, extra_data_parameter = data_dict) #class definition is stored in var
        self.classes_imported.add(class_name) #add the name of class selected to a set
        self.variables[varname] = var #store in dictionary with varname as key and var as value
        
        
        # script representation
        if class_name not in self.script_imports:
            self.script_imports = self.script_imports + var.add_import() + '\n'
        self.script_variables = self.script_variables + var.add_instance() + '\n'
        
        
    def register_classes_available(self):
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
    
    @on_trait_change ('script_imports, script_variables, script_execution')
    def update_script(self):
        self.script = self.script_imports + self.script_variables + self.script_execution
    
    def calc_comp(self, varname, msname):
        self.script_execution = self.script_execution + self.variables[varname].add_calc(msname) + '\n'
        
if __name__ == "__main__":
    fcm = FuelCycleModel()
    fcm.add_instance("sr1","Storage")
    fcm.add_instance("ms1","MassStream",{922350:1.0})
    fcm.calc_comp("sr1","ms1")  
    print fcm.script
