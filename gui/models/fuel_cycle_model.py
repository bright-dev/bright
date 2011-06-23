from traits.api import HasTraits, Str, Dict, Set
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


    def _variables_changed(self):
        pass 


    def add_instance(self, varname, class_name):
        # in-memory representation
        var = self.classes_available[class_name](var=varname)
        self.classes_imported.add(class_name)
        self.variables[varname] = var
       
        # script representation
        self.script_imports = self.script_imports + var.add_import() + '\n'


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

if __name__ == "__main__":
    fcm = FuelCycleModel()
    fcm.add_instance("sr1","Storage")
   
    print fcm.script_imports
