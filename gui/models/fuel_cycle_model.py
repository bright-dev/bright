from traits.api import HasTraits, Str, Dict
from class_models.class_model import ClassModel


import os
import re

class FuelCycleModel(HasTraits):

    def __init__(self, *args, **kwargs):
        super(FuelCycleModel, self).__init__(*args, **kwargs)
        self.register_classes_available()
       

    script = Str
    variables = Dict
    classes_available = Dict
    classes_imported = Dict


    def _variables_changed(self):
        pass 


    def add_instance(self, class_name):
        pass       

    def register_classes_available(self):
        '''\os, exec
        dirlist = os.listdir('class_models')
        
        for i in dirlist:
            if "model.py" in i:
                self.classes_available[i.split('_model.py')[0]] = i'''
        localdict = {}
        
        dirlist = os.listdir('class_models')
        
        for i in dirlist:
            match = re.search('(.+_model).py$', i)
            if match:
                exec('from class_models.{0} import *'.format(match.group(1)), {}, localdict)
                for key, value in localdict.items():
                    if issubclass(value, ClassModel):
                       self.classes_available[key] = value 

if __name__ == "__main__":
    fcm = FuelCycleModel()
   # fcm.add_instance("Storage")
   
    print fcm.classes_available
