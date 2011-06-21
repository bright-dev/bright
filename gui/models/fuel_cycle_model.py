from traits.api import HasTraits, Str, Dict

import os


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
        #        

    def register_classes_available(self):
        #os, exec
        list = os.listdir('class_modles')
        
        for i in list:
            if "model" in i:
                classes_available[i.split('_models.py')[0]] = i

if __name__ == "__main__":
    fcm = FuelCycleModel()
    fcm.add_instance("Storage")

    print fcm.script
