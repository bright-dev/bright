from traits.api import HasTraits, Str, Dict



class FuelCycleModel(HasTraits):

    def __init__(self, *args, **kwargs):
        super(FuelCycleModel, self).__init__(*args, **kwargs)
        # Do registration

    script = Str

    variables = Dict
    classes_available = Dict
    classes_imported = Dict


    def _variables_changed(self):
        pass 


    def add_instance(self, class_name):
        # 



if __name__ == "__main__":
    fcm = FuelCycleModel()
    fcm.add_instance("Storage")

    print fcm.script
