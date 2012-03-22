from traits.api import HasTraits, Str, Dict


class ClassModel(HasTraits):

    name = Str
    var = Str
    class_name = Str
    extra_data_parameter = Dict
    
    
    import_line = Str("from bright.{import_name} import {name}")
    instance_line = Str("{var} = {classname}()")
    calc_line = Str("{var}.calc({ms})")
    
    def add_import(self):
        return self.import_line.format(name=self.name, import_name = self.name.lower())

    def add_instance(self):
        return self.instance_line.format(var = self.var, classname = self.class_name, mass_dict = self.extra_data_parameter)

    def add_calc(self, varname, ms_attribute = None):
        if ms_attribute is not None:
            varname = varname + "." + ms_attribute
        return self.calc_line.format(var = self.var, ms = varname)
    
    def __str__(self):
        return self.var
    
    def __repr__(self):
        return self.__str__()
