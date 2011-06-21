from traits.api import HasTraits, Str, Dict


class ClassModel(HasTraits):

    name = Str
    var = Str
    class_name = Str

    import_line = Str("from bright import {name}")
    instance_line = Str("{var} = {classname}()")

    def add_import(self):
        return import_line.format(name=self.name)

    def add_instance(self):
        return instance_line.format(var = self.var, classname = self.class_name)
