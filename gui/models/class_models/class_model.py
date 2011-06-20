from traits.api import HasTraits, Str


class ClassModel(HasTraits):

    name = Str

    import_line = Str("from bright import {name}")


    def add_import(self):
        return import_line.format(name=self.name)
