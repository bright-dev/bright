from class_model import ClassModel, Str, Dict

class Material(ClassModel):
    
    name = Str("Material")
    var = Str("ms")
    class_name = Str("Material")


    import_line = Str("from pyne.material import {name}")
    instance_line = Str("{var} = {classname}({mass_dict})")


