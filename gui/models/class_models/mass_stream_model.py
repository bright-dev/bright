from class_model import ClassModel, Str, Dict

class MassStream(ClassModel):
    
    name = Str("MassStream")
    var = Str("ms")
    class_name = Str("MassStream")


    import_line = Str("from mass_stream import {name}")
    instance_line = Str("{var} = {classname}({mass_dict})")


