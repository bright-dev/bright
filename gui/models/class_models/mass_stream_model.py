from class_model import ClassModel, Str

class MassStream(ClassModel):
    
    name = Str("MassStream")
    var = Str("ms")
    class_name = Str("MassStream")
   


    instance_line = Str("{var} = {classname}({mass_dict})")

