from traits.api import Str
from class_model import ClassModel


class Reactor(ClassModel):
    
    name = Str("Reactor")
    var = Str("rx")
    class_name = Str("Reactor")
