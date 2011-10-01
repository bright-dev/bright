from traits.api import Str
from class_model import ClassModel


class Enrichment(ClassModel):
    
    name = Str("Enrichment")
    var = Str("enr")
    class_name = Str("Enrichment")
