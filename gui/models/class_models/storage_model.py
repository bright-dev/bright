from traits.api import Str
from class_model import ClassModel


class Storage(ClassModel):

    name = Str("Storage")
    var = Str("strg")
    class_name = Str("Storage")

