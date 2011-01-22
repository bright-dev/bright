"""Provides a UI application on top of char.  Currently, this is provided as a way of visualizing the
cross-section database output.  It the future, it may also drive the char system and be able to spwan 
and monitor runs."""

from enthought.traits.api import HasTraits, Float, File, Str, Int, Array
from enthought.traits.ui.api import View, Item, HGroup, VGroup

import tables as tb


class Application(HasTraits):
    """Front-facing char application."""

    rx_h5_file = File


    traits_view = View('rx_h5_file')
