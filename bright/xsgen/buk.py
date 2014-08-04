"""Plugin that runs burnup-criticality calculation.
"""
from __future__ import print_function

from bright.xsgen.plugins import Plugin

class XSGenPlugin(Plugin):

    def execute(self, rc):
        for state in rc.states:
            lib = rc.engine.generate(state)
            for writer in rc.writers:
                writer.write(state, lib)
