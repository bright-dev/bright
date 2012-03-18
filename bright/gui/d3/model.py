import ast

from bright.gui.models import fuel_cycle_model

import d3_js

fcscript = """from bright import bright_config
from bright.Storage import Storage
from pyne.material import Material
from bright.Reactor import Reactor
from bright.Enrichment import Enrichment
LWRUEnrich = Enrichment()
candu0 = Storage()
LWRUReactor = Reactor()
uranium_mine0 = Material({922340: 5.5e-05, 922380: 0.992745, 922350: 0.0072})
LWRUStorage = Storage()
LWRUNu = Material({})
LWRUStorage.calc(LWRUNu)
LWRUEnrich.calc(LWRUNu)
LWRUReactor.calc(LWRUEnrich)
"""

def main():
    fcm = fuel_cycle_model.FuelCycleModel()
    stg = fuel_cycle_model.ScriptToGraphParser()
    stg.graph_from_script = fcm.graph
    astrep = ast.parse(fcscript)
    stg.visit(astrep)

    for i, (node, attrs) in enumerate(fcm.graph.node.items()):
        attrs['REC'] = i
        attrs['Label'] = "{}_{}".format(i, node)

    d3_js.export_d3_js(fcm.graph, files_dir="fcmodel", graphname="fcm", 
                       group="REC", node_labels=True)


if __name__ == "__main__":
    main()
