import ast

from bright.gui.models import fuel_cycle_model

import d3_js

fctemplate = """from bright import bright_config
from bright.Storage import Storage
from pyne.material import Material
from bright.Reactor import Reactor
from bright.Enrichment import Enrichment
{mines}
{enrs}
LWRUReactor = Reactor()
LWRUStorage = Storage()
{calc_enrs}
LWRUReactor.calc(LWRUEnrich)
"""

def main():
    # setup mines
    mines = ['CrowButte', 'SmithRanchHighland', 'Vasquez', 'AltaMesa']
    mtemplate = "{} = Material({{922340: 5.5e-05, 922380: 0.992745, 922350: 0.0072}})"
    mstr = "\n".join([mtemplate.format(mine) for mine in mines])

    # setup enrichment facilities
    enrs = ['Paducah', 'NationalEnrichmentFacility']
    estr = "\n".join(["{} = Enrichment()".format(enr) for enr in enrs])

    # setup enrs calc
    cenrs = ["{}.calc({})".format(enrs[i%len(enrs)], mine) for i, mine in enumerate(mines)]
    cestr = "\n".join(cenrs)

    # fill in script template
    fcscript = fctemplate.format(mines=mstr, enrs=estr, calc_enrs=cestr)

    # convert script to graph
    fcm = fuel_cycle_model.FuelCycleModel()
    stg = fuel_cycle_model.ScriptToGraphParser()
    stg.graph_from_script = fcm.graph
    astrep = ast.parse(fcscript)
    stg.visit(astrep)

    # color code graph
    for node, attrs in fcm.graph.node.items():
        if node in mines:
            i = 1
        elif node in enrs:
            i = 2
        else:
            i = 0

        attrs['REC'] = i

    # write to d3
    d3_js.export_d3_js(fcm.graph, files_dir="us_once_through", graphname="us_once_through", 
                       group="REC", node_labels=True)


if __name__ == "__main__":
    main()
