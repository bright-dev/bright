import ast

from bright.gui.models import fuel_cycle_model

import d3_js

fctemplate = """from bright import bright_config
from pyne.material import Material
from bright.Enrichment import Enrichment
from bright.Reactor import Reactor
from bright.Storage import Storage
{mines}
{enrs}
{rxs}
{stors}
{calc_enrs}
{calc_rxs}
{calc_stors}
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

    # setup reactors
    rxs = ['Farley1', 'Farley2', 'PaloVerde1', 'PaloVerde2', 'PaloVerde3', 
           'ArkansasNuclearOne1', 'ArkansasNuclearOne2', 'CalvertCliffs1', 
           'CalvertCliffs2', 'Pilgrim1', 'Brunswick1', 'Brunswick2', 'Harris1', 
           'HBRobinson2', 'Perry1', 'Braidwood1', 'Braidwood2', 'Byron1', 'Byron2', 
           'Dresden1', 'Dresden2', 'Dresden3', 'LaSalleCounty1', 'LaSalleCounty2', 
           'QuadCities1', 'QuadCities2', 'Zion1', 'Zion2', 'IndianPoint1', 
           'IndianPoint2', 'BigRockPoint', 'Palisades', 'LaCrosse', 'EnricoFermi2', 
           'Catawba1', 'Catawba2', 'McGuire1', 'McGuire2', 'Oconee1', 'Oconee2', 
           'Oconee3', 'BeaverValley1', 'BeaverValley2', 'CrystalRiver3', 'StLucie1', 
           'StLucie2', 'TurkeyPoint3', 'TurkeyPoint4', 'ThreeMileIsland1', 'OysterCreek', 
           'Hatch1', 'Hatch2', 'Vogtle1', 'Vogtle2', 'RiverBend1', 'SouthTexas1', 
           'SouthTexas2', 'Clinton1', 'DuaneArnold', 'WolfCreek1', 'Waterford3', 
           'MaineYankee', 'GrandGulf1', 'CooperStation', 'NineMilePoint1', 
           'NineMilePoint2', 'Millstone1', 'Millstone2', 'Millstone3', 'Monticello', 
           'PrairieIsland1', 'PrairieIsland2', 'FortCalhoun', 'DiabloCanyon1', 
           'DiabloCanyon2', 'HumboldtBay', 'Susquehanna1', 'Susquehanna2', 'Limerick1', 
           'Limerick2', 'PeachBottom2', 'PeachBottom3', 'Trojan', 'Fitzpatrick', 
           'IndianPoint3', 'HopeCreek', 'Salem1', 'Salem2', 'REGinna', 'RanchoSeco', 
           'Summer', 'SanOnofre1', 'SanOnofre2', 'SanOnofre3', 'BrownsFerry1', 
           'BrownsFerry2', 'BrownsFerry3', 'Sequoyah1', 'Sequoyah2', 'WattsBar1', 
           'ComanchePeak1', 'ComanchePeak2', 'DavisBesse', 'Callaway', 'NorthAnna1', 
           'NorthAnna2', 'Surry1', 'Surry2', 'Columbia', 'PointBeach1', 'PointBeach2', 
           'Kewaunee', 'YankeeRowe', 'HaddamNeck', 'Cook1', 'Cook2', 'Seabrook', 
           'VermontYankee', 'GEVallecitos']
    rstr = "\n".join(["{} = Reactor()".format(rx) for rx in rxs])

    # setup rxs calc
    crxs = ["{}.calc({})".format(rx, enrs[i%len(enrs)]) for i, rx in enumerate(rxs)]
    crstr = "\n".join(crxs)

    # setup storage
    rx_stor = {rx: (rx[:-1] if rx[-1].isdigit() else rx) + "Storage" for rx in rxs}
    stors = set(rx_stor.values())
    sstr = "\n".join(["{} = Storage()".format(stor) for stor in stors])

    # setup rxs calc
    cstors = ["{}.calc({})".format(rx_stor[rx], rx) for rx in rxs]
    csstr = "\n".join(cstors)

    # fill in script template
    fcscript = fctemplate.format(mines=mstr, enrs=estr, calc_enrs=cestr, rxs=rstr, 
                                 calc_rxs=crstr, stors=sstr, calc_stors=csstr,)
    #for i, s in enumerate(fcscript.splitlines()):
    #    print i+1, s

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
        elif node in rxs:
            i = 3
        elif node in stors:
            i = 4
        else:
            i = 0

        attrs['REC'] = i

    # write to d3
    d3_js.export_d3_js(fcm.graph, files_dir="us_once_through", graphname="us_once_through", 
                       group="REC", node_labels=True, height=720)


if __name__ == "__main__":
    main()
