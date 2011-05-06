""" Defines custom BlockCanvas application object and runs it if executed as __main__. """

from enthought.contexts.api import DataContext, MultiContext

from enthought.block_canvas.app.app import Application
from enthought.block_canvas.class_tools.class_library import ClassLibrary   
from enthought.block_canvas.function_tools.function_search import FunctionSearch
from enthought.block_canvas.function_tools.function_library import FunctionLibrary


if __name__ == '__main__':
    
    code = "from bright_view.storage import Storage\n" \
        "from bright_view.lwr import LightWaterReactor\n" \
        "from bright_view.mass_stream import NaturalUranium\n" \
        "from bright_view.enrichment import Enrichment\n" \
        "nu_out_1 = NaturalUranium(nu_in_1)\n" \
        "IsosOut_1 = Enrichment(nu_out_1, enrichment_1, name_1)\n" \
        "IsosOut_2 = LightWaterReactor(IsosIn_2, burnup_1, batches_1, radius_1, IsosOut_1, flux_1, fuel_density_1, coolant_density_1, pnl_1, use_disadvantage_1, hydrogen_rescale_1, lattice_type_1, open_slots_1, total_slots_1, name_2)\n" \
        "IsosOut_3 = Storage(IsosIn_3, IsosOut_2, name_3)"

    modules = ['bright_view']

    class_library = ClassLibrary(modules=modules)
    func_library = FunctionLibrary(modules=modules)
    func_search = FunctionSearch(all_functions=func_library.functions)
    
    app = Application(
        code=code,
        data_context=DataContext(name='data'),
        class_library=class_library, 
        function_library=func_library, 
        function_search=func_search,
        )

    app.configure_traits()
