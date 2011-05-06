import os
import BriPy
import tables as tb

def _get_lwr_data_path():
    data_dir = os.getenv("BRIGHT_DATA")
    lwr_data = data_dir + "/LWR.h5"
    return lwr_data

def _init_comp():
    """Initializes Bright for a fuel cycle component."""
    # Load isos2track
    lwr_data = _get_lwr_data_path()
    BriPy.load_isos2track_hdf5(lwr_data)

    # Write to hdf5 only
    BriPy.write_text(False)
    BriPy.write_hdf5(True)


def _get_comp_pass_number(natural_name):
    output_filename = BriPy.output_filename()

    if not os.path.isfile(output_filename):
        return 0

    of = tb.openFile(output_filename, 'r')
    mass = getattr(of.root, natural_name).IsosIn.Mass

    pass_num = len(mass)
    of.close()

    return pass_num
