class NCode(object):
    """A neutronics code class"""
    name    = ""
    run_str = ""

    place_remote_files = []
    """A dummy list for putting the files needed for running the neutronics code on the remote server."""

    fetch_remote_files = []
    """A dummy list for getting the output files from running the neutronics code on the remote server."""

    def make_input(self):
        """A dummy method for making the neutronics code input."""
        pass

    def run_script_fill_values(self):
        """A dummy method for filling the run script."""
        pass

    def run(self):
        """A dummy method for running the neutronics code."""
        pass

    def parse(self):
        """A dummy method for parsing the neutronics code."""
        pass

