import warnings
import os
import sys

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from isoname import *
    from MassStream import *
    from FCComps import *

# Path to user's home directory
user_home = os.path.expanduser('~')

if '.BRIGHT_DATA' in os.listdir(user_home):
    if 'decay.h5' in os.listdir(user_home + "/.BRIGHT_DATA/"):
        os.environ["BRIGHT_DATA"] = user_home + "/.BRIGHT_DATA/"
else:
    # Find where 'decay.h5' exists on PYTHONPATH and set this as BRIGHT_DATA
    BreakOut = False
    for n in range(len(sys.path)-1, -1, -1):
        for root, dirs, files in os.walk(sys.path[n]):
            if "decay.h5" in files:
                os.environ["BRIGHT_DATA"] = root
                BreakOut = True
                break

        if BreakOut:
            break

bright_start()
