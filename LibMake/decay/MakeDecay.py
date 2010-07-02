#! /usr/bin/env python

from tables import *
import isoname
import math

class DecayIso(IsDescription):
    fromiso     = Int32Col(pos = 0)  		#C-type int
    halflife    = Float64Col(dflt=1.0, pos = 1)	#C-type double
    decayconst  = Float64Col(dflt=1.0, pos = 2)	#C-type double
    toiso       = Int32Col(dflt=1.0, pos = 3)	#C-type int
    branchratio = Float64Col(dflt=1.0, pos = 4)	#C-type double

# Open a file in "w"rite mode
decayfile = openFile("decay.h5", mode = "w")
# Get the HDF5 root group
root = decayfile.root

#initialize decay table
decaytbl = decayfile.createTable(root, "Decay", DecayIso)

# Now, fill the table:
row = decaytbl.row

lib = open("Decay.LIB", 'r')			
for line in lib:
	ls = line.split()
	ls0name = isoname.LLZZZM_2_aazzzm(ls[0])
       	row['fromiso']	  = int(ls0name) #set fromiso

	#set halflife
	if ls[1] == "inf":
		hl = 10.0**300
	else:
	       	hl = float(ls[1])
       	row['halflife'] = hl

	#set decay constants
	row['decayconst'] = math.log(2.0) / hl
	
	#set to iso
	if ls[2] == "None":
		ti = 0
	else:
		ti = int(isoname.LLZZZM_2_aazzzm(ls[2], quiet = True))
	row['toiso'] = ti

	#set branch ratio
	row['branchratio'] = float(ls[3])

        # This injects the Record values
        row.append()
lib.close()

# Flush the table buffer
decaytbl.flush()

# Finally, close the file (this also will flush all the remaining buffers!)
decayfile.close()



