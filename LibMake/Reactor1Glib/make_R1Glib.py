from __future__ import print_function
import tables
import os
import isoname

#Reactor = "FR"
#phi = 2.0*(10**15)

#Reactor = "LWR"
#phi = 4.0*(10**14)

Reactor = "MOX"
phi = 4.0*(10**14)

#Initializes Burnup Parameters
F = []		#Fluence in [n/kb]
BUi_F_ = {}	#Burnup [MWd/kgIHM]
pi_F_ = {}	#Production rate [n/s]
di_F_ = {}	#Destruction rate [n/s]
Tij_F_ = {}	#Transformation Matrix [kg_i/kgIHM]

def toFloatList(lyst, multby = 1.0 ):
	"Converts list values to float and scales values by multby value."
	fltList = []
	if multby == 1.0:
		for el in lyst:
			fltList.append( float(el) )
	else:
		for el in lyst:
			fltList.append( float(el) * multby )
	return fltList

def loadLib(dirlib = "libs/"):
	"Loads Apporiate Libraries for Reactor and makes them into Burnup Parameters [F, pi(F), di(F), BUi(F), Tij(F)]."

	#Loads Fluence, Burnup, Production and Destruction Rates
	pdbfile = open(dirlib + "ProDesBuTable.txt", 'r')
	for line in pdbfile:
		ls = line.split()
		if len(ls) == 0:
			continue
		elif ls[0] == "Time" and len(F) == 0:
			for el in ls[1:]:
				F.append(float(el) * phi * 24.0 * 3600.0 * 1000.0 / (10.0**24) )
		elif ls[0][6:] == "BURNUP":
			BUi_F_[ls[0][:6]] = toFloatList( ls[1:] )
		elif ls[0][6:] == "PRODUC":
			pi_F_[ls[0][:6]] = toFloatList( ls[1:] )
		elif ls[0][6:] == "DESTRU":
			di_F_[ls[0][:6]] = toFloatList( ls[1:] )
		else:
			pass
	pdbfile.close()

	#Loads Transformation Matrix
	finlibdir = os.listdir(dirlib)
	for f in finlibdir:
		if f[-4:] == ".TAB":
			Tij_F_[ f[:-4] ] = {}
			tijFile = open(dirlib + f, 'r')
			for line in tijFile:
				ls = line.split()
				if len(ls) == 0:
					continue
				elif ls[0] == "Isotope":
					continue
				else:
					Tij_F_[ f[:-4] ][ ls[0] ] = toFloatList(ls[1:], 10.0**-3)
			tijFile.close()
		else:
			pass

	return

def write_H5(h5name = "reactor.h5"):		
	"Writes the reactor to an HDF5 File"

	libfile = tables.openFile(h5name, mode = "w", title = '{0} One Group Library'.format(Reactor))
	lbr = libfile.root

	libfile.createArray(lbr, "Fluence", F, "Time Integrated Flux [n/kb]") 	

	FromIso = [] 
	ToIso = [] 

	#Create Groups
	libfile.createGroup(lbr, "Burnup")
	libfile.createGroup(lbr, "Production")
	libfile.createGroup(lbr, "Destruction")
	libfile.createGroup(lbr, "Transmutation")

	#Fill Groups
	for iso in BUi_F_.keys():
		iso_LL =  isoname.zzaaam_2_LLAAAM(int(iso))
		libfile.createArray("/Burnup", iso_LL, BUi_F_[iso], "Burnup BU(F) for {0} [MWd/kgIHM]".format(iso_LL)) 	
		FromIso.append(iso_LL)

	for iso in pi_F_.keys():
		iso_LL =  isoname.zzaaam_2_LLAAAM(int(iso))
		libfile.createArray("/Production", iso_LL, pi_F_[iso], "Neutron Production Rate p(F) for {0}".format(iso_LL)) 	

	for iso in di_F_.keys():
		iso_LL =  isoname.zzaaam_2_LLAAAM(int(iso))
		libfile.createArray("/Destruction", iso_LL, di_F_[iso], "Neutron Destruction Rate d(F) for {0}".format(iso_LL)) 	

	for i in Tij_F_.keys():
		i_LL =  isoname.zzaaam_2_LLAAAM(int(i))
		libfile.createGroup("/Transmutation", i_LL)
		for j in Tij_F_[i].keys():
			j_LL =  isoname.zzaaam_2_LLAAAM(int(j))
			libfile.createArray("/Transmutation/" + i_LL, j_LL, Tij_F_[i][j], "Transmutation Matrix T(F) from {0} to {1}".format(i_LL, j_LL)) 	
			if int(i) == 922350:
				ToIso.append(j_LL)

	libfile.createArray(lbr, "FromIso_LL", FromIso, "Initial Loading Nuclide List (LL)")
	libfile.createArray(lbr, "FromIso_zz", isoname.LLAAAM_2_zzaaam_List(FromIso), "Initial Loading Nuclide List (zz)")
	libfile.createArray(lbr, "ToIso_LL",   ToIso,   "Transmutation Nuclide List (LL)")
	libfile.createArray(lbr, "ToIso_zz",   isoname.LLAAAM_2_zzaaam_List(ToIso),   "Transmutation Nuclide List (zz)")

	libfile.close()
	return

loadLib("libs/" + Reactor + "/")
write_H5(Reactor + ".h5")
