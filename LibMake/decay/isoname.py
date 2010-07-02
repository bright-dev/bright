# Converts between naming conventions for nuclides.
# aazzzm is for numerals only (923350)
# LLZZZM is for letters  as well (U-235)
import string

LLaadic = {"AC": "89", "AL": "13", "AM": "95", "SB": "51", "AR": "18", "AS": "33", "AT": "85", "BA": "56", "BK": "97", "BE": "04", "BI": "83", "BH": "107", "B": "05", "BR": "35", "CD": "48", "CS": "55", "CA": "20", "CF": "98", "C": "06", "CE": "58", "CL": "17", "CR": "24", "CO": "27", "CU": "29", "CM": "96", "DS": "110", "DB": "105", "DY": "66", "ES": "99", "ER": "68", "EU": "63", "FM": "100", "F": "09", "FR": "87", "GD": "64", "GA": "31", "GE": "32", "AU": "79", "HF": "72", "HS": "108", "HE": "02", "HO": "67", "H": "01", "IN": "49", "I": "53", "IR": "77", "FE": "26", "KR": "36", "LA": "57", "LR": "103", "PB": "82", "LI": "03", "LU": "71", "MG": "12", "MN": "25", "MT": "109", "MD": "101", "HG": "80", "MO": "42", "ND": "60", "NE": "10", "NP": "93", "NI": "28", "NB": "41", "N": "07", "NO": "102", "OS": "76", "O": "08", "PD": "46", "P": "15", "PT": "78", "PU": "94", "PO": "84", "K": "19", "PR": "59", "PM": "61", "PA": "91", "RA": "88", "RN": "86", "RE": "75", "RH": "45", "RG": "111", "RB": "37", "RU": "44", "RF": "104", "SM": "62", "SC": "21", "SG": "106", "SE": "34", "SI": "14", "AG": "47", "NA": "11", "SR": "38", "S": "16", "TA": "73", "TC": "43", "TE": "52", "TB": "65", "TL": "81", "TH": "90", "TM": "69", "SN": "50", "TI": "22", "W": "74", "U": "92", "V": "23", "XE": "54", "YB": "70", "Y": "39", "ZN": "30", "ZR": "40"}

aaLLdic = {}
for key in LLaadic.keys():
	aaLLdic[LLaadic[key]] = key
	if int(LLaadic[key])  < 10:
		aaLLdic[LLaadic[key][1]] = key

class NotANuclide(Exception):
	def __init__(self, nucwas, nucnow, quiet = False):
		if not quiet:
			print "Not a Nuclide: " + nucwas + " --> " + nucnow

def LLZZZM_2_aazzzm(nuc, quiet = False):
	newnuc = ""
	nucstr = string.upper(nuc)
	nucstr = nucstr.strip('-')

	znum = '000'
	zzz = nucstr.strip('ABCDEFGHIJKLMNOPQRSTUVWXYZ-')
	if len(zzz) == 1:
		znum = '00' + zzz
	elif len(zzz) == 2:
		znum = '0' + zzz
	else:
		znum = zzz

	if nucstr[-1] == "M":
		newnuc = znum + "1"
	elif nucstr[-1] in ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]:
		newnuc = znum + "0"
	else:
		raise NotANuclide(nucstr, newnuc, quiet)

	LL = nucstr[:-1].strip('0123456789-')

	if LL in LLaadic.keys():
		newnuc = LLaadic[LL] + newnuc
	else:
		newnuc = 'aa' + newnuc
		raise NotANuclide(nucstr, newnuc, quiet)
	return newnuc

def aazzzm_2_LLZZZM(nuc, quiet=False):
	newnuc = ""
	nucstr = str(nuc)
	nucstr = nucstr.strip('-')
	
	try:
		if nucstr[-1] == "1":
			newnuc = str(int(nucstr[-4:-1])) + "M"
		elif nucstr[-1] == "0":
			newnuc = str(int(nucstr[-4:-1]))
		else:
			raise NotANuclide(nucstr, newnuc, quiet)
	except:
		raise NotANuclide(nucstr, newnuc, quiet)

	if nucstr[:-4] in aaLLdic.keys():
		newnuc = aaLLdic[nucstr[:-4]] + newnuc
	else:
		newnuc = 'LL' + newnuc
		raise NotANuclide(nucstr, newnuc, quiet)
	return newnuc

def LLZZZM_2_aazzzm_List(nuclist):
	newnuclist = []
	for nuc in nuclist:
		newnuclist.append( LLZZZM_2_aazzzm(nuc) )
	return newnuclist

def aazzzm_2_LLZZZM_List(nuclist):
	newnuclist = []
	for nuc in nuclist:
		newnuclist.append( aazzzm_2_LLZZZM(nuc) )
	return newnuclist
	
def mixed_2_aazzzm_List(nuclist):
	"Converts a list of mixed aazzzm and LLZZZM types to just the aazzzm type."
	newnuclist = []
	for nuc in nuclist:
		newnucname = string.upper(nuc)
		newnucname = newnucname.strip('-')
		try:
			aazzzm_2_LLZZZM(nuc, quiet = True)
		except NotANuclide:
			try:
				newnucname = LLZZZM_2_aazzzm(nuc, quiet = True)
			except NotANuclide:
				raise NotANuclide(nuc, newnucname, quiet = False)
		newnuclist.append( newnucname )
	return newnuclist

def mixed_2_LLZZZM_List(nuclist):
	"Converts a list of mixed aazzzm and LLZZZM types to just the aazzzm type."
	newnuclist = []
	for nuc in nuclist:
		newnucname = string.upper(nuc)
		newnucname = newnucname.strip('-')
		try:
			LLZZZM_2_aazzzm(nuc, quiet = True)
		except NotANuclide:
			try:
				newnucname = aazzzm_2_LLZZZM(nuc, quiet = True)
			except NotANuclide:
				raise NotANuclide(nuc, newnucname, quiet = False)
		newnuclist.append( newnucname )
	return newnuclist
