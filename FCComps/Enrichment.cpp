// Enrichment Component Class

#include "Enrichment.h"

/************************************************/
/*** Enrichment Component Class and Functions ***/
/************************************************/

/***************************/
/*** Protected Functions ***/
/***************************/

void Enrichment::initialize(SepEffDict sed)
{
    //Initializes the enrichment component. with specific separation efficiencies.
}

/*******************************/
/*** Enrichment Constructors ***/
/*******************************/

Enrichment::Enrichment ()
{
    //Enrichment Fuel Cycle Component. Perfomrs Multi-Component Enrichment
}

Enrichment::Enrichment (std::string n) : FCComp (rep_p2track, n)
{
    //Enrichmenting Fuel Cycle Component.  Applies Separation Efficiencies.
    initialize();
}


/************************/
/*** Public Functions ***/
/************************/

void Enrichment::setParams ()
{
    ParamsIn["Mass"]  = IsosIn.mass;
    ParamsOut["Mass"] = IsosOut.mass;	
}

MassStream Enrichment::doCalc ()
{
    //Does the Enrichmenting
    CompDict incomp  = IsosIn.multByMass();
    CompDict outcomp;
    for (CompIter i = incomp.begin(); i != incomp.end(); i++)
    {
        outcomp[i->first] = (i->second) * sepeff[i->first];
    }
    IsosOut = MassStream (outcomp);
    return IsosOut;
}

MassStream Enrichment::doCalc (CompDict incomp)
{
    //Does the Enrichmenting
    //incomp = input component dictionary of all nuclides. Standard CompDict object. Assigns this to IsosIn.
    IsosIn = MassStream (incomp);
    CompDict outcomp;
    for (CompIter i = incomp.begin(); i != incomp.end(); i++)
    {
        outcomp[i->first] = (i->second) * sepeff[i->first];
    }
    IsosOut = MassStream (outcomp);
    return IsosOut;
}

MassStream Enrichment::doCalc (MassStream instream)
{
    //Does the Enrichmenting
    //instream = input stream of all nuclides. Standard MassStream object.
    IsosIn = instream;
    return doCalc();
}


double Enrichment::PoverF(double x_F, double x_P, double x_W)
{
    //Product over Feed Enrichment Ratio
    return ((x_F - x_W)/(x_P - x_W));
}

double Enrichment::WoverF(x_F, x_P, x_W)
{
    //Waste over Feed Enrichment Ratio
    return ((x_F - x_P)/(x_W - x_P));
}

double Enrichment::alphastar_i(double alpha_0, double M_i, double Mstar)
{
       	return pow(alpha_0, (Mstar - M_i));
}

double Enrichment::Ei(double alpha_0, double M_i, double Mstar)
{
	double astar_i = alphastar_i(alpha0, M_i, Mstar);
    return ((astar_i - 1.0) / (1.0 - pow(astar_i, -N) ));
}

def Si(alpha0, Mi, Mstar):
	astari = alphastari(alpha0, Mi, Mstar)
       	return ((astari - 1.0)/(astari**(M+1) - 1.0))

def FindNM(alpha0, Mstar, compF, j, xjP, xjW, N0, M0):
	#this is the order of exactness.  
	#gives the degree to which N and M are solved for
	ooe = 7.0	

	PoF = PoverF(compF[j], xjP, xjW)
	WoF = WoverF(compF[j], xjP, xjW)
	astarj = alphastari(alpha0, float(j), Mstar)
	N = N0
	M = M0

#	print N, M, "<----"
	
	lhsP = PoF*xjP/compF[j]
	rhsP = (astarj**(M+1.0) - 1.0) / (astarj**(M+1.0) - astarj**(-N))
	lhsW = WoF*xjW/compF[j]
	rhsW = (1.0 - astarj**(-N)) / (astarj**(M+1.0) - astarj**(-N))

	n = 1.0
	
	while 10**(-ooe) <  abs(lhsP - rhsP) and 10**(-ooe) <  abs(lhsW - rhsW):
		if 10**(-ooe) <  abs(lhsW - rhsW):
			M = M - (lhsW - rhsW)*M
			rhsW = (1.0 - astarj**(-N)) / (astarj**(M+1.0) - astarj**(-N))
		if 10**(-ooe) <  abs(lhsP - rhsP):
			N = N + (lhsP - rhsP)*N
			rhsP = (astarj**(M+1.0) - 1.0) / (astarj**(M+1.0) - astarj**(-N))
		if N < 10**(-ooe):
			N = N0 + 1.0*n
			M = M0 + 1.0*n
			n = n + 1.0
			#print "n is now ", n
		if M < 10**(-ooe):
			N = N0 + 1.0*n
			M = M0 + 1.0*n
			n = n + 1.0
			#print "n is now ", n

#		print N, M
#	print N, M, "     ---->"
	return [N, M]
  
def xiP(comp, N, M, alpha0, Mstar, compF, j,  xjP, xjW):
	astari = alphastari(alpha0, float(comp), Mstar)
	return compF[comp]*(astari**(M+1.0) - 1.0) / (astari**(M+1.0) - astari**(-N)) / PoverF(compF[j], xjP, xjW)

def xiW(comp, N, M, alpha0, Mstar, compF, j,  xjP, xjW):
	astari = alphastari(alpha0, float(comp), Mstar)
	return compF[comp]*(1.0 - astari**(-N)) / (astari**(M+1.0) - astari**(-N)) / WoverF(compF[j], xjP, xjW)

#this func takes a given initial guess number of enriching and stripping stages for a given composition
#of fuel with a given jth key component, knowing the values that are desired in both Product and Waste streams.
#Having this it solves for what the actual N and M stage numbers are and also what the product and waste streams 
#compositions are.  It returns precisely these.

def SolveNM(alpha0, Mstar, compF, j,  xjP, xjW, N0, M0):
	NM = FindNM(alpha0, Mstar, compF, j,  xjP, xjW, N0, M0)
#	print "-------------"
#	print "N = ", NM[0]
#	print "M = ", NM[1]

	compP = {}
	compW = {}
	compstr = ""
	compFstr = ""
	compWstr = ""
	compPstr = ""
	for comp in compF.keys():
		compP[comp] = xiP(comp, NM[0], NM[1], alpha0, Mstar, compF, j,  xjP, xjW)
		compW[comp] = xiW(comp, NM[0], NM[1], alpha0, Mstar, compF, j,  xjP, xjW)
		compstr = compstr + comp + "\t"

#	compPstr = makecompstr(compP)
#	compFstr = makecompstr(compF)
#	compWstr = makecompstr(compW)
	
#	print "A Num  " + '\t' + compstr
#	print "Product" + '\t' + compPstr
#	print "Feed   " + '\t' + compFstr
#	print "Waste  " + '\t' + compWstr

	return [NM, compP, compW]


#This func actually solve the whole system of equations.  It uses SolveNM to find the roots for the enriching
# andstripping stage numbers.  It then checks to see if the product and waste streams add up to 100% like they
# should.  If they don't then it trys other values of N and M varied by Newton's Method.

def Comp2UnitySecant(alpha0, Mstar, compF, j,  xjP, xjW, N0, M0):
        #this is the order of exactness.
        #gives the degree to which N and M are solved for
        ooe = 5.0
	#Is the hisorty of N and M that has been input
	HistoryTable =[]
	counter = 0

        N = N0
        M = M0
	NLast = N0 + 1.0
	MLast = M0 + 1.0
	Ntemp = 0.0
	Mtemp = 0.0
	snm = SolveNM(alpha0, Mstar, compF, j,  xjP, xjW, N, M)
	snmLast = SolveNM(alpha0, Mstar, compF, j,  xjP, xjW, NLast, MLast)

	N = snm[0][0]
	M = snm[0][1]
	NLast = snmLast[0][0]
	MLast = snmLast[0][1]

        sumP = 0.0
        sumW = 0.0
	sumPLast = 0.0
	sumWLast = 0.0
        for comp in compF.keys():
	        sumP = sumP + snm[1][comp]
        	sumW = sumW + snm[2][comp]
	        sumPLast = sumPLast + snmLast[1][comp]
        	sumWLast = sumWLast + snmLast[2][comp]

        while 10.0**(-ooe) < abs(1.0 - sumP) and 10.0**(-ooe) < abs(1.0 - sumW):
#		print "---------------------"

                if 10.0**(-ooe) > abs(1.0 - sumP):
                        N = snm[0][0]
                else:
			N = snm[0][0] + (1.0 - sumP)*((snm[0][0] - snmLast[0][0])/(sumP - sumPLast))
			NLast = snm[0][0]
			if N < 0.0:
				N = (snm[0][0] + N0)/2.0
#				print "N < 0"

                if 10.0**(-ooe) > abs(1.0 - sumW):
                        M = snm[0][1]
                else:
			M = snm[0][1] + (1.0 - sumW)*((snm[0][1] - snmLast[0][1])/(sumW - sumWLast))
			MLast = snm[0][1]
			if M < 0.0:
				M = (snm[0][1] + M0)/2.0
#				print "M < 0"

		if InfiniteLoopChckr(HistoryTable, [N,M]):
#			print "Infinite Loop Found and Exception Taken"
			raise Exception('Infinite Loop Found')

		if counter > 10000:
#			print "Secant Counter Hit Limit....Breaking..."
#			break
			raise Exception('Counter Limit Hit')
		else:
			counter = counter + 1

		snmLast = SolveNM(alpha0, Mstar, compF, j,  xjP, xjW, NLast, MLast)
#		print "===="
                snm = SolveNM(alpha0, Mstar, compF, j,  xjP, xjW, N, M)
        	sumP = 0.0
	        sumW = 0.0
		sumPLast = 0.0
		sumWLast = 0.0
        	for comp in compF.keys():
		        sumP = sumP + snm[1][comp]
        		sumW = sumW + snm[2][comp]
		        sumPLast = sumPLast + snmLast[1][comp]
	        	sumWLast = sumWLast + snmLast[2][comp]
	return snm

def Comp2UnityOther(alpha0, Mstar, compF, j,  xjP, xjW, N0, M0):
        #this is the order of exactness.
        #gives the degree to which N and M are solved for
        ooe = 5.0
	#Is the hisorty of N and M that has been input
	HistoryTable =[]

	snm = []
        N = N0
        M = M0
	snm = SolveNM(alpha0, Mstar, compF, j,  xjP, xjW, N, M)
        sumP = 0.0
        sumW = 0.0
	for comp in compF.keys():
		sumP = sumP + snm[1][comp]
		sumW = sumW + snm[2][comp]
        while 10.0**(-ooe) < abs(1.0 - sumP) and 10.0**(-ooe) < abs(1.0 - sumW):
                if 10.0**(-ooe) > abs(1.0 - sumP):
                        N = snm[0][0]
                else:
			N = snm[0][0] - (1.0 - sumP)/(1.0 + sumP)

                if 10.0**(-ooe) > abs(1.0 - sumW):
                        M = snm[0][1]
                else:
			M = snm[0][1] + (1.0 - sumW)/(1.0 + sumW)

		#Note this infinite loop checker does not raise an exception
		#Thus the exception cannot be caught by a try statement and then another
		#root finding Comp2Unity method tried.
		#This simply returns whatever the vaule of the of snm is at the time.
		#This is fine for M* = 235.1 for U, since M* won't be optimized here as an outlier 
		#and since it is looping it is probably signalling around some actual value.
		if InfiniteLoopChckr(HistoryTable, [N,M]):
#			print "Infinite Loop Found and Broken"
			break

                snm = SolveNM(alpha0, Mstar, compF, j,  xjP, xjW, N, M)
                sumP = 0.0
                sumW = 0.0
                for comp in compF.keys():
                        sumP = sumP + snm[1][comp]
                        sumW = sumW + snm[2][comp]
	return snm

#The FindMstar finds a value of Mstar by minimzing the seperative power.  The intial vaules needed are the
#same as for the other functions with the note that Mstar -> Mstar0 where Mstar0 is some intial guess at
#what Mstar might be.

#Comes from eq 32 in Houston Wood paper. Divided by G.
def eq28denom(comp, alpha0, Mstar, j):
	astari = alphastari(alpha0, float(comp), Mstar)
	return math.log(alpha0**(Mstar - float(j)))*((astari - 1.0)/(astari + 1.0))

#Finds the total flow rate over the feed flow rate
#UnityMeth is a string that indicates which Comp2Unity root finding method should be used
#UnityMeth uses the Secant Method by default.
def LoverF(alpha0, Mstar, compF, j,  xjP, xjW, N0, M0, k):
	gotValue = False

	#N,M,Product,Waste
	NMPW = []
	try:
#		print "Secant"
		NMPW = Comp2UnitySecant(alpha0, Mstar, compF, j,  xjP, xjW, N0, M0)
		gotValue = True
	except:
		try:
#			print "Other"
			NMPW = Comp2UnityOther(alpha0, Mstar, compF, j,  xjP, xjW, N0, M0)
			gotValue = True
		except:
			pass

	if gotValue:
		PoF = PoverF(compF[j], xjP, xjW)
		WoF = WoverF(compF[j], xjP, xjW)
		#Matched Ratios
		RP = NMPW[1][j] / NMPW[1][k]
		RF = compF[j] / compF[k]
		RW = NMPW[2][j] / NMPW[2][k]

		LtotalOverF = 0.0
		SWUoverF = 0.0

		for comp in compF.keys():
			tempNumerator = (PoF*NMPW[1][comp]*math.log(RP) + WoF*NMPW[2][comp]*math.log(RW) - compF[comp]*math.log(RF))
#			print tempNumerator
			LtotalOverF = LtotalOverF + tempNumerator / eq28denom(comp, alpha0, Mstar, j)
			SWUoverF = SWUoverF + tempNumerator

		NMPW.append(LtotalOverF)
		#The -1 term is put in because otherwise SWUoverF represents the SWU one would be undoing if you were to
		#deenrich the whole process.  The SWU to enrich is -1x this num.  This is because of the value function used.
		#This is the SWU for 1 kg of feed material
		NMPW.append(-1 * SWUoverF)
		#This is the SWU for 1 kg of Product material
		NMPW.append(-1 * SWUoverF / PoF)
		return NMPW
	else:
		return False

#This is the final function that actually solves for an optimized M* that makes the cascade!
def MstarOptimize(alpha0, Mstar0, compF, j,  xjP, xjW, N0, M0, k):
	#History table that has Mstar, LoF, and slope between this point and the last one 
	hist = []

	#ooe = order of exactness
	#xpn = exponent index
	ooe = 7.0
	xpn = 1.0

	MstarLast = Mstar0
	LoverFLast = LoverF(alpha0, MstarLast, compF, j,  xjP, xjW, N0, M0, k)
	MstarNow = Mstar0 + 0.1
	LoverFNow = LoverF(alpha0, MstarNow, compF, j,  xjP, xjW, N0, M0, k)

	m = slope(MstarNow, LoverFNow[3], MstarLast, LoverFLast[3])
	Sign = m / abs(m)
	del m

	if Sign < 0:
		hist.insert(0, [MstarLast, LoverFLast] )
		hist.insert(0, [MstarNow, LoverFNow] )
	elif 0 < Sign:
		hist.insert(0, [MstarNow, LoverFNow] )
		hist.insert(0, [MstarLast, LoverFLast] )
	else:
		print "Slope for intial conditions is already Zero, so you already have minimized M*"

#	print str(hist[1][0]) + '\t' + str(hist[1][1][3]) + '\t'
#	print str(hist[0][0]) + '\t' + str(hist[0][1][3]) + '\t'

	while xpn < ooe:
		MstarNow = hist[0][0] - Sign*10**(-xpn)

#		print "Got Here"
		LoverFNow = LoverF(alpha0, MstarNow, compF, j,  xjP, xjW, N0, M0, k)

		hist.insert(0, [MstarNow, LoverFNow] )
	
#		print str(hist[0][0]) + '\t' + str(hist[0][1][3]) + '\t',

		if hist[1][1][3] < hist[0][1][3]:
			tempMstar = hist[0][0] - Sign*10**(-xpn)
			tempLoverF = LoverF(alpha0, tempMstar, compF, j,  xjP, xjW, N0, M0, k)
			tempm = slope(MstarNow, LoverFNow[3], tempMstar, tempLoverF[3])
			if tempm == 0.0:
				hist.insert(0, [tempMstar, tempLoverF] )
#				print str(hist[0][0]) + '\t' + str(hist[0][1][3]) + '\t'
				break
			tempSign = tempm / abs(tempm)
			if Sign == tempSign:
				pass
			else:
				xpn = xpn + 1
				hist.pop(0)

				tempMstar = hist[0][0] + Sign*10**(-xpn)
				tempLoverF = LoverF(alpha0, tempMstar, compF, j,  xjP, xjW, N0, M0, k)
				tempm = slope(MstarNow, LoverFNow[3], tempMstar, tempLoverF[3])
				if tempm == 0.0:
					hist.insert(0, [tempMstar, tempLoverF] )
#					print str(hist[0][0]) + '\t' + str(hist[0][1][3]) + '\t'
					break
				Sign = tempm / abs(tempm)
#			print Sign, tempSign,
			del tempMstar, tempLoverF, tempm, tempSign

#		print str(hist[0][0]) + '\t' + str(hist[0][1][3]) + '\t'

	result = hist[0][1]
	result.append(hist[0][0])
	return result
	#Def of result: List with elements,
	# [ [N,M], xP, xW, L/F, M*  ]

#The Following Few Functions yeild the cascade That is optimized and emriches fuel to a certain Max Burnup

#Gives the maximum burnup of fuel with u235 and u236 enrichements and a certain number of batches.
#All variables should be input as floats.   
def MaxBurn(u235, u236, batches):
	a = 3.479397448490798
	b = 11.048652045162228
	c = -0.002046586356345106
	d = 1.6011191339360904
	f = -0.12958835040198202
	g = -2.483552758121467
	h = 1.9153885037054987
	j = 12.797726339337736
	k = 7.975144395645545
	m = 0.741332522502401

	batchcoeff = 2.0*batches/(batches + 1.0)
	part1 = b*((u235 - h)**f)/(u236 - g)
	part2 = -d*(u236 - g)
	part3 = c*((u236 - g)**a)/(u235 - h)
	part4 = j*(u235 - h)**m
	paren = part1 + part2 + part3 + part4 + k

	return batchcoeff*paren

def GetMaxBurn(CasOut, batches):
	return MaxBurn(CasOut['235']*100, CasOut['236']*100, batches)

def CascadeForBurnGoal(alpha0, Mstar0, compF, j,  xjP0, xjW, N0, M0, k, BurnGoal, batches):
	ooe = 4.0

	xjPLast = xjP0
	MaxBurnLast = GetMaxBurn( MstarOptimize(alpha0, Mstar0, compF, j,  xjPLast, xjW, N0, M0, k)[1], batches)
	print xjPLast, MaxBurnLast
	xjPNow = xjP0 + 0.001
	CascadeNow = MstarOptimize(alpha0, Mstar0, compF, j,  xjPNow, xjW, N0, M0, k)
	MaxBurnNow = GetMaxBurn( CascadeNow[1], batches )
	print xjPNow, MaxBurnNow

	while 10**(-ooe) < abs(MaxBurnNow - BurnGoal):
		tempxjP = xjPNow
		tempMaxBurn = MaxBurnNow
	
		xjPNow = xjPNow + (BurnGoal - MaxBurnNow)*(xjPNow - xjPLast)/(MaxBurnNow - MaxBurnLast)
		print xjPNow,
		CascadeNow = MstarOptimize(alpha0, Mstar0, compF, j,  xjPNow, xjW, N0, M0, k)
		MaxBurnNow = GetMaxBurn( CascadeNow[1], batches )
		print MaxBurnNow
	
		xjPLast = tempxjP
		MaxBurnLast = tempMaxBurn

	return CascadeNow

#This set of Function Blends two fuels together but enriches one first, 
#It finds the enrichment needed to get to a certain Burnup
def dictMult(dict, scalar):
	tempdict = {}
	for k in dict.keys():
		tempdict[k] = dict[k] * scalar
	return tempdict

def xFmix(xF1, kg1, xF2, kg2):
	tempxF = {}
	tempxF1 = dictMult(xF1, kg1)
	tempxF2 = dictMult(xF2, kg2)

	for k in tempxF1.keys():
		if k in tempxF2.keys():
			tempxF[k] = tempxF1[k] + tempxF2[k]
		else:
			tempxF[k] = tempxF1[k]
	for k in tempxF2.keys():
		if k in tempxF1.keys():
			pass
		else:
			tempxF[k] = tempxF2[k]

	tempxF = dictMult(tempxF, 1.0/(kg1 + kg2) )
	return tempxF

def CascadeForBurnGoalBlended(alpha0, Mstar0, compF, j,  xjP0, xjW, N0, M0, k, BurnGoal, batches, kg, kgBlend, compFBlend):
	ooe = 5.0

	xjPLast = xjP0
	CascadeLast = MstarOptimize(alpha0, Mstar0, compF, j,  xjPLast, xjW, N0, M0, k)
	xFmixedLast = xFmix( CascadeLast[1], kg, compFBlend, kgBlend)
	MaxBurnLast = GetMaxBurn( xFmixedLast, batches)
	print xjPLast, MaxBurnLast
	xjPNow = xjP0 + 0.001
	CascadeNow = MstarOptimize(alpha0, Mstar0, compF, j,  xjPNow, xjW, N0, M0, k)
	xFmixedNow = xFmix( CascadeNow[1], kg, compFBlend, kgBlend)
	MaxBurnNow = GetMaxBurn( xFmixedNow, batches )
	print xjPNow, MaxBurnNow

	while 10**(-ooe) < abs(MaxBurnNow - BurnGoal):
		tempxjP = xjPNow
	        tempMaxBurn = MaxBurnNow
	
        	xjPNow = xjPNow + (BurnGoal - MaxBurnNow)*(xjPNow - xjPLast)/(MaxBurnNow - MaxBurnLast)
	        print xjPNow,
		CascadeNow = MstarOptimize(alpha0, Mstar0, compF, j,  xjPNow, xjW, N0, M0, k)
		xFmixedNow = xFmix( CascadeNow[1], kg, compFBlend, kgBlend)
		MaxBurnNow = GetMaxBurn( xFmixedNow, batches )
	        print MaxBurnNow

        	xjPLast = tempxjP
	        MaxBurnLast = tempMaxBurn

	return CascadeNow

#This functuion finds the mass of a certain amount of a composition needed to obtain a given burnup when the composition is enriched 
#a set amount.  This differs from the function abopve in which the mass of the composition is known and the enrichment weight percent 
#is solved for to obtain a given burnup.  Initial Conditions should be set close to the value.  
def CascadeForBurnGoalBlendedMass(alpha0, Mstar0, compF, j,  xjP, xjW, N0, M0, k, BurnGoal, batches, kg0, kgBlend, compFBlend):
	ooe = 5.0

	Cascade = MstarOptimize(alpha0, Mstar0, compF, j,  xjP, xjW, N0, M0, k)

	kgLast = kg0 
	xFmixedLast = xFmix( Cascade[1], kgLast, compFBlend, kgBlend)
	MaxBurnLast = GetMaxBurn( xFmixedLast, batches)
	print kgLast, MaxBurnLast
	kgNow = kg0 + 0.01
	xFmixedNow = xFmix( Cascade[1], kgNow, compFBlend, kgBlend)
	MaxBurnNow = GetMaxBurn( xFmixedNow, batches )
	print kgNow, MaxBurnNow

	while 10**(-ooe) < abs(MaxBurnNow - BurnGoal):
		tempkg = kgNow
	        tempMaxBurn = MaxBurnNow
	
        	kgNow = kgNow + (BurnGoal - MaxBurnNow)*(kgNow - kgLast)/(MaxBurnNow - MaxBurnLast)
	        print kgNow,
		xFmixedNow = xFmix( Cascade[1], kgNow, compFBlend, kgBlend)
		MaxBurnNow = GetMaxBurn( xFmixedNow, batches )
	        print MaxBurnNow

        	kgLast = tempkg
	        MaxBurnLast = tempMaxBurn

	Cascade.append(kgNow)
	return Cascade
