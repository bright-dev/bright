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

double Enrichment::get_alphastar_i (double M_i, double Mstar)
{
    //M_i is the mass of the ith isotope    
    return pow(alpha_0, (Mstar - M_i));
}

double Enrichment::get_Ei (double M_i, double Mstar)
{
    double alphastar_i = get_alphastar_i(M_i, Mstar);
    return ((alphastar_i - 1.0) / (1.0 - pow(alphastar_i, -N) ));
};

double Enrichment::get_Si (double Mi, double Mstar)
{
    double alphastar_i = get_alphastar_i(M_i, Mstar);
    return ((alphastari - 1.0)/(pow(alphastari, M+1) - 1.0));
};

void Enrichment::FindNM(double Mstar)
{
    //This give the order-of-exactness to which N and M are solved for.
    double ooe = 7.0;
    double tolerance = pow(10.0, -ooe);

    double PoF = PoverF(IsosIn[j], xP_j, xW_j);
    double WoF = WoverF(IsosIn[j], xP_j, xW_j);
    alphastar_j = get_alphastar_i(isoname::nuc_weight(j), Mstar);
    N = N0;
    M = M0;

    double lhsP = PoF * xP_j / IsosIn[j];
    double rhsP = (pow(alphastar_j, M+1.0) - 1.0) / (pow(alphastar_j, M+1.0) - pow(alphastar_j, -N));
    double lhsW = WoF * xW_j / IsosIn[j];
    double rhsW = (1.0 - pow(alphastar_j, -N)) / (pow(alphastar_j, M+1.0) - pow(alphastar_j, -N))

    double n = 1.0;
	while (tolerance < fabs(lhsP - rhsP) && tolerance < fabs(lhsW - rhsW))
    {
		if (tolerance < fabs(lhsW - rhsW))
        {
			M = M - ((lhsW - rhsW) * M);
			rhsW = (1.0 - pow(alphastar_j, -N)) / (pow(alphastar_j, M+1.0) - pow(alphastar_j, -N));
        };

		if (tolerance < abs(lhsP - rhsP))
        {
			N = N + ((lhsP - rhsP) * N);
			rhsP = (pow(alphastar_j, M+1.0) - 1.0) / (pow(alphastar_j, M+1.0) - pow(alphastar_j, -N));
        };

		if (N < tolerance)
        {
			N = N0 + (1.0 * n);
			M = M0 + (1.0 * n);
			n = n + 1.0;
        };

		if (M < tolerance)
			N = N0 + (1.0 * n);
			M = M0 + (1.0 * n);
			n = n + 1.0;
    };
};
  
double Enrichment::xP_i(int i, double Mstar)
{
    double alphastar_i = get_alphastar_i(isoname::nuc_weight(i), Mstar);
    double numerator = IsosIn[i]*(pow(alphastar_i, M+1.0) - 1.0);
    double denominator = (pow(alphastar_i, M+1.0) - pow(alphastar_i, -N)) / PoverF(IsosIn[j], xP_j, xW_j);
    return numerator / denominator;
};

double Enrichment::xW_i(int i, double Mstar)
{
    double alphastar_i = get_alphastar_i(isoname::nuc_weight(i), Mstar);
    double numerator = IsosIn[i] * (1.0 - pow(alphastar_i, -N));
	double denominator = (pow(alphastar_i, M+1.0) - pow(alphastar_i, -N)) / WoverF(IsosIn[j], xP_j, xW_j);
    return numerator / denominator;
};


def Enrichment::SolveNM(double Mstar)
{
    //This function takes a given initial guess number of enriching and stripping stages 
    //for a given composition of fuel with a given jth key component, knowing the values 
    //that are desired in both Product and Waste streams.  Having this it solves for what 
    //the actual N and M stage numbers are and also what the product and waste streams 
    //compositions are.  It returns precisely these.

    FindNM(Mstar);

    CompDict compP;
    CompDict compW;

    for (CompIter i = IsosIn.comp.begin(); i != IsosIn.comp.end(); i++)
    {
        compP[i->first] = xP_i(i->first, Mstar);
        compW[i->first] = xW_i(i->first, Mstar);
    };

    IsosOut  = MassStream(compP);
    IsosTail = MassStream(compW);

    return;
};


void Enrichment::Comp2UnitySecant(double Mstar)
{
    //This function actually solves the whole system of equations.  It uses SolveNM 
    //to find the roots for the enriching and stripping stage numbers.  It then 
    //checks to see if the product and waste streams add up to 100% like they
    //should.  If they don't then it trys other values of N and M varied by Newton's Method.

    //This give the order-of-exactness to which N and M are solved for.
    double ooe = 5.0;
    double tolerance = pow(10.0, -ooe);

	//Is the hisorty of N and M that has been input
	std::list<double> historyN;
	std::list<double> historyM;

    //Start iteration Counter
	int counter = 0;

    //Set first two points
	double lastN = N0 + 1.0;
	double lastM = M0 + 1.0;
    double currN = N0;
    double currM = M0;

    //Initialize 'last' point
    N = lastN;
    M = lastM;
    SolveNM(Mstar);
    double massLastP = IsosOut.mass;
    double massLastW = IsosTail.mass;
    histortN.push_back(N);
    histortN.push_back(M);

    //Initialize 'current' point
    N = currN;
    M = currM;
    SolveNM(Mstar);
    double massCurrP = IsosOut.mass;
    double massCurrW = IsosTail.mass;
    histortN.push_back(N);
    histortN.push_back(M);

    //My guess is that what we are checkin here is that the isotopic compositions
    //make sense with abs(1.0 - massCurrP) rather than calculatign the 
    //relative product to watse mass streams.
    double tempN = 0.0;
    double tempM = 0.0;
    while (tolerance < fabs(1.0 - massCurrP) && tolerance < fabs(1.0 - massCurrW))
    {
        if (tolerance <= fabs(1.0 - massCurrP)
        {
            //Make a new guess for N
            tempN = currN;
            currN = currN + (1.0 - massCurrP)*((currN - lastN)/(massCurrP - massLastP));
			lastN = tempN;

            //If the new value of N is less than zero, reset.
			if (currN < 0.0)
            {
				currN = (tempN + N0)/2.0;
            };
        };

        if (tolerance <= fabs(1.0 - massCurrW)
        {
            //Make a new guess for M
            tempM = currM;
            currM = currM + (1.0 - massCurrW)*((currM - lastM)/(massCurrW - massLastW));
            lastM = tempM;

            //If the new value of M is less than zero, reset.
            if (M < 0.0)
            {
                currM = (tempM + M0)/2.0;
            };
        };

        //Check for infinite loops
        for (int h = 0; h < historyN.size(), h++)
        {
            if (historyN[h] == currN && historyM[h] == currM)
            {
                throw EnrichmentInfiniteLoopError();
            };

            if (150 <= historyN.size())
            {
                historyN.pop_front();
                historyM.pop_front();
            };
        };
        historyN.push_back(N);
        historyM.push_back(M);

        if (10000 < counter)
        {
            throw EnrichmentIterationLimit();
        }
        else
        {
            counter = counter + 1;
        };

        //Calculate new isotopics for valid (N, M)        
        massLastP = massCurrP;
        massLastW = massCurrW;

        N = currN;
        M = currM;
        SolveNM(Mstar);
        massCurrP = IsosOut.mass;
        massCurrW = IsosTail.mass;
    };

	return;
};

void Enrichment::Comp2UnityOther(double Mstar)
{
    //This give the order-of-exactness to which N and M are solved for.
    double ooe = 5.0;
    double tolerance = pow(10.0, -ooe);

	//Is the hisorty of N and M that has been input
	std::list<double> historyN;
	std::list<double> historyM;

    //Initial point
    N = N0;
    M = M0;
	SolveNM(Mstar);
    double massP = IsosOut.mass;
    double massW = IsosTails.mass;

    while (tolerance < fabs(1.0 - massP) && tolerance < fabs(1.0 - massW) )
    {
        if (tolerance <= fabs(1.0 - massP))
        {
            N = N - (1.0 - massP)/(1.0 + massP);
        };

        if (tolerance <= fabs(1.0 - massW))
        {
			M = M + (1.0 - massW)/(1.0 + massW)
        };

		//Note this infinite loop checker does not raise an exception
		//Thus the exception cannot be caught by a try statement and then another
		//root finding Comp2Unity method tried.
		//This simply leaves N and M in whatever their current state is at the time.
		//This is fine for M* = 235.1 for U, since M* won't be optimized here as an outlier 
		//and since it is looping it is probably signalling around some actual value.

        //Check for infinite loops
        for (int h = 0; h < historyN.size(), h++)
        {
            if (historyN[h] == currN && historyM[h] == currM)
            {
                //throw EnrichmentInfiniteLoopError(); //Possible future use.
                return;
            };

            if (150 <= historyN.size())
            {
                historyN.pop_front();
                historyM.pop_front();
            };
        };
        historyN.push_back(N);
        historyM.push_back(M);

        //Calculate new masses
        SolveNM(Mstar);
        massP = IsosOut.mass;
        massW = IsosTail.mass;
    };

	return;
};

double Enrichment::deltaU_i_OverG(int i, double Mstar)
{
    //Solves for a stage separative power relevant to the ith component
    //per unit of flow G.  This is taken from Equation 31 divided by G 
    //from the paper "Wood, Houston G., Borisevich, V. D. and Sulaberidze, G. A.,
    //'On a Criterion Efficiency for Multi-Isotope Mixtures Separation', 
    //Separation Science and Technology, 34:3, 343 - 357"
    //To link to this article: DOI: 10.1081/SS-100100654
    //URL: http://dx.doi.org/10.1081/SS-100100654

    double alphastar_i = get_alphastar_i(isoname::nuc_weight(i), Mstar);
	return log(pow(alpha_0, (Mstar - isoname::nuc_weight(j))) * ((alphastar_i - 1.0)/(alphastar_i + 1.0));
};

void Enrichment::LoverF(double Mstar)
{
    //This function finds the total flow rate (L) over the feed flow rate (F)

	bool compConverged = false;

	try
    {
        //Try secant method first
		Comp2UnitySecant(Mstar);
		compConverged = true;
    }
	catch (...)
    {
		try
        {
            //Then try other cr8zy method
			Comp2UnityOther(Mstar);
    		compConverged = true;
        }
		catch (...)
        {
            //Nol other methods to try!
    		compConverged = false;
        };
    };

	if (compConverged)
    {
		double PoF = PoverF(IsosIn[j], xP_j, xW_j);
		double WoF = WoverF(IsosIn[j], xP_j, xW_j);

		//Matched Flow Ratios
		double RF = IsosIn.comp[j]   / IsosIn.comp[k];
		double RP = IsosOut.comp[j]  / IsosOut.comp[k];
		double RW = IsosTail.comp[j] / IsosTail.comp[k];

		double LtotalOverF = 0.0;
		double SWUoverF = 0.0;
        double tempNumerator = 0.0; 

		for (CompIter i = IsosIn.comp.begin(); i != IsosIn.comp.end(); i++)
        {
			tempNumerator = (PoF*IsosOut.comp[i->first]*log(RP) + WoF*IsosTail.comp[i->first]*log(RW) - IsosIn.comp[i->first]*log(RF));
			LtotalOverF = LtotalOverF + (tempNumerator / deltaU_i_OverG(i->first, Mstar));
			SWUoverF = SWUoverF + tempNumerator;
        };

        //Assign flow rates
        TotalPerFeed = LtotalOverF;

        //The -1 term is put in the SWU calculation because otherwise SWUoverF 
        //represents the SWU that would be undone if you were to deenrich the 
        //whole process.  Thus the SWU to enrich is -1x this number.  This is 
        //a by-product of the value function used as a constraint.
        SWUperFeed    = -1 * SWUoverF;          //This is the SWU for 1 kg of Feed material.
		SWUperProduct = -1 * SWUoverF / PoF;	//This is the SWU for 1 kg of Product material.

        //Assign Isotopic streams the proper masses.
        IsosOut.mass  = IsosIn.mass * PoF;
        IsosTail.mass = IsosIn.mass * WoF;
    };

    return;
};

#The FindMstar finds a value of Mstar by minimzing the seperative power.  The intial vaules needed are the
#same as for the other functions with the note that Mstar -> Mstar0 where Mstar0 is some intial guess at
#what Mstar might be.

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

