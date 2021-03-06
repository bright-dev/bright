// Enrichment Component Class

#include "bright_enrichment.h"



std::string bright::enr_p2t [] = {"MassFeed", "MassProduct", "MassTails", "N", "M", \
                                  "Mstar", "TotalPerFeed", "SWUperFeed", "SWUperProduct"};
std::set<std::string> bright::enr_p2track (enr_p2t, enr_p2t+9);

/************************************************/
/*** Enrichment Component Class and Functions ***/
/************************************************/

/***************************/
/*** Protected Functions ***/
/***************************/

void bright::Enrichment::initialize(EnrichmentParameters ep)
{
  // Initializes the enrichment component. 
  alpha_0 = ep.alpha_0;
  Mstar_0 = ep.Mstar_0;
  j       = ep.j;
  k       = ep.k;
  N0      = ep.N0;
  M0      = ep.M0;
  xP_j    = ep.xP_j;
  xW_j    = ep.xW_j;
};


/*******************************/
/*** Enrichment Constructors ***/
/*******************************/

bright::Enrichment::Enrichment(std::string n) : bright::FCComp(enr_p2track, n)
{
  // Enrichmenting Fuel Cycle Component.  Applies Separation Efficiencies.
  initialize(bright::UraniumEnrichmentDefaults);
}


bright::Enrichment::Enrichment(bright::EnrichmentParameters ep, std::string n) : bright::FCComp(enr_p2track, n)
{
  // Enrichmenting Fuel Cycle Component.  Applies Separation Efficiencies.
  initialize(ep);
};

bright::Enrichment::~Enrichment()
{
}


/************************/
/*** Public Functions ***/
/************************/

void bright::Enrichment::calc_params()
{
  params_prior_calc["MassFeed"]  = mat_feed.mass;
  params_after_calc["MassFeed"] = 0.0;	

  params_prior_calc["MassProduct"]  = 0.0;
  params_after_calc["MassProduct"] = mat_prod.mass;	

  params_prior_calc["MassTails"]  = 0.0;
  params_after_calc["MassTails"] = mat_tail.mass;	

  params_prior_calc["N"]  = N;
  params_after_calc["N"] = N;	

  params_prior_calc["M"]  = M;
  params_after_calc["M"] = M;	

  params_prior_calc["Mstar"]  = Mstar;
  params_after_calc["Mstar"] = Mstar;	

  params_prior_calc["TotalPerFeed"]  = TotalPerFeed;
  params_after_calc["TotalPerFeed"] = TotalPerFeed;	

  params_prior_calc["SWUperFeed"]  = SWUperFeed;
  params_after_calc["SWUperFeed"] = 0.0;	

  params_prior_calc["SWUperProduct"]  = 0.0;
  params_after_calc["SWUperProduct"] = SWUperProduct;	
};


pyne::Material bright::Enrichment::calc()
{
  // Does the Enriching
  MstarOptimize();
  return mat_prod;
};


pyne::Material bright::Enrichment::calc(pyne::comp_map incomp)
{
  // Does the Enriching
  // incomp = input component dictionary of all nuclides. Standard pyne::comp_map object. Assigns this to mat_feed.
  mat_feed = pyne::Material (incomp);
  return calc();
}

pyne::Material bright::Enrichment::calc (pyne::Material mat)
{
  // Does the Enrichmenting
  // mat = input stream of all nuclides. Standard pyne::Material object.
  mat_feed = mat;
  return calc();
}


double bright::Enrichment::PoverF(double x_F, double x_P, double x_W)
{
  // Product over Feed Enrichment Ratio
  return ((x_F - x_W)/(x_P - x_W));
}

double bright::Enrichment::WoverF(double x_F, double x_P, double x_W)
{
  // Waste over Feed Enrichment Ratio
  return ((x_F - x_P)/(x_W - x_P));
}


double bright::Enrichment::get_alphastar_i (double M_i)
{
  // M_i is the mass of the ith isotope    
  return pow(alpha_0, (Mstar - M_i));
}

double bright::Enrichment::get_Ei (double M_i)
{
  double alphastar_i = get_alphastar_i(M_i);
  return ((alphastar_i - 1.0) / (1.0 - pow(alphastar_i, -N) ));
};


double bright::Enrichment::get_Si (double M_i)
{
  double alphastar_i = get_alphastar_i(M_i);
  return ((alphastar_i - 1.0)/(pow(alphastar_i, M+1) - 1.0));
};

void bright::Enrichment::FindNM()
{
  // This give the order-of-exactness to which N and M are solved for.
  double ooe = 7.0;
  double tolerance = pow(10.0, -ooe);

  double PoF = PoverF(mat_feed.comp[j], xP_j, xW_j);
  double WoF = WoverF(mat_feed.comp[j], xP_j, xW_j);
  double alphastar_j = get_alphastar_i(pyne::atomic_mass(j));

  // Save original state of N & M
  double origN = N;
  double origM = M;

  if (2 < bright::verbosity)
    std::cout << "    <---- N = " << N << "\tM = " << M << "\n";     

  double lhsP = PoF * xP_j / mat_feed.comp[j];
  double rhsP = (pow(alphastar_j, M+1.0) - 1.0) / (pow(alphastar_j, M+1.0) - pow(alphastar_j, -N));
  double lhsW = WoF * xW_j / mat_feed.comp[j];
  double rhsW = (1.0 - pow(alphastar_j, -N)) / (pow(alphastar_j, M+1.0) - pow(alphastar_j, -N));

  double n = 1.0;
  while (tolerance < fabs(lhsP - rhsP) && tolerance < fabs(lhsW - rhsW))
  {
    if (tolerance < fabs(lhsP - rhsP))
    {
      N = N - ((lhsP - rhsP) * N);
      rhsP = (pow(alphastar_j, M+1.0) - 1.0) / (pow(alphastar_j, M+1.0) - pow(alphastar_j, -N));
    };

    if (tolerance < fabs(lhsW - rhsW))
    {
      M = M - ((lhsW - rhsW) * M);
      rhsW = (1.0 - pow(alphastar_j, -N)) / (pow(alphastar_j, M+1.0) - pow(alphastar_j, -N));
    };

    // print summary
    if (3 < bright::verbosity)
    {
      std::cout << "            N = " << N << "\tlhsP = " << lhsP << "\trhsP = " << rhsP << "\n";
      std::cout << "            M = " << M << "\tlhsW = " << lhsW << "\trhsW = " << rhsW << "\n";
    };

    if (N < tolerance)
    {
      N = origN + n;
      M = origM + n;
      n = n + 1.0;

      if (2 < bright::verbosity)
        std::cout << "          N set n equal to " << n << "\n";
    };

    if (M < tolerance)
    {
      N = origN + n;
      M = origM + n;
      n = n + 1.0;

      if (2 < bright::verbosity)
        std::cout << "          M set n equal to " << n << "\n";
    };

    if (2 < bright::verbosity)
      std::cout << "    ----- N = " << N << "\tM = " << M << "\n";     
  };

  if (2 < bright::verbosity)
    std::cout << "    ----> N = " << N << "\tM = " << M << "\n";     
  return; 
};
  

double bright::Enrichment::xP_i(int i)
{
  double alphastar_i = get_alphastar_i(pyne::atomic_mass(i));
  double numerator = mat_feed.comp[i]*(pow(alphastar_i, M+1.0) - 1.0);
  double denominator = (pow(alphastar_i, M+1.0) - pow(alphastar_i, -N)) / PoverF(mat_feed.comp[j], xP_j, xW_j);
  return numerator / denominator;
};


double bright::Enrichment::xW_i(int i)
{
  double alphastar_i = get_alphastar_i(pyne::atomic_mass(i));
  double numerator = mat_feed.comp[i] * (1.0 - pow(alphastar_i, -N));
	double denominator = (pow(alphastar_i, M+1.0) - pow(alphastar_i, -N)) / WoverF(mat_feed.comp[j], xP_j, xW_j);
  return numerator / denominator;
};


void bright::Enrichment::SolveNM()
{
  //This function takes a given initial guess number of enriching and stripping stages 
  //for a given composition of fuel with a given jth key component, knowing the values 
  //that are desired in both Product and Waste streams.  Having this it solves for what 
  //the actual N and M stage numbers are and also what the product and waste streams 
  //compositions are.  It returns precisely these.

  FindNM();

  pyne::comp_map compP;
  pyne::comp_map compW;

  for (pyne::comp_iter i = mat_feed.comp.begin(); i != mat_feed.comp.end(); i++)
  {
    compP[i->first] = xP_i(i->first);
    compW[i->first] = xW_i(i->first);
  };

  mat_prod  = pyne::Material(compP);
  mat_tail = pyne::Material(compW);

  return;
};


void bright::Enrichment::Comp2UnitySecant()
{
  // This function actually solves the whole system of equations.  It uses SolveNM 
  // to find the roots for the enriching and stripping stage numbers.  It then 
  // checks to see if the product and waste streams meet their target enrichments
  // for the jth component like they should.  If they don't then it trys other values 
  // of N and M varied by Newton's Method.  Rinse and repeat as needed.

  // This give the order-of-exactness to which N and M are solved for.
  double ooe = 7.0;
  double tolerance = pow(10.0, -ooe);

  // Is the hisorty of N and M that has been input
  std::vector<double> historyN;
  std::vector<double> historyM;

  // Start iteration Counter
  int counter = 0;

  // Set first two points
  double lastN = N0 + 1.0;
  double lastM = M0 + 1.0;
  double currN = N0;
  double currM = M0;

  // Initialize 'last' point
  N = lastN;
  M = lastM;
  SolveNM();
  double lastxP_j = mat_prod.comp[j];
  double lastxW_j = mat_tail.comp[j];
  historyN.push_back(N);
  historyN.push_back(M);

  // Initialize 'current' point
  N = currN;
  M = currM;
  SolveNM();
  double currxP_j = mat_prod.comp[j];
  double currxW_j = mat_tail.comp[j];
  historyN.push_back(N);
  historyN.push_back(M);

  // My guess is that what we are checkin here is that the isotopic compositions
  // make sense with abs(1.0 - massCurrP) rather than calculatign the 
  // relative product to watse mass streams.
  double tempCurrN = 0.0;
  double tempCurrM = 0.0;
  double tempLastN = 0.0;
  double tempLastM = 0.0;

  while (tolerance < fabs(xP_j - currxP_j) || tolerance < fabs(xW_j - currxW_j))
  {
    if (1 < bright::verbosity)
      std::cout << "--------------------\n";

    if (tolerance <= fabs(xP_j - currxP_j))
    {
      // Make a new guess for N
      tempCurrN = currN;
      tempLastN = lastN;
      currN = currN + (xP_j - currxP_j)*((currN - lastN)/(currxP_j - lastxP_j));
      lastN = tempCurrN;

      // If the new value of N is less than zero, reset.
      if (currN < 0.0)
      {
        currN = (tempCurrN + tempLastN)/2.0;
        if (1 < bright::verbosity)
          std::cout << "    N < 0, resetting.\n";
      };
    };

    if (tolerance <= fabs(xW_j - currxW_j))
    {
      // Make a new guess for M
      tempCurrM = currM;
      tempLastM = lastM;
      currM = currM + (xW_j - currxW_j)*((currM - lastM)/(currxW_j - lastxW_j));
      lastM = tempCurrM;

      // If the new value of M is less than zero, reset.
      if (M < 0.0)
      {
        currM = (tempCurrM + tempLastM)/2.0;
        if (1 < bright::verbosity)
          std::cout << "    M < 0, resetting.\n";
      };
    };

    // Check for infinite loops
    for (int h = 0; h < historyN.size(); h++)
    {
      if (historyN[h] == currN && historyM[h] == currM)
      {
        if (1 < bright::verbosity)
          std::cout << "~~~ Infinite loop found and exception thrown! ~~~.\n";
        throw EnrichmentInfiniteLoopError();
      };
    };

    if (150 <= historyN.size())
    {
      historyN.erase(historyN.begin());
      historyM.erase(historyM.begin());
    };
    historyN.push_back(N);
    historyM.push_back(M);

    if (10000 < counter)
    {
      if (1 < bright::verbosity)
        std::cout << "~~~ Secant method counter limit hit! ~~~.\n";
      throw EnrichmentIterationLimit();
    }
    else
    {
      counter = counter + 1;
    };

    // Calculate new isotopics for valid (N, M)        
    lastxP_j = currxP_j;
    lastxW_j = currxW_j;

    N = currN;
    M = currM;
    SolveNM();
    currxP_j = mat_prod.comp[j];
    currxW_j = mat_tail.comp[j];

    if (1 < bright::verbosity)
    {
      std::cout << "Product Mass: " << currxP_j << "\tWaste Mass: " << currxW_j << "\n";
      std::cout << "====================\n";
    };
  };

  return;
};


// I have serious doubts that this works...
void bright::Enrichment::Comp2UnityOther()
{
  // This give the order-of-exactness to which N and M are solved for.
  double ooe = 5.0;
  double tolerance = pow(10.0, -ooe);

  // Is the hisorty of N and M that has been input
  std::vector<double> historyN;
  std::vector<double> historyM;

  // Initial point
  N = N0;
  M = M0;
  SolveNM();
  double massP = mat_prod.mass;
  double massW = mat_tail.mass;

  while (tolerance < fabs(1.0 - massP) && tolerance < fabs(1.0 - massW) )
  {
    if (tolerance <= fabs(1.0 - massP))
      N = N - (1.0 - massP)/(1.0 + massP);

    if (tolerance <= fabs(1.0 - massW))
      M = M + (1.0 - massW)/(1.0 + massW);

    // Note this infinite loop checker does not raise an exception
    // Thus the exception cannot be caught by a try statement and then another
    // root finding Comp2Unity method tried.
    // This simply leaves N and M in whatever their current state is at the time.
    // This is fine for M* = 235.1 for U, since M* won't be optimized here as an outlier 
    // and since it is looping it is probably signalling around some actual value.

    // Check for infinite loops
    for (int h = 0; h < historyN.size(); h++)
    {
      if (historyN[h] == N && historyM[h] == M)
      {
        //throw EnrichmentInfiniteLoopError(); //Possible future use.
        return;
      };

      if (150 <= historyN.size())
      {
        historyN.erase(historyN.begin());
        historyM.erase(historyM.begin());
      };
    };
    historyN.push_back(N);
    historyM.push_back(M);

    // Calculate new masses
    SolveNM();
    massP = mat_prod.mass;
    massW = mat_tail.mass;
  };

  return;
};


double bright::Enrichment::deltaU_i_OverG(int i)
{
  // Solves for a stage separative power relevant to the ith component
  // per unit of flow G.  This is taken from Equation 31 divided by G 
  // from the paper "Wood, Houston G., Borisevich, V. D. and Sulaberidze, G. A.,
  // 'On a Criterion Efficiency for Multi-Isotope Mixtures Separation', 
  // Separation Science and Technology, 34:3, 343 - 357"
  // To link to this article: DOI: 10.1081/SS-100100654
  // URL: http://dx.doi.org/10.1081/SS-100100654

  double alphastar_i = get_alphastar_i(pyne::atomic_mass(i));
  return log(pow( alpha_0, (Mstar - pyne::atomic_mass(j)) )) * ((alphastar_i - 1.0)/(alphastar_i + 1.0));
};


void bright::Enrichment::LoverF()
{
  // This function finds the total flow rate (L) over the feed flow rate (F)
  bool compConverged = false;

  try
  {
    // Try secant method first
    if (0 < bright::verbosity)
      std::cout << "Attempting Secant Method in L/F Calculation...\n";
    Comp2UnitySecant();
    compConverged = true;
  }
  catch (...)
  {
    try
    {
      // Then try other cr8zy method
      if (0 < bright::verbosity)
        std::cout << "Attempting Another Method in L/F Calculation...\n";
      Comp2UnityOther();
    	compConverged = true;
    }
		catch (...)
    {
      // No other methods to try!
      compConverged = false;
    };
  };

  if (compConverged)
  {
    double PoF = PoverF(mat_feed.comp[j], xP_j, xW_j);
    double WoF = WoverF(mat_feed.comp[j], xP_j, xW_j);

    // Matched Flow Ratios
    double RF = mat_feed.comp[j]   / mat_feed.comp[k];
    double RP = mat_prod.comp[j]  / mat_prod.comp[k];
    double RW = mat_tail.comp[j] / mat_tail.comp[k];

    double LtotalOverF = 0.0;
    double SWUoverF = 0.0;
    double tempNumerator = 0.0; 

    for (pyne::comp_iter i = mat_feed.comp.begin(); i != mat_feed.comp.end(); i++)
    {
      tempNumerator = (PoF*mat_prod.comp[i->first]*log(RP) + WoF*mat_tail.comp[i->first]*log(RW) - mat_feed.comp[i->first]*log(RF));
      LtotalOverF = LtotalOverF + (tempNumerator / deltaU_i_OverG(i->first));
      SWUoverF = SWUoverF + tempNumerator;
    };

    if (0 < bright::verbosity)
      std::cout << "    L/F = " << LtotalOverF << "\n";        

    // Assign flow rates
    TotalPerFeed = LtotalOverF;

    // The -1 term is put in the SWU calculation because otherwise SWUoverF   
    // represents the SWU that would be undone if you were to deenrich the 
    // whole process.  Thus the SWU to enrich is -1x this number.  This is 
    // a by-product of the value function used as a constraint.
    SWUperFeed    = -1 * SWUoverF;          //This is the SWU for 1 kg of Feed material.
    SWUperProduct = -1 * SWUoverF / PoF;	//This is the SWU for 1 kg of Product material.

    // Assign Isotopic streams the proper masses.
    mat_prod.mass  = mat_feed.mass * PoF;
    mat_tail.mass = mat_feed.mass * WoF;
  };

  return;
};


void bright::Enrichment::MstarOptimize()
{
  // The MstarOptimize function finds a value of Mstar by minimzing the seperative power.  
  // Note that Mstar0 represents an intial guess at what Mstar might be.
  // This is the final function that actually solves for an optimized M* that makes the cascade!

  // History table that has Mstar, LoF, and slope between this point and the last one 
  // hist = []

  // This give the order-of-exactness to which M* is solved for.
  double ooe = 7.0;
  double tolerance = pow(10.0, -ooe);

  // xpn is the exponential index index
  double xpn = 1.0;

  // Initialize 'last' point
  double lastMstar = Mstar_0;
  Mstar = lastMstar;
  LoverF();
  double lastLoverF = TotalPerFeed;

  // Initialize 'current' point
  double currMstar = Mstar_0 + 0.1;
  Mstar = currMstar;
  LoverF();
  double currLoverF = TotalPerFeed;

  double m = pyne::slope(currMstar, currLoverF, lastMstar, lastLoverF);
  double m_sign = m / fabs(m);

  double tempMstar;
  double tempLoverF;
  double tempm;
  double tempm_sign;

  if (0.0 < m_sign)
  {
    tempMstar  = lastMstar;
    tempLoverF = lastLoverF;
    lastMstar  = currMstar;
    lastLoverF = currLoverF;
    currMstar  = tempMstar;
    currLoverF = tempLoverF;
  };

  // print points
  if (0 < bright::verbosity)
  {
    std::cout << "Last Point: M* = " << lastMstar << "\tL/F = " << lastLoverF << "\n";
    std::cout << "Curr Point: M* = " << currMstar << "\tL/F = " << currLoverF << "\n";
  };

  // Start iterations.    
  while (xpn < ooe)
  {
    // Check that parameters are still well-formed
    if ( isnan(currMstar) || isnan(currLoverF) || isnan(lastMstar) || isnan(lastLoverF) )
      throw EnrichmentIterationNaN();

    lastMstar  = currMstar;
    lastLoverF = currLoverF;

    currMstar = currMstar - (m_sign * pow(10.0, -xpn));
    Mstar = currMstar;
    LoverF();
    currLoverF = TotalPerFeed;

    // print Point
    if (0 < bright::verbosity)
      std::cout << "Next Point: M* = " << currMstar << "\tL/F = " << currLoverF << "\n";

    if (lastLoverF < currLoverF)
    {
      tempMstar = currMstar - (m_sign * pow(10.0, -xpn));
      Mstar = tempMstar;
      LoverF();
      tempLoverF = TotalPerFeed; 

      tempm = pyne::slope(currMstar, currLoverF, tempMstar, tempLoverF);
      if (tempm == 0.0)
      {
        lastMstar  = currMstar;
        lastLoverF = currLoverF;
        currMstar  = tempMstar;
        currLoverF = tempLoverF;

        // print Point
        if (0 < bright::verbosity)
          std::cout << "Next Point: M* = " << currMstar << "\tL/F = " << currLoverF << "\n";
        break;
      };

      tempm_sign = tempm / fabs(tempm);
      if (m_sign != tempm_sign)
      {
        xpn = xpn + 1;

        tempMstar = lastMstar + (m_sign * pow(10.0, -xpn));
        Mstar = tempMstar;
        LoverF();
        tempLoverF = TotalPerFeed;
        tempm = pyne::slope(lastMstar, lastLoverF, tempMstar, tempLoverF);

        if (tempm == 0.0)
        {
          currMstar  = tempMstar;
          currLoverF = tempLoverF;

          // print Point
          if (0 < bright::verbosity)
            std::cout << "Next Point: M* = " << currMstar << "\tL/F = " << currLoverF << "\n";
          break;
        };

        m_sign = tempm / fabs(tempm);
      };
    };
  };

  Mstar        = currMstar;
  TotalPerFeed = currLoverF;
  return;
};
