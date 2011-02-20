// Mass Stream Test File
//Compile with:
//	g++ -o MassStreamtry MassStreamtry.cpp ../FCComps/genlib.cpp ../FCComps/MassStream.cpp ../FCComps/isoname.cpp 

#include "../genlib.h"
#include "../MassStream.h"
#include "../isoname.h"
#include <map>
#include <iostream>
#include <string>
#include <set>

using namespace std;

void Print( string s)
{
	cout << s << "\n";
}

int main()
{
	CompDict LEU;
	LEU[922350] = 0.25;
	LEU[922380] = 4.75;
	MassStream msLEU (LEU, -1, "LEU From CompDict");
	msLEU.Print();
	Print("");

	MassStream msDU ("MassStreamtry01.txt", -1, "DU from File (with char*)");
	msDU.Print();
	Print("");

	string msf = "MassStreamtry02.txt"; 
	MassStream msSNF (msf, -1, "LWR-SNF from File (with string)");
	msSNF.Print();
	Print("");

	Print("Checking substrem from int set:");
	int intnuc [] = {922380, 94, 2, 57};
	set<int> intset (intnuc, intnuc+4);
	MassStream msIntSub = msSNF.getSubStream(intset, "Integer SubSet");
	msIntSub.Print();
	Print("");

	Print("Checking substrem from string set:");
	string strnuc [] = {"922350", "94242", "CS137", "57", "CF", "TQ"};
	set<string> strset (strnuc, strnuc+6);
	MassStream msStrSub = msSNF.getSubStream(strset, "String SubSet");
	msStrSub.Print();
	Print("");

	Print("Checking Uranium sub-stream:");
	MassStream msU = msSNF.getU("Uranium Sub-Stream");
	msU.Print();
	Print("");

	Print("Checking Plutonium sub-stream:");
	MassStream msPU = msSNF.getPU("Plutonium Sub-Stream");
	msPU.Print();
	Print("");

	Print("Checking Lanthanide sub-stream:");
	MassStream msLAN = msSNF.getLAN("Lanthanide Sub-Stream");
	msLAN.Print();
	Print("");

	Print("Checking Actinide sub-stream:");
	MassStream msACT = msSNF.getACT("Actinide Sub-Stream");
	msACT.Print();
	Print("");

	Print("Checking Transuranic sub-stream:");
	MassStream msTRU = msSNF.getTRU("Transuranic Sub-Stream");
	msTRU.Print();
	Print("");

	Print("Checking Minor Actinides sub-stream:");
	MassStream msMA = msSNF.getMA("Minor Actinide Sub-Stream");
	msMA.Print();
	Print("");

	Print("Checking Fission Product sub-stream:");
	MassStream msFP = msSNF.getFP("Fission Product Sub-Stream");
	msFP.Print();
	Print("");

	Print("*************************************************");
	Print("*** Now let's try overloading some operators! ***");
	Print("*************************************************");
	Print("");

	Print("Left double Addition:");
	MassStream msLplus  = msLEU + 2.0;
	msLplus.Print(); 
	Print("");
	
	Print("Right double Addition:");
	MassStream msRplus  = 4 + msLEU;
	msRplus.Print(); 
	Print("");

	Print("Right MassStream Addition:");
	MassStream msRMS  = msSNF + msLEU;
	msRMS.Print(); 
	Print("");

	Print("Left MassStream Addition:");
	MassStream msLMS  = msLEU + msSNF;
	msLMS.Print(); 
	Print("");

	Print("Left double Multiplication:");
	MassStream msLmult  = msLEU * 2.0;
	msLmult.Print(); 
	Print("");

	Print("Right double Multiplication:");
	MassStream msRmult  = 4 * msLEU;
	msRmult.Print(); 
	Print("");

	Print("Left double Division:");
	MassStream msLdiv  = msLEU / 2.0;
	msLdiv.Print(); 
	Print("");

	return 0;
}
