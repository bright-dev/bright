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

void print( string s)
{
	cout << s << "\n";
}

int main()
{
	CompDict LEU;
	LEU[922350] = 0.25;
	LEU[922380] = 4.75;
	MassStream msLEU (LEU, -1, "LEU From CompDict");
	msLEU.print_ms();
	print("");

	MassStream msDU ("MassStreamtry01.txt", -1, "DU from File (with char*)");
	msDU.print_ms();
	print("");

	string msf = "MassStreamtry02.txt"; 
	MassStream msSNF (msf, -1, "LWR-SNF from File (with string)");
	msSNF.print_ms();
	print("");

	print("Checking substrem from int set:");
	int intnuc [] = {922380, 94, 2, 57};
	set<int> intset (intnuc, intnuc+4);
	MassStream msIntSub = msSNF.get_sub_stream(intset, "Integer SubSet");
	msIntSub.print_ms();
	print("");

	print("Checking substrem from string set:");
	string strnuc [] = {"922350", "94242", "CS137", "57", "CF", "TQ"};
	set<string> strset (strnuc, strnuc+6);
	MassStream msStrSub = msSNF.get_sub_stream(strset, "String SubSet");
	msStrSub.print_ms();
	print("");

	print("Checking Uranium sub-stream:");
	MassStream msU = msSNF.get_u("Uranium Sub-Stream");
	msU.print_ms();
	print("");

	print("Checking Plutonium sub-stream:");
	MassStream msPU = msSNF.get_pu("Plutonium Sub-Stream");
	msPU.print_ms();
	print("");

	print("Checking Lanthanide sub-stream:");
	MassStream msLAN = msSNF.get_lan("Lanthanide Sub-Stream");
	msLAN.print_ms();
	print("");

	print("Checking Actinide sub-stream:");
	MassStream msACT = msSNF.get_act("Actinide Sub-Stream");
	msACT.print_ms();
	print("");

	print("Checking Transuranic sub-stream:");
	MassStream msTRU = msSNF.get_tru("Transuranic Sub-Stream");
	msTRU.print_ms();
	print("");

	print("Checking Minor Actinides sub-stream:");
	MassStream msMA = msSNF.get_ma("Minor Actinide Sub-Stream");
	msMA.print_ms();
	print("");

	print("Checking Fission Product sub-stream:");
	MassStream msFP = msSNF.get_fp("Fission Product Sub-Stream");
	msFP.print_ms();
	print("");

	print("*************************************************");
	print("*** Now let's try overloading some operators! ***");
	print("*************************************************");
	print("");

	print("Left double Addition:");
	MassStream msLplus  = msLEU + 2.0;
	msLplus.print_ms(); 
	print("");
	
	print("Right double Addition:");
	MassStream msRplus  = 4 + msLEU;
	msRplus.print_ms(); 
	print("");

	print("Right MassStream Addition:");
	MassStream msRMS  = msSNF + msLEU;
	msRMS.print_ms(); 
	print("");

	print("Left MassStream Addition:");
	MassStream msLMS  = msLEU + msSNF;
	msLMS.print_ms(); 
	print("");

	print("Left double Multiplication:");
	MassStream msLmult  = msLEU * 2.0;
	msLmult.print_ms(); 
	print("");

	print("Right double Multiplication:");
	MassStream msRmult  = 4 * msLEU;
	msRmult.print_ms(); 
	print("");

	print("Left double Division:");
	MassStream msLdiv  = msLEU / 2.0;
	msLdiv.print_ms(); 
	print("");

	return 0;
}
