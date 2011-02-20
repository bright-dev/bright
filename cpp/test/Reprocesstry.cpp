// FCComp Test File
//Compile with:
//	g++ -o Reprocesstry Reprocesstry.cpp ../FCComps/Reprocess.cpp ../FCComps/FCComp.cpp ../FCComps/MassStream.cpp ../FCComps/isoname.cpp ../FCComps/genlib.cpp

#include "../Reprocess.h"
#include <map>
#include <iostream>
#include <string>
#include <set>

using namespace std;

void print( string s)
{
	cout << s << "\n";
}

void printSepEff(Reprocess r)
{
	print("\tSeparation Efficiencies:");
	for (SepEffIter i = r.sepeff.begin(); i != r.sepeff.end(); i++)
	{
		print("\t\t" + to_str(i->first) + "\t" + to_str(i->second) );
	}
}

int main()
{
	print("Test Initialize Empty Reprocessing Facility:");
	Reprocess repempty;
	print("\tName\t\t\t" + repempty.name);
	print("\tPass Number\t\t" + to_str(repempty.PassNum));
	print("\tFigure Directory\t" + repempty.figdir);
	print("");

	int itrack_arr [] = {922350, 922360, 922380, 942390, 952421};
	set<int> itrack (itrack_arr, itrack_arr+5); 	

	print("Test Initialize SepEffDict Reprocessing Facility:");
	SepEffDict sed; 
	sed[92] = 0.999;
	sed[942390] = 0.99;
	Reprocess rsed (sed, itrack, "Reprocess_SED");
	print("\tName\t\t\t" + rsed.name);
	print("\tPass Number\t\t" + to_str(rsed.PassNum));
	print("\tFigure Directory\t" + rsed.figdir);
	printSepEff(rsed);
	print("");

	print("Test Initialize String SepEffDict Reprocessing Facility:");
	map<string, double> ssed; 
	ssed["PU"] = 0.99;
	ssed["952421"] = 0.99;
	ssed["92236"] = 0.999;
	ssed["U235"] = 0.999;
	ssed["HYIS"] = 0.99; //should fail, but not break the run
	Reprocess rssed (ssed, itrack, "Reprocess_SSED");
	print("\tName\t\t\t" + rssed.name);
	print("\tPass Number\t\t" + to_str(rssed.PassNum));
	print("\tFigure Directory\t" + rssed.figdir);
	printSepEff(rssed);
	print("");

	print("Test doCalc(CompDict):");
	CompDict cd;
	cd[922350] = 10.0;
	cd[922360] = 1.0;
	cd[922380] = 50.0;
	cd[942390] = 20.0;
	cd[952421] = 5.0;
	rssed.doCalc(cd).print_ms();
	rssed.writeout();
	print("\t\t\tDone!");

	print("Test doCalc(MassStream):");
	rssed.doCalc(rssed.ms_prod).print_ms();
	rssed.writeout();
	print("\t\t\tDone!");

	print("Test doCalc():");
	rssed.ms_feed = rssed.ms_prod;
	rssed.doCalc().print_ms();
	rssed.writeout();
	print("\t\t\tDone!");
	
	return 0;
}
