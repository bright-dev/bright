// FCComp Test File
//Compile with:
//	g++ -o Reprocesstry Reprocesstry.cpp ../FCComps/Reprocess.cpp ../FCComps/FCComp.cpp ../FCComps/MassStream.cpp ../FCComps/isoname.cpp ../FCComps/genlib.cpp

#include "../Reprocess.h"
#include <map>
#include <iostream>
#include <string>
#include <set>

using namespace std;

void Print( string s)
{
	cout << s << "\n";
}

void PrintSepEff(Reprocess r)
{
	Print("\tSeparation Efficiencies:");
	for (SepEffIter i = r.sepeff.begin(); i != r.sepeff.end(); i++)
	{
		Print("\t\t" + to_str(i->first) + "\t" + to_str(i->second) );
	}
}

int main()
{
	Print("Test Initialize Empty Reprocessing Facility:");
	Reprocess repempty;
	Print("\tName\t\t\t" + repempty.name);
	Print("\tPass Number\t\t" + to_str(repempty.PassNum));
	Print("\tFigure Directory\t" + repempty.figdir);
	Print("");

	int itrack_arr [] = {922350, 922360, 922380, 942390, 952421};
	set<int> itrack (itrack_arr, itrack_arr+5); 	

	Print("Test Initialize SepEffDict Reprocessing Facility:");
	SepEffDict sed; 
	sed[92] = 0.999;
	sed[942390] = 0.99;
	Reprocess rsed (sed, itrack, "Reprocess_SED");
	Print("\tName\t\t\t" + rsed.name);
	Print("\tPass Number\t\t" + to_str(rsed.PassNum));
	Print("\tFigure Directory\t" + rsed.figdir);
	PrintSepEff(rsed);
	Print("");

	Print("Test Initialize String SepEffDict Reprocessing Facility:");
	map<string, double> ssed; 
	ssed["PU"] = 0.99;
	ssed["952421"] = 0.99;
	ssed["92236"] = 0.999;
	ssed["U235"] = 0.999;
	ssed["HYIS"] = 0.99; //should fail, but not break the run
	Reprocess rssed (ssed, itrack, "Reprocess_SSED");
	Print("\tName\t\t\t" + rssed.name);
	Print("\tPass Number\t\t" + to_str(rssed.PassNum));
	Print("\tFigure Directory\t" + rssed.figdir);
	PrintSepEff(rssed);
	Print("");

	Print("Test doCalc(CompDict):");
	CompDict cd;
	cd[922350] = 10.0;
	cd[922360] = 1.0;
	cd[922380] = 50.0;
	cd[942390] = 20.0;
	cd[952421] = 5.0;
	rssed.doCalc(cd).Print();
	rssed.writeout();
	Print("\t\t\tDone!");

	Print("Test doCalc(MassStream):");
	rssed.doCalc(rssed.IsosOut).Print();
	rssed.writeout();
	Print("\t\t\tDone!");

	Print("Test doCalc():");
	rssed.IsosIn = rssed.IsosOut;
	rssed.doCalc().Print();
	rssed.writeout();
	Print("\t\t\tDone!");
	
	return 0;
}
