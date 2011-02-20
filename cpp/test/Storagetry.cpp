// FCComp Test File
//Compile with:
//	h5c++ -o Storagetry Storagetry.cpp ../FCComps/Storage.cpp ../FCComps/FCComp.cpp ../FCComps/MassStream.cpp ../FCComps/isoname.cpp ../FCComps/genlib.cpp
//  -OR-
//	g++ -lhdf5 -lz -lm -o Storagetry Storagetry.cpp ../FCComps/Storage.cpp ../FCComps/FCComp.cpp ../FCComps/MassStream.cpp ../FCComps/isoname.cpp ../FCComps/genlib.cpp

#include "../Storage.h"
#include "../MassStream.h"
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
	print("Test Initialize Empty Reprocessing Facility:");
	Storage sempty;
	print("\tName\t\t\t" + sempty.name);
	print("\tPass Number\t\t" + to_str(sempty.PassNum));
	print("\tFigure Directory\t" + sempty.figdir);
	print("");

	int itrack_arr [] = {882270, 862270, 922350, 922360, 922380, 942390, 952421, 822070};
	set<int> itrack (itrack_arr, itrack_arr+8); 	

	print("Test Initialize Reprocessing Facility:");
	Storage s (itrack, "Storage");
	print("\tName\t\t\t" + s.name);
	print("\tPass Number\t\t" + to_str(s.PassNum));
	print("\tFigure Directory\t" + s.figdir);
	print("");

	CompDict samplemass;
	samplemass[922350] = 10.0;
	samplemass[882270] = 0.0;
	samplemass[862270] = 1.0;
        samplemass[922360] = 1.0;
        samplemass[922380] = 50.0;
        samplemass[942390] = 20.0;
        samplemass[952421] = 5.0;

	print("Try Decaying Implicit Mass for 1000 days.");
	s.decay_time = 1000.0 * 24.0 * 3600.0;
	s.ms_feed = MassStream (samplemass);
	s.doCalc().print_ms();
	print("");

	print("Try Decaying CompDict Mass for 2000 days.");
	s.decay_time = 2000.0 * 24.0 * 3600.0;
	s.doCalc(samplemass).print_ms();
	print("");

	print("Try Decaying Mass Stream Mass for 3000 days.");
	s.decay_time = 3000.0 * 24.0 * 3600.0;
	MassStream ms (samplemass);
	s.doCalc(ms).print_ms();
	print("");

	double t;
	print("Try Decaying Implicit Mass for t = 4000 days.");
	t = 4000.0 * 24.0 * 3600.0;
	s.ms_feed = MassStream (samplemass);
	s.doCalc(t).print_ms();
	print("");
		
	print("Try Decaying CompDict Mass for t = 5000 days.");
	t = 5000.0 * 24.0 * 3600.0;
	s.doCalc(samplemass, t).print_ms();
	print("");
		
	print("Try Decaying CompDict Mass for t = 6000 days.");
	t = 6000.0 * 24.0 * 3600.0;
	s.doCalc(ms, t).print_ms();
	print("");

	print("Try Decaying Mass From a File for t = 7000 days.");
	t = 7000.0 * 24.0 * 3600.0;
	MassStream ms2 ("MassStreamtry02.txt");
	set<int> itrack2;
	for (CompIter ci = ms2.comp.begin(); ci != ms2.comp.end(); ci++ )
	{
		itrack2.insert(ci->first);
	}
	Storage s2 (itrack2, "Storage");
	s2.doCalc(ms2, t).print_ms();


	return 0;
}
