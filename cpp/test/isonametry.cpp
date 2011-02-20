// compile with g++ -o isonametry isonametry.cpp ../FCComps/isoname.cpp ../FCComps/genlib.cpp

#include "../genlib.h"
#include "../isoname.h"
#include <iostream>

using namespace std;
using namespace isoname;

void print( string s)
{
	cout << s << "\n";
}

void printGroupSize(ElementGroup eg, string name)
{
	cout << "\t" << name << " set size = " << eg.size() << "\n";
}

void TestNuclide(string l, int z, int m)
{
	print("Testing " + l + ":");
	print("\tCurrent Form of " + l + " is " + CurrentForm(l) );
	print("\tCurrent Form of " + to_str(z) + " is " + CurrentForm(z) );
	print("\tCurrent Form of " + to_str(m) + " is " + CurrentForm(m) );
	print("");
	print("\t" + l + " can change into " + to_str(LLAAAM_2_zzaaam(l)) + " or " + to_str(LLAAAM_2_MCNP(l)) );
	print("\t" + to_str(z) + " can change into " + zzaaam_2_LLAAAM(z) + " or " + to_str(zzaaam_2_MCNP(z)) );
	print("\t" + to_str(m) + " can change into " + MCNP_2_LLAAAM(m) + " or " + to_str(MCNP_2_zzaaam(m)) );
	print("");
	cout << "\tMixed " << l <<  " converts to " << mixed_2_LLAAAM(l) << ", " << mixed_2_zzaaam(l) << ", and "  << mixed_2_MCNP(l) << "\n";
	cout << "\tMixed " << z <<  " converts to " << mixed_2_LLAAAM(z) << ", " << mixed_2_zzaaam(z) << ", and "  << mixed_2_MCNP(z) << "\n";
	cout << "\tMixed " << m <<  " converts to " << mixed_2_LLAAAM(m) << ", " << mixed_2_zzaaam(m) << ", and "  << mixed_2_MCNP(m) << "\n";
}

int main()
{
	print( "Check LLzz[BE] = " + to_str(LLzz["BE"]) ); 
	print( "Check LLzz[B] = " + to_str(LLzz["B"]) ); 


	print( "Check LLzz[U] = " + to_str(LLzz["U"]) ); 
	print( "Check zzLL[92] = " + zzLL[92] ); 
	print("");

	print("Check Elemental Group Arrays:");
	print("\tLAN_Array = {" + LAN_Array[0] + ", " + LAN_Array[1] + ", ...}" );
	print("\tACT_Array = {" + ACT_Array[0] + ", " + ACT_Array[1] + ", ...}" );
	print("\tTRU_Array = {" + TRU_Array[0] + ", " + TRU_Array[1] + ", ...}" );
	print("\tMA_Array  = {" + MA_Array[0] + ", " +  MA_Array[1] + ", ...}" );
	print("\tFP_Array  = {" + FP_Array[0] + ", " +  FP_Array[1] + ", ...}" );
	print("");

	print("Check Elemental Group Sets:");
	printGroupSize(LAN, "LAN" );
	printGroupSize(ACT, "ACT" );
	printGroupSize(TRU, "TRU" );
	printGroupSize(MA,  "MA" );
	printGroupSize(FP,  "FP" );
	print("");

	print("Testing Exceptions:");
	try
	{
		throw NotANuclide;
	}
	catch (exception& e)
	{
		print("\tCaught NotANuclide");
	}
	try
	{
		throw IndeterminateNuclideForm;
	}
	catch (exception& e)
	{
		print("\tCaught IndeterminateNuclideForm");
	}
	print("");

	string LLZZZM;
	int zzaaam, MCNP;

	LLZZZM = "U235";
	zzaaam = 922350;
	MCNP = 92235;

	TestNuclide(LLZZZM, zzaaam, MCNP);
	print("");

	LLZZZM = "Am-242m";
	zzaaam = 952421;
	MCNP = 95242;

	TestNuclide(LLZZZM, zzaaam, MCNP);
	return 0;
}
