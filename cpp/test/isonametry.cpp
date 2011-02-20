// compile with g++ -o isonametry isonametry.cpp ../FCComps/isoname.cpp ../FCComps/genlib.cpp

#include "../genlib.h"
#include "../isoname.h"
#include <iostream>

using namespace std;
using namespace isoname;

void Print( string s)
{
	cout << s << "\n";
}

void PrintGroupSize(ElementGroup eg, string name)
{
	cout << "\t" << name << " set size = " << eg.size() << "\n";
}

void TestNuclide(string l, int z, int m)
{
	Print("Testing " + l + ":");
	Print("\tCurrent Form of " + l + " is " + CurrentForm(l) );
	Print("\tCurrent Form of " + to_str(z) + " is " + CurrentForm(z) );
	Print("\tCurrent Form of " + to_str(m) + " is " + CurrentForm(m) );
	Print("");
	Print("\t" + l + " can change into " + to_str(LLAAAM_2_zzaaam(l)) + " or " + to_str(LLAAAM_2_MCNP(l)) );
	Print("\t" + to_str(z) + " can change into " + zzaaam_2_LLAAAM(z) + " or " + to_str(zzaaam_2_MCNP(z)) );
	Print("\t" + to_str(m) + " can change into " + MCNP_2_LLAAAM(m) + " or " + to_str(MCNP_2_zzaaam(m)) );
	Print("");
	cout << "\tMixed " << l <<  " converts to " << mixed_2_LLAAAM(l) << ", " << mixed_2_zzaaam(l) << ", and "  << mixed_2_MCNP(l) << "\n";
	cout << "\tMixed " << z <<  " converts to " << mixed_2_LLAAAM(z) << ", " << mixed_2_zzaaam(z) << ", and "  << mixed_2_MCNP(z) << "\n";
	cout << "\tMixed " << m <<  " converts to " << mixed_2_LLAAAM(m) << ", " << mixed_2_zzaaam(m) << ", and "  << mixed_2_MCNP(m) << "\n";
}

int main()
{
	Print( "Check LLzz[BE] = " + to_str(LLzz["BE"]) ); 
	Print( "Check LLzz[B] = " + to_str(LLzz["B"]) ); 


	Print( "Check LLzz[U] = " + to_str(LLzz["U"]) ); 
	Print( "Check zzLL[92] = " + zzLL[92] ); 
	Print("");

	Print("Check Elemental Group Arrays:");
	Print("\tLAN_Array = {" + LAN_Array[0] + ", " + LAN_Array[1] + ", ...}" );
	Print("\tACT_Array = {" + ACT_Array[0] + ", " + ACT_Array[1] + ", ...}" );
	Print("\tTRU_Array = {" + TRU_Array[0] + ", " + TRU_Array[1] + ", ...}" );
	Print("\tMA_Array  = {" + MA_Array[0] + ", " +  MA_Array[1] + ", ...}" );
	Print("\tFP_Array  = {" + FP_Array[0] + ", " +  FP_Array[1] + ", ...}" );
	Print("");

	Print("Check Elemental Group Sets:");
	PrintGroupSize(LAN, "LAN" );
	PrintGroupSize(ACT, "ACT" );
	PrintGroupSize(TRU, "TRU" );
	PrintGroupSize(MA,  "MA" );
	PrintGroupSize(FP,  "FP" );
	Print("");

	Print("Testing Exceptions:");
	try
	{
		throw NotANuclide;
	}
	catch (exception& e)
	{
		Print("\tCaught NotANuclide");
	}
	try
	{
		throw IndeterminateNuclideForm;
	}
	catch (exception& e)
	{
		Print("\tCaught IndeterminateNuclideForm");
	}
	Print("");

	string LLZZZM;
	int zzaaam, MCNP;

	LLZZZM = "U235";
	zzaaam = 922350;
	MCNP = 92235;

	TestNuclide(LLZZZM, zzaaam, MCNP);
	Print("");

	LLZZZM = "Am-242m";
	zzaaam = 952421;
	MCNP = 95242;

	TestNuclide(LLZZZM, zzaaam, MCNP);
	return 0;
}
