#include "../genlib.h"
#include <iostream>

using namespace std;

void Print( string s)
{
	cout << s << "\n";
}

int main()
{
	string n = "1";
	int nn = to_int(n);
	string nnn = to_str(nn);
	Print(nnn);

	string t = "Np-237";
	Print(t);
	Print( ToUpper(t) ); 
	Print( ToLower(t) ); 
	Print( Strip(t, "-") ); 
	Print( MultiStrip(t, "-"+digits) ); 
	Print( LastChar(t) ); 
	Print( SubFromEnd(t, -3,3) ); 

	Print( to_str(ChainGreaterCompare(0, 10, 40)) );
	Print( to_str(ChainGreaterCompare(100, 10, 40)) );

	Print( to_str(SubInString("High-Five", "Five")) );
	Print( to_str(SubInString("High-Five", "no")) );

	return 0;
}
