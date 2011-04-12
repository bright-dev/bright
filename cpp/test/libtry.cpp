#include "../h5wrap.h"
#include <iostream>

using namespace std;

void print( string s)
{
	cout << s << "\n";
}

int main()
{
	string n = "1";
	int nn = to_int(n);
	string nnn = to_str(nn);
	print(nnn);

	string t = "Np-237";
	print(t);
	print( ToUpper(t) ); 
	print( ToLower(t) ); 
	print( Strip(t, "-") ); 
	print( MultiStrip(t, "-"+digits) ); 
	print( LastChar(t) ); 
	print( SubFromEnd(t, -3,3) ); 

	print( to_str(ChainGreaterCompare(0, 10, 40)) );
	print( to_str(ChainGreaterCompare(100, 10, 40)) );

	print( to_str(SubInString("High-Five", "Five")) );
	print( to_str(SubInString("High-Five", "no")) );

	return 0;
}
