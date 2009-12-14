#include <iostream>
#include <cstdlib>
#include "Sscma.h" 

using namespace std; 

namespace extendedleaps {

void msg(const string& s)
{
	cout << s;
}

void errmsg(const string& s)
{
	cerr << s;
	exit(1);	
}

}

