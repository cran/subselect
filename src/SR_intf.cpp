#include <iostream>
#include "R.h"
#include "Sscma.h" 

namespace extendedleaps {

void msg(const std::string& s)
{
	Rprintf("%s",s.c_str());
}

void errmsg(const std::string& s)
{
	error("\nError: %s",s.c_str());
}

}

