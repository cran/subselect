#include "Sscma.h" 
#include "R.h"


void leapsnbnds::msg(const std::string& s)
{
	Rprintf("%s",s.c_str());
}

void leapsnbnds::errmsg(const std::string& s)
{
	error("\nError: %s",s.c_str());
}
