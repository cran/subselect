#include "R.h"

#ifndef interfacef
void prmtend(const char *);
void errmsg(char *);
#define interfacef
#endif


void errmsg(char *s)
{
	error("\nError: %s",s);
}

void prmtend(const char *s)
{
	error("Error, when trying to allocate memory for %s",s);
}





