#include <cmath>
#include <cassert>
#include "Sscma.h"
#include "Subsets.h"

namespace leapsnbnds  {
long unsigned sbsetind,maxsbst;
}                                                  
extern sbset    **sbsarr;

using namespace leapsnbnds;

sbset *csbset(vind n,vind* v,real c)
{
	assert(sbsetind < maxsbst);
	sbset    *s = sbsarr[sbsetind++];
	s->nvar_ = n; 
	for (vind i=0;i<n;i++) s->actvar_[i] = v[i];
	s->crt_ = c;
	return s;
}

void dsbset(sbset *s)
{
	assert(sbsetind > 0);
	(sbsarr[s->pos]=sbsarr[--sbsetind])->pos = s->pos;
	(sbsarr[sbsetind]=s)->pos = sbsetind;
	return;
}

sbset::sbset(vind p,unsigned n)
:  pos(p), nvar_(n)
{
	actvar_ = new vind[nvar_];
}

sbset::~sbset(void)
{
	delete[] actvar_;
}
