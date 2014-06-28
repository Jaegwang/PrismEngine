
#pragma once

#include "MATH_DEFINITION.h"

static TV3 RandomVector()
{
	return TV3(((TS)rand()/(TS)RAND_MAX-(TS)0.5)*(TS)2, ((TS)rand()/(TS)RAND_MAX-(TS)0.5)*(TS)2, ((TS)rand()/(TS)RAND_MAX-(TS)0.5)*(TS)2);
}