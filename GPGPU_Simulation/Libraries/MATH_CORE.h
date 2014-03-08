
#pragma once

#include "MATH_DEFINITION.h"

static Vec3 RandomVector()
{
	return Vec3(((FLT)rand()/(FLT)RAND_MAX-(FLT)0.5)*(FLT)2, ((FLT)rand()/(FLT)RAND_MAX-(FLT)0.5)*(FLT)2, ((FLT)rand()/(FLT)RAND_MAX-(FLT)0.5)*(FLT)2);
}