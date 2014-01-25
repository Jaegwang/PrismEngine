
#pragma once

#include "VECTOR3_T.h"


static T MPMSplineKernel(const T d, const T one_over_h) restrict(cpu,amp)
{
	const T x = d * one_over_h;
	const T abs_x = ABS(x);

	if      ((T)0 <= abs_x && abs_x < (T)1) return (T)0.5*POW3(abs_x) - POW2(x) + (T)0.6666;
	else if ((T)1 <= abs_x && abs_x < (T)2) return (T)-0.1666*POW3(abs_x) + POW2(x) - (T)2*abs_x + (T)1.3333;
	else                                    return (T)0;	
}

static T MPMSplineKernel(const Vec3T d, const T one_over_h) restrict(cpu, amp)
{
	return MPMSplineKernel(d.x, one_over_h)*MPMSplineKernel(d.y, one_over_h)*MPMSplineKernel(d.z, one_over_h);
}

static T MPMSplineKernel(const Vec3T d, const T one_over_dx, const T one_over_dy, const T one_over_dz) restrict(cpu, amp)
{
	return MPMSplineKernel(d.x, one_over_dx)*MPMSplineKernel(d.y, one_over_dy)*MPMSplineKernel(d.z, one_over_dz);
}