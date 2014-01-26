
#pragma once

#include "VECTOR3_T.h"


/*
static T MPMSplineKernel(const T d, const T one_over_h) restrict(cpu,amp)
{
	const T x = d * one_over_h;
	const T abs_x = ABS(x);

	if      ((T)0 <= abs_x && abs_x < (T)1) return (T)0.5*POW3(abs_x) - POW2(x) + (T)0.6666;
	else if ((T)1 <= abs_x && abs_x < (T)2) return (T)-0.1666*POW3(abs_x) + POW2(x) - (T)2*abs_x + (T)1.3333;
	else                                    return (T)0;	
}

static T MPMSplineKernelGradient(const T d, const T one_over_h) restrict(cpu, amp)
{
	const T x = d * one_over_h;
	const T abs_x = ABS(x);

	if      ((T)0 <= abs_x && abs_x < (T)1) return (T)1.5 * POW2(x) - (T)2 * x;
	else if ((T)1 <= abs_x && abs_x < (T)2) return (T)-0.5 * POW2(x) + (T)2 * x - (T)2;
	else                                    return (T)0;
}
*/

static T MPMSplineKernel(const T d, const T one_over_h) restrict(cpu, amp)
{
	const T q = d*one_over_h;

	//TODO: 1. check coefficients. 2. optimize.
	if (q < (T)-2) return (T)0;
	else if (q < (T)-1) return (T)1 / 6 * q*q*q + q*q + (T)2 * q + (T)4 / (T)3;
	else if (q < (T)0)  return -(T)0.5*q*q*q - q*q + (T)2 / (T)3;
	else if (q < (T)1)  return (T)0.5*q*q*q - q*q + (T)2 / (T)3;
	else if (q < (T)2)  return -(T)1 / 6 * q*q*q + q*q - (T)2 * q + (T)4 / (T)3;
	else return (T)0;
}

static T MPMSplineKernelGradient(const T d, const T one_over_h) restrict(cpu, amp)
{
	const T q = d*one_over_h;

	//TODO: 1. check coefficients. 2. optimize.
	if (q < (T)-2) return (T)0;
	else if (q < (T)-1) return  ((T)0.5*q*q + (T)2 * q + (T)2);
	else if (q < (T)0)  return (-(T)1.5*q*q - (T)2 * q);
	else if (q < (T)1)  return  ((T)1.5*q*q - (T)2 * q);
	else if (q < (T)2)  return (-(T)0.5*q*q + (T)2 * q - (T)2);
	else return (T)0;
}

static T MPMSplineKernel(const Vec3T d, const T one_over_h) restrict(cpu, amp)
{
	return MPMSplineKernel(d.x, one_over_h)*MPMSplineKernel(d.y, one_over_h)*MPMSplineKernel(d.z, one_over_h);
}

static T MPMSplineKernel(const Vec3T d, const T one_over_dx, const T one_over_dy, const T one_over_dz) restrict(cpu, amp)
{
	return MPMSplineKernel(d.x, one_over_dx)*MPMSplineKernel(d.y, one_over_dy)*MPMSplineKernel(d.z, one_over_dz);
}


static Vec3T MPMSplineKernelGradient(const Vec3T d, const T one_over_h) restrict(cpu, amp)
{
	return Vec3T(MPMSplineKernelGradient(d.x, one_over_h), MPMSplineKernelGradient(d.y, one_over_h), MPMSplineKernelGradient(d.z, one_over_h));
}

static Vec3T MPMSplineKernelGradient(const Vec3T d, const T one_over_dx, const T one_over_dy, const T one_over_dz) restrict(cpu, amp)
{
	return Vec3T(MPMSplineKernelGradient(d.x, one_over_dx), MPMSplineKernelGradient(d.y, one_over_dy), MPMSplineKernelGradient(d.z, one_over_dz));
}