
#pragma once

#include "MATH_CORE.h"

static FLT QuadBSplineKernel(const FLT d, const FLT one_over_h)
{
	const FLT q = d*one_over_h;

	if(q < -(FLT)3/(FLT)2) return 0;
	else if(q <= -(FLT)1/(FLT)2) return (FLT)0.5*q*q + (FLT)3/(FLT)2*q + (FLT)9/(FLT)8; 
	else if(q <=  (FLT)1/(FLT)2) return -q*q + (FLT)3/(FLT)4;
	else if(q <=  (FLT)3/(FLT)2) return (FLT)0.5*q*q - (FLT)3/(FLT)2*q + (FLT)9/(FLT)8; 
	else return 0;
}

static FLT QuadBSplineKernelGradient(const FLT d, const FLT one_over_h)
{
	const FLT q = d*one_over_h;

	if(q < -(FLT)3/(FLT)2) return 0;
	else if(q <= -(FLT)1/(FLT)2) return (FLT)q + (FLT)3/(FLT)2;
	else if(q <=  (FLT)1/(FLT)2) return -(FLT)2*q;
	else if(q <=  (FLT)3/(FLT)2) return (FLT)q - (FLT)3/(FLT)2; 
	else return 0;
}

static FLT QuadBSplineKernel(const Vec3 d, const FLT one_over_dx, const FLT one_over_dy, const FLT one_over_dz)
{
	return QuadBSplineKernel(d.x, one_over_dx)*QuadBSplineKernel(d.y, one_over_dy)*QuadBSplineKernel(d.z, one_over_dz);
}

static Vec3 QuadBSplineKernelGradient(const Vec3 d, const FLT one_over_dx, const FLT one_over_dy, const FLT one_over_dz)
{ 
	return Vec3(QuadBSplineKernelGradient(d.x, one_over_dx), QuadBSplineKernelGradient(d.y, one_over_dy), QuadBSplineKernelGradient(d.z, one_over_dz));
}


/*
static FLT MPMSplineKernel(const FLT d, const FLT one_over_h) 
{
	const FLT x = d * one_over_h;
	const FLT abs_x = ABS(x);

	if      ((FLT)0 <= abs_x && abs_x < (FLT)1) return (FLT)0.5*POW3(abs_x) - POW2(x) + (FLT)0.6666;
	else if ((FLT)1 <= abs_x && abs_x < (FLT)2) return (FLT)-0.1666*POW3(abs_x) + POW2(x) - (FLT)2*abs_x + (FLT)1.3333;
	else                                    return (FLT)0;	
}

static FLT MPMSplineKernelGradient(const FLT d, const FLT one_over_h)
{
	const FLT x = d * one_over_h;
	const FLT abs_x = ABS(x);

	if      ((FLT)0 <= abs_x && abs_x < (FLT)1) return (FLT)1.5 * POW2(x) - (FLT)2 * x;
	else if ((FLT)1 <= abs_x && abs_x < (FLT)2) return (FLT)-0.5 * POW2(x) + (FLT)2 * x - (FLT)2;
	else                                    return (FLT)0;
}
*/

static FLT MPMSplineKernel(const FLT d, const FLT one_over_h)
{
	const FLT q = d*one_over_h;

	//TODO: 1. check coefficients. 2. optimize.
	if (q < (FLT)-2) return (FLT)0;
	else if (q < (FLT)-1) return (FLT)1 / 6 * q*q*q + q*q + (FLT)2 * q + (FLT)4 / (FLT)3;
	else if (q < (FLT)0)  return -(FLT)0.5*q*q*q - q*q + (FLT)2 / (FLT)3;
	else if (q < (FLT)1)  return (FLT)0.5*q*q*q - q*q + (FLT)2 / (FLT)3;
	else if (q < (FLT)2)  return -(FLT)1 / 6 * q*q*q + q*q - (FLT)2 * q + (FLT)4 / (FLT)3;
	else return (FLT)0;
}

static FLT MPMSplineKernelGradient(const FLT d, const FLT one_over_h)
{
	const FLT q = d*one_over_h;

	//TODO: 1. check coefficients. 2. optimize.
	if (q < (FLT)-2) return (FLT)0;
	else if (q < (FLT)-1) return  ((FLT)0.5*q*q + (FLT)2 * q + (FLT)2);
	else if (q < (FLT)0)  return (-(FLT)1.5*q*q - (FLT)2 * q);
	else if (q < (FLT)1)  return  ((FLT)1.5*q*q - (FLT)2 * q);
	else if (q < (FLT)2)  return (-(FLT)0.5*q*q + (FLT)2 * q - (FLT)2);
	else return (FLT)0;
}

static FLT MPMSplineKernel(const Vec3 d, const FLT one_over_h)
{
	return MPMSplineKernel(d.x, one_over_h)*MPMSplineKernel(d.y, one_over_h)*MPMSplineKernel(d.z, one_over_h);
}

static FLT MPMSplineKernel(const Vec3 d, const FLT one_over_dx, const FLT one_over_dy, const FLT one_over_dz)
{
	return MPMSplineKernel(d.x, one_over_dx)*MPMSplineKernel(d.y, one_over_dy)*MPMSplineKernel(d.z, one_over_dz);
}


static Vec3 MPMSplineKernelGradient(const Vec3 d, const FLT one_over_h)
{
	return Vec3(MPMSplineKernelGradient(d.x, one_over_h), MPMSplineKernelGradient(d.y, one_over_h), MPMSplineKernelGradient(d.z, one_over_h));
}

static Vec3 MPMSplineKernelGradient(const Vec3 d, const FLT one_over_dx, const FLT one_over_dy, const FLT one_over_dz)
{
	return Vec3(MPMSplineKernelGradient(d.x, one_over_dx), MPMSplineKernelGradient(d.y, one_over_dy), MPMSplineKernelGradient(d.z, one_over_dz));
}