
#pragma once

#include "MATH_CORE.h"

static TS QuadBSplineKernel(const TS d, const TS one_over_h)
{
	const TS q = d*one_over_h;

	if(q < -(TS)3/(TS)2) return 0;
	else if(q <= -(TS)1/(TS)2) return (TS)0.5*q*q + (TS)3/(TS)2*q + (TS)9/(TS)8; 
	else if(q <=  (TS)1/(TS)2) return -q*q + (TS)3/(TS)4;
	else if(q <=  (TS)3/(TS)2) return (TS)0.5*q*q - (TS)3/(TS)2*q + (TS)9/(TS)8; 
	else return 0;
}

static TS QuadBSplineKernel(const TS q)
{
	if(q < -(TS)3/(TS)2) return 0;
	else if(q <= -(TS)1/(TS)2) return (TS)0.5*q*q + (TS)3/(TS)2*q + (TS)9/(TS)8; 
	else if(q <=  (TS)1/(TS)2) return -q*q + (TS)3/(TS)4;
	else if(q <=  (TS)3/(TS)2) return (TS)0.5*q*q - (TS)3/(TS)2*q + (TS)9/(TS)8; 
	else return 0;
}

static TS QuadBSplineKernelGradient(const TS d, const TS one_over_h)
{
	const TS q = d*one_over_h;

	if(q < -(TS)3/(TS)2) return 0;
	else if(q <= -(TS)1/(TS)2) return (TS)q + (TS)3/(TS)2;
	else if(q <=  (TS)1/(TS)2) return -(TS)2*q;
	else if(q <=  (TS)3/(TS)2) return (TS)q - (TS)3/(TS)2; 
	else return 0;
}

static TS QuadBSplineKernel(const TV3 d, const TS one_over_dx, const TS one_over_dy, const TS one_over_dz)
{
	return QuadBSplineKernel(d.x, one_over_dx)*QuadBSplineKernel(d.y, one_over_dy)*QuadBSplineKernel(d.z, one_over_dz);
}

static TV3 QuadBSplineKernelGradient(const TV3 d, const TS one_over_dx, const TS one_over_dy, const TS one_over_dz)
{ 
	return TV3(QuadBSplineKernelGradient(d.x, one_over_dx), QuadBSplineKernelGradient(d.y, one_over_dy), QuadBSplineKernelGradient(d.z, one_over_dz));
}

static void QuadBSplineKernel(const TV3& d, const TV3& one_over_h, TV3& w)
{
	w.x = QuadBSplineKernel(d.x, one_over_h.x);
	w.y = QuadBSplineKernel(d.y, one_over_h.y);
	w.z = QuadBSplineKernel(d.z, one_over_h.z);
}

static void QuadBSplineKernelGradient(const TV3& d, const TV3& one_over_h, TV3& w)
{ 
	w.x = QuadBSplineKernelGradient(d.x, one_over_h.x);
	w.y = QuadBSplineKernelGradient(d.y, one_over_h.y);
	w.z = QuadBSplineKernelGradient(d.z, one_over_h.z);
}

/*
static TS MPMSplineKernel(const TS d, const TS one_over_h) 
{
	const TS x = d * one_over_h;
	const TS abs_x = ABS(x);

	if      ((TS)0 <= abs_x && abs_x < (TS)1) return (TS)0.5*POW3(abs_x) - POW2(x) + (TS)0.6666;
	else if ((TS)1 <= abs_x && abs_x < (TS)2) return (TS)-0.1666*POW3(abs_x) + POW2(x) - (TS)2*abs_x + (TS)1.3333;
	else                                    return (TS)0;	
}

static TS MPMSplineKernelGradient(const TS d, const TS one_over_h)
{
	const TS x = d * one_over_h;
	const TS abs_x = ABS(x);

	if      ((TS)0 <= abs_x && abs_x < (TS)1) return (TS)1.5 * POW2(x) - (TS)2 * x;
	else if ((TS)1 <= abs_x && abs_x < (TS)2) return (TS)-0.5 * POW2(x) + (TS)2 * x - (TS)2;
	else                                    return (TS)0;
}
*/

static TS MPMSplineKernel(const TS d, const TS one_over_h)
{
	const TS q = d*one_over_h;

	//TODO: 1. check coefficients. 2. optimize.
	if (q < (TS)-2) return (TS)0;
	else if (q < (TS)-1) return (TS)1 / 6 * q*q*q + q*q + (TS)2 * q + (TS)4 / (TS)3;
	else if (q < (TS)0)  return -(TS)0.5*q*q*q - q*q + (TS)2 / (TS)3;
	else if (q < (TS)1)  return (TS)0.5*q*q*q - q*q + (TS)2 / (TS)3;
	else if (q < (TS)2)  return -(TS)1 / 6 * q*q*q + q*q - (TS)2 * q + (TS)4 / (TS)3;
	else return (TS)0;
}

static TS MPMSplineKernelGradient(const TS d, const TS one_over_h)
{
	const TS q = d*one_over_h;

	//TODO: 1. check coefficients. 2. optimize.
	if (q < (TS)-2) return (TS)0;
	else if (q < (TS)-1) return  ((TS)0.5*q*q + (TS)2 * q + (TS)2);
	else if (q < (TS)0)  return (-(TS)1.5*q*q - (TS)2 * q);
	else if (q < (TS)1)  return  ((TS)1.5*q*q - (TS)2 * q);
	else if (q < (TS)2)  return (-(TS)0.5*q*q + (TS)2 * q - (TS)2);
	else return (TS)0;
}

static TS MPMSplineKernel(const TV3 d, const TS one_over_h)
{
	return MPMSplineKernel(d.x, one_over_h)*MPMSplineKernel(d.y, one_over_h)*MPMSplineKernel(d.z, one_over_h);
}

static TS MPMSplineKernel(const TV3 d, const TS one_over_dx, const TS one_over_dy, const TS one_over_dz)
{
	return MPMSplineKernel(d.x, one_over_dx)*MPMSplineKernel(d.y, one_over_dy)*MPMSplineKernel(d.z, one_over_dz);
}


static TV3 MPMSplineKernelGradient(const TV3 d, const TS one_over_h)
{
	return TV3(MPMSplineKernelGradient(d.x, one_over_h), MPMSplineKernelGradient(d.y, one_over_h), MPMSplineKernelGradient(d.z, one_over_h));
}

static TV3 MPMSplineKernelGradient(const TV3 d, const TS one_over_dx, const TS one_over_dy, const TS one_over_dz)
{
	return TV3(MPMSplineKernelGradient(d.x, one_over_dx), MPMSplineKernelGradient(d.y, one_over_dy), MPMSplineKernelGradient(d.z, one_over_dz));
}