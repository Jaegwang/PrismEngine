#pragma once
#include "FIELD.h"

enum BOUNDARY_CELL
{
	BND_NULL,
	BND_WALL,
	BND_FULL,
};

class PROJECTION_METHOD
{
public:
	
	void ComputeDivergence(FIELD<Vec3>* vel, FIELD<int>* bnd, FIELD<FLT>* div);
	void ComputePressure(FIELD<FLT>* div, FIELD<int>* bnd, FIELD<FLT>* press, FIELD<FLT>* press_temp, const int itr);
	void ComputeVelocity(FIELD<FLT>* press, FIELD<int>* bnd, FIELD<Vec3>* vel);

	void Jacobi(FIELD<Vec3>* vel, FIELD<int>* bnd, FIELD<FLT>* div, FIELD<FLT>* press, FIELD<FLT>* press_temp, const int itr)
	{
		ComputeDivergence(vel, bnd, div);

		ComputePressure(div, bnd, press, press_temp, itr);

		ComputeVelocity(press, bnd, vel);
	}
};
