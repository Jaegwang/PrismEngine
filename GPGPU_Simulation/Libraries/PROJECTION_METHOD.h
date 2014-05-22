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
	
	void DetermineDivergence(const FIELD<int>* bnd, const FIELD<Vec3>* vel, FIELD<FLT>* div);

	void DeterminePressure(const FIELD<int>* bnd, const FIELD<FLT>* div, FIELD<FLT>* press, FIELD<FLT>* press_temp, const int itr=20);

	void DetermineVelocity(const FIELD<int>* bnd, const FIELD<FLT>* press, FIELD<Vec3>* vel);


	void Jacobi(const FIELD<int>* bnd, FIELD<Vec3>* vel, FIELD<FLT>* div, FIELD<FLT>* press, FIELD<FLT>* press_temp, const int itr)
	{
		DetermineDivergence(bnd, vel, div);

		DeterminePressure(bnd, div, press, press_temp, itr);

		DetermineVelocity(bnd, press, vel);
	}
};