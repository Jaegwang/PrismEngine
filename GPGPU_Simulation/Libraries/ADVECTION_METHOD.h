#pragma once
#include "FIELD.h"

class ADVECTION_METHOD
{
public:

	template<class TT>
	void SemiLagrangian(FIELD<Vec3>* vel, const FLT dt, FIELD<TT>* in_field, FIELD<TT>* out_field);
};