#pragma once
#include "FIELD.h"

class ADVECTION_METHOD
{
public:

	template<class TT>
	void SemiLagrangian(FIELD<TV3>* vel, const TS dt, FIELD<TT>* in_field, FIELD<TT>* out_field);
};