
#include "ADVECTION_METHOD.h"
#include "GENERAL_MACRO.h"

template<class TT>
void ADVECTION_METHOD::SemiLagrangian(FIELD<TV3>* vel, const TS dt, FIELD<TT>* in_field, FIELD<TT>* out_field)
{
	GRID grid = out_field->Grid();

	FOR_EACH_PARALLER(k, 0, grid.k_res_-1)
	{
		for(int j=0; j<grid.j_res_; j++) for(int i=0; i<grid.i_res_; i++)	
		{
			TV3 cell_center = grid.CellCenterPosition(i,j,k);
			TV3 velocity = vel->Get(cell_center);

			out_field->Set(i,j,k, in_field->Get(cell_center-velocity*dt));
		}
	}
}

template void ADVECTION_METHOD::SemiLagrangian(FIELD<TV3>* vel, const TS dt, FIELD<float>* in_field, FIELD<float>* out_field);
template void ADVECTION_METHOD::SemiLagrangian(FIELD<TV3>* vel, const TS dt, FIELD<double>* in_field, FIELD<double>* out_field);
template void ADVECTION_METHOD::SemiLagrangian(FIELD<TV3>* vel, const TS dt, FIELD<TV3>* in_field, FIELD<TV3>* out_field);
