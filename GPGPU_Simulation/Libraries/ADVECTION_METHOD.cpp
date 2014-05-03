
#include "ADVECTION_METHOD.h"

template<class TT>
void ADVECTION_METHOD::SemiLagrangian(FIELD<Vec3>* vel, const FLT dt, FIELD<TT>* in_field, FIELD<TT>* out_field)
{
	GRID grid = out_field->Grid();

	#pragma omp parallel for
	for(int k=0; k<grid.k_res_; k++) for(int j=0; j<grid.j_res_; j++) for(int i=0; i<grid.i_res_; i++)
	{
		Vec3 cell_center = grid.CellCenterPosition(i,j,k);
		Vec3 velocity = vel->Get(cell_center);

		out_field->Set(i,j,k, in_field->Get(cell_center-velocity*dt));
	}
}

template void ADVECTION_METHOD::SemiLagrangian(FIELD<Vec3>* vel, const FLT dt, FIELD<float>* in_field, FIELD<float>* out_field);
template void ADVECTION_METHOD::SemiLagrangian(FIELD<Vec3>* vel, const FLT dt, FIELD<double>* in_field, FIELD<double>* out_field);
template void ADVECTION_METHOD::SemiLagrangian(FIELD<Vec3>* vel, const FLT dt, FIELD<Vec3>* in_field, FIELD<Vec3>* out_field);