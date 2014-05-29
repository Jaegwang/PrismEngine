#pragma once

#include <omp.h>
#include "FIELD.h"

template<class TT>
void RasterizeParticleToField(FIELD<TT>& val_field, FIELD<FLT>& weight_field, const ARRAY<Vec3>& pos_array, const ARRAY<TT>& val_array)
{
	int size = pos_array.Size();
	int i,j,k;

	int i_start, i_end;
	int j_start, j_end;
	int k_start, k_end;

	GRID grid = val_field.Grid();

	#pragma omp parallel for
	for(int p=0; p<grid.ijk_res_; p++)
	{
		val_field.Set(p, TT());
		weight_field.Set(p, (FLT)0);
	}

//	#pragma omp parallel for
	for(int p=0; p<size; p++)
	{
		Vec3 pos = pos_array.Get(p);
		TT   val = val_array.Get(p);

		grid.CellCenterIndex(pos,i,j,k);		

		i_start = MAX(0, i-1);
		j_start = MAX(0, j-1);
		k_start = MAX(0, k-1);

		i_end = MIN(grid.i_res_-1, i+1);
		j_end = MIN(grid.j_res_-1, j+1);
		k_end = MIN(grid.k_res_-1, k+1);

		for(int i=i_start; i<=i_end; i++) for(j=j_start; j<=j_end; j++) for(k=k_start; k<=k_end; k++)
		{
			int ix = grid.Index3Dto1D(i,j,k);

			Vec3 cell_center = grid.CellCenterPosition(i,j,k);
			Vec3 deviation = cell_center-pos;

			FLT w = QuadBSplineKernel(deviation.x*grid.one_over_dx_)*
				    QuadBSplineKernel(deviation.y*grid.one_over_dx_)*
					QuadBSplineKernel(deviation.z*grid.one_over_dx_);

//			#pragma omp critical
			{
				weight_field.Set(ix, weight_field.Get(ix)+w);
				val_field.Set(ix, val_field.Get(ix)+val*w);
			}
		}
	}	

	#pragma omp parallel for
	for(int p=0; p<grid.ijk_res_; p++)
	{
		FLT weight = weight_field.Get(p);
		
		if(weight >= (FLT)1e-06)
			val_field.Set(p, val_field.Get(p)/weight);
		else
			val_field.Set(p, TT());
	}
}


template<class TT>
void RasterizeFieldToParticle(ARRAY<Vec3>& pos_array, ARRAY<TT>& val_array, const FIELD<TT>& val_field)
{
	int size = pos_array.Size();
	int i,j,k;

	int i_start, i_end;
	int j_start, j_end;
	int k_start, k_end;

	GRID grid = val_field.Grid();

	#pragma omp parallel for
	for(int p=0; p<size; p++)
	{
		Vec3 pos = pos_array.Get(p);
		grid.CellCenterIndex(pos,i,j,k);		

		i_start = MAX(0, i-1);
		j_start = MAX(0, j-1);
		k_start = MAX(0, k-1);

		i_end = MIN(grid.i_res_-1, i+1);
		j_end = MIN(grid.j_res_-1, j+1);
		k_end = MIN(grid.k_res_-1, k+1);

		TT val_p = TT();
		FLT w_p = (FLT)0;

		for(int i=i_start; i<=i_end; i++) for(j=j_start; j<=j_end; j++) for(k=k_start; k<=k_end; k++)
		{
			int ix = grid.Index3Dto1D(i,j,k);

			TT val_g = val_field.Get(ix);

			Vec3 cell_center = grid.CellCenterPosition(i,j,k);
			Vec3 deviation = cell_center-pos;

			FLT w = QuadBSplineKernel(deviation.x*grid.one_over_dx_)*
				    QuadBSplineKernel(deviation.y*grid.one_over_dx_)*
					QuadBSplineKernel(deviation.z*grid.one_over_dx_);

			val_p += val_g * w;
			w_p += w;
		}

		if(w_p > (FLT)1e-06) 
			val_array.Set(p, val_p/w_p);
		else 
			val_array.Set(p, TT());
	}	
}

template void RasterizeParticleToField(FIELD<Vec3>& val_field, FIELD<FLT>& weight_field, const ARRAY<Vec3>& pos_array, const ARRAY<Vec3>& val_array);
template void RasterizeFieldToParticle(ARRAY<Vec3>& pos_array, ARRAY<Vec3>& val_array, const FIELD<Vec3>& val_field);