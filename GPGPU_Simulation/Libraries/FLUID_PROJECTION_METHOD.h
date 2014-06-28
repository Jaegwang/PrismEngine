#pragma once 

#include "GRID.h"

class FLUID_PROJECT_METHOD
{
public:

	FLUID_PROJECT_METHOD();
	~FLUID_PROJECT_METHOD();

public:

	void Project(GRID& grid, const TV3* velocity_arr, TS* divergence_arr, TS** pressure_arr, TS** pressure_arr_temp)
	{
		#pragma omp parallel for
		for(int p = 0; p < grid.ijk_res_; p++) 
		{
			(*pressure_arr)[p] = (TS)0;
			(*pressure_arr_temp)[p] = (TS)0;
		}	
	}

	void ComputeDivergence(GRID& grid, const TV3* velocity_arr, TS* divergence_arr)
	{
		#pragma omp parallel for
		for(int p = 0; p < grid.ijk_res_; p++) 
		{
			int i, j, k;
			grid.Index1Dto3D(p, i, j, k);

			TS div = (TS)0;
				
			if(grid.IsGhostCell(i, j, k) == false)
			{
				div += (velocity_arr[p+1].x - velocity_arr[p-1].x)*grid.one_over_dx_*(TS)0.5;
				div += (velocity_arr[p+grid.i_res_].y - velocity_arr[p-grid.i_res_].y)*grid.one_over_dy_*(TS)0.5;
				div += (velocity_arr[p+grid.ij_res_].z - velocity_arr[p+grid.ij_res_].z)*grid.one_over_dz_*(TS)0.5;

				divergence_arr[p] = div;
			}
			else
			{
				divergence_arr[p] = (TS)0;
			}
		}	
	}

	void ComputePressureJacobi(GRID& grid, TS* divergence_arr, TS** pressure_arr, TS** pressure_arr_ghost)
	{
		TS dxdydz = POW2(grid.dx_*grid.dy_*grid.dz_);
		TS dxdy   = POW2(grid.dx_*grid.dy_);
		TS dydz   = POW2(grid.dy_*grid.dz_);
		TS dxdz   = POW2(grid.dx_*grid.dz_);

		TS one_over_h = (dxdy + dydz + dxdz)*(TS)0.5;

		#pragma omp parallel for
		for(int p = 0; p < grid.ijk_res_; p++) 
		{
			int i, j, k;
			grid.Index1Dto3D(p, i, j, k);
				
			if(grid.IsGhostCell(i, j, k) == false)
			{
				TS p_i1 = *pressure_arr[p+1];
				TS p_i0 = *pressure_arr[p-1];

				TS p_j1 = *pressure_arr[p+grid.i_res_];
				TS p_j0 = *pressure_arr[p-grid.i_res_];

				TS p_k1 = *pressure_arr[p+grid.ij_res_];
				TS p_k0 = *pressure_arr[p-grid.ij_res_];

				(*pressure_arr_ghost)[p] = (dydz*(p_i1+p_i0) + dxdz*(p_j1+p_j0) + dxdy*(p_k1+p_k0) - divergence_arr[p]*dxdydz) * one_over_h;
			}
			else
			{
				(*pressure_arr_ghost)[p] = (TS)0;
			}
		}

		TS* temp = *pressure_arr_ghost;
		*pressure_arr_ghost = *pressure_arr;
		*pressure_arr = temp;
	}
};