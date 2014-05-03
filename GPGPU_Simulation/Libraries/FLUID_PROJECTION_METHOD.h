#pragma once 

#include "GRID.h"

class FLUID_PROJECT_METHOD
{
public:

	FLUID_PROJECT_METHOD();
	~FLUID_PROJECT_METHOD();

public:

	void Project(GRID& grid, const Vec3* velocity_arr, FLT* divergence_arr, FLT** pressure_arr, FLT** pressure_arr_temp)
	{
		#pragma omp parallel for
		for(int p = 0; p < grid.ijk_res_; p++) 
		{
			(*pressure_arr)[p] = (FLT)0;
			(*pressure_arr_temp)[p] = (FLT)0;
		}	
	}

	void ComputeDivergence(GRID& grid, const Vec3* velocity_arr, FLT* divergence_arr)
	{
		#pragma omp parallel for
		for(int p = 0; p < grid.ijk_res_; p++) 
		{
			int i, j, k;
			grid.Index1Dto3D(p, i, j, k);

			FLT div = (FLT)0;
				
			if(grid.IsGhostCell(i, j, k) == false)
			{
				div += (velocity_arr[p+1].x - velocity_arr[p-1].x)*grid.one_over_dx_*(FLT)0.5;
				div += (velocity_arr[p+grid.i_res_].y - velocity_arr[p-grid.i_res_].y)*grid.one_over_dy_*(FLT)0.5;
				div += (velocity_arr[p+grid.ij_res_].z - velocity_arr[p+grid.ij_res_].z)*grid.one_over_dz_*(FLT)0.5;

				divergence_arr[p] = div;
			}
			else
			{
				divergence_arr[p] = (FLT)0;
			}
		}	
	}

	void ComputePressureJacobi(GRID& grid, FLT* divergence_arr, FLT** pressure_arr, FLT** pressure_arr_ghost)
	{
		FLT dxdydz = POW2(grid.dx_*grid.dy_*grid.dz_);
		FLT dxdy   = POW2(grid.dx_*grid.dy_);
		FLT dydz   = POW2(grid.dy_*grid.dz_);
		FLT dxdz   = POW2(grid.dx_*grid.dz_);

		FLT one_over_h = (dxdy + dydz + dxdz)*(FLT)0.5;

		#pragma omp parallel for
		for(int p = 0; p < grid.ijk_res_; p++) 
		{
			int i, j, k;
			grid.Index1Dto3D(p, i, j, k);
				
			if(grid.IsGhostCell(i, j, k) == false)
			{
				FLT p_i1 = *pressure_arr[p+1];
				FLT p_i0 = *pressure_arr[p-1];

				FLT p_j1 = *pressure_arr[p+grid.i_res_];
				FLT p_j0 = *pressure_arr[p-grid.i_res_];

				FLT p_k1 = *pressure_arr[p+grid.ij_res_];
				FLT p_k0 = *pressure_arr[p-grid.ij_res_];

				(*pressure_arr_ghost)[p] = (dydz*(p_i1+p_i0) + dxdz*(p_j1+p_j0) + dxdy*(p_k1+p_k0) - divergence_arr[p]*dxdydz) * one_over_h;
			}
			else
			{
				(*pressure_arr_ghost)[p] = (FLT)0;
			}
		}

		FLT* temp = *pressure_arr_ghost;
		*pressure_arr_ghost = *pressure_arr;
		*pressure_arr = temp;
	}
};