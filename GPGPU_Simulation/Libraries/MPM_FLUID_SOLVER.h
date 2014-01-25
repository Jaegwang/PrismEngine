
#pragma once

#include "KERNEL_FUNCTIONS.h"
#include "PARTICLE_MANAGER_3D.h"

#include "GL\glut.h"

class MPM_FLUID_SOLVER
{
public:

	GRID_UNIFORM_3D grid_;

	PARTICLE_MANAGER_3D particle_manager_;

	T mass_;
	T rest_density_;

	T* density_field_;

public:

	MPM_FLUID_SOLVER() : density_field_(0), mass_(1), rest_density_(1)
	{}

	~MPM_FLUID_SOLVER()
	{
		if (density_field_) delete[] density_field_;
	}

	void Initialize(const Vec3T min, const Vec3T max, const int i_res, const int j_res, const int k_res, const int ghost_width, const int num_pts_res)
	{
		grid_.Initialize(min, max, i_res, j_res, k_res, ghost_width);

		particle_manager_.Initialize(grid_, num_pts_res);		

		mass_ = 1;
		rest_density_ = 1;

		density_field_ = new T[grid_.ijk_res_];

		memset((void*)density_field_, 0, sizeof(T)*grid_.ijk_res_);

		
		particle_manager_.InitializeParticleArray();
		particle_manager_.RebuildParticleDataStructure();
	}

	void AdvanceTimeStep(const T spf, const int steps)
	{
		const T dt = spf / (T)steps;

		for (int i = 0; i < steps; i++)
		{


		}
	}

	void RasterizeDensityParticlesToGrid()
	{
		array_view<T, 1> density_field_view(grid_.ijk_res_, density_field_);

		array_view<GRID_UNIFORM_3D, 1> grid_view(1, &grid_);

		array_view<T, 1> pts_mass(1, &mass_);
		array_view<T, 1> rest_density(1, &rest_density_);

		array_view<int, 1> num_pts_view(grid_.ijk_res_, particle_manager_.num_pts_cell_);
		array_view<int, 1> start_idx_view(grid_.ijk_res_, particle_manager_.start_idx_cell_);

		array_view<Vec3T, 1> position_arr_view(particle_manager_.num_of_pts_, particle_manager_.position_array_);

		extent<1> ext(grid_.ijk_res_);

		parallel_for_each(ext, [=](index<1> idx) restrict(amp)
		{
			int i, j, k;
			grid_view[0].Index1Dto3D(idx[0], i, j, k);

			int ix = grid_view[0].Index3Dto1D(i, j, k);

			if (grid_view[0].IsGhostCell(i, j, k) == true)
			{
//				density_field_view[ix] = rest_density[0];
				density_field_view[ix] = (T)0;
				return;
			}
			
			Vec3T cell_center = grid_view[0].CellCenterPosition(i, j, k);

			T mass_weighted = (T)0;

			int start_l = i-1, start_m = j-1, start_n = k-1;
			int end_l = i+1, end_m = j+1, end_n = k+1;

			for(int n = start_n; n <= end_n; n++) for(int m = start_m; m <= end_m; m++) for(int l = start_l; l <= end_l; l++)
			{
				int n_ix = grid_view[0].Index3Dto1D(l, m, n);

				int num = num_pts_view[n_ix];
				int b_ix = start_idx_view[n_ix];

				for (int p = 0; p < num; p++)
				{
					const Vec3T& pos = position_arr_view[b_ix + p];

					T w = MPMSplineKernel(pos - cell_center, grid_view[0].one_over_dx_, grid_view[0].one_over_dy_, grid_view[0].one_over_dz_);

					mass_weighted += pts_mass[0] * w;
				}
			}

			if(mass_weighted > FLT_EPSILON) density_field_view[ix] = mass_weighted;
			else density_field_view[ix] = (T)0;

		});

		density_field_view.synchronize();
	}


	void RenderDensityField()
	{
		glDisable(GL_LIGHTING);

		glPushMatrix();
		glColor3f(0, 1, 0);

		glBegin(GL_POINTS);
		for (int c = 0; c < grid_.ijk_res_; c++)
		{
			int i, j, k;
			grid_.Index1Dto3D(c, i, j, k);

			Vec3T cell_center = grid_.CellCenterPosition(i, j, k);

			T density = density_field_[c];

			if (density > FLT_EPSILON)
			{
				glVertex3f(cell_center.x, cell_center.y, cell_center.z);
			}
		}
		glEnd();

		glPopMatrix();
	}
};