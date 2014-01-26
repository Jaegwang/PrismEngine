
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

	T stiffness_;

	T*     density_field_;
	Vec3T* force_field_;

public:

	MPM_FLUID_SOLVER() : density_field_(0), force_field_(0), mass_(1), rest_density_(1)
	{}

	~MPM_FLUID_SOLVER()
	{
		if (density_field_) delete[] density_field_;
		if (force_field_) delete[] force_field_;
	}

	void Initialize(const Vec3T min, const Vec3T max, const int i_res, const int j_res, const int k_res, const int ghost_width, const int num_pts_res)
	{
		grid_.Initialize(min, max, i_res, j_res, k_res, ghost_width);

		particle_manager_.Initialize(grid_, num_pts_res);		

		mass_ = 1;
		rest_density_ = 10;
		stiffness_ = 0.00001;

		density_field_ = new T[grid_.ijk_res_];
		force_field_ = new Vec3T[grid_.ijk_res_];

		memset((void*)density_field_, 0, sizeof(T)*grid_.ijk_res_);
		memset((void*)force_field_, 0, sizeof(Vec3T)*grid_.ijk_res_);

		SourceFormSphere(Vec3T(0.5, 0.5, 0.5), Vec3T(0.0, 0.0, 0.0), 0.1, 1000);
	}

	void AdvanceTimeStep(const T spf, const int steps)
	{
		const T dt = spf / (T)steps;

		for (int i = 0; i < steps; i++)
		{

			AdvectParticles(dt);

			RasterizeParticleDensityToGrid();
			ComputeParticleDenistyFromGrid();

			RasterizeParticForceToGrid();

//			ApplyForce(dt);
		}
	}

	void AdvectParticles(const T dt)
	{
		if(particle_manager_.num_of_pts_ == 0) return;

		array_view<Vec3T, 1> pos_arr_view(particle_manager_.num_of_pts_, particle_manager_.position_array_);
		array_view<Vec3T, 1> vel_arr_view(particle_manager_.num_of_pts_, particle_manager_.velocity_array_);

		array_view<T, 1> dt_view(1, &(T)dt);

		parallel_for_each(pos_arr_view.extent, [=](index<1> idx) restrict(amp)
		{
			pos_arr_view[idx[0]] += vel_arr_view[idx[0]] * dt_view[0];
		});

		pos_arr_view.synchronize();
		particle_manager_.RebuildParticleDataStructure();
	}

	void SourceFormSphere(const Vec3T& pos, const Vec3T& vel, const T rad, const int num)
	{	
		Vec3T* pos_arr = particle_manager_.position_array_;
		Vec3T* vel_arr = particle_manager_.velocity_array_;

		int num_pts = 0;
		while (num_pts < num)
		{
			Vec3T new_pos = pos + RandomVector()*rad;

			T dist = (new_pos - pos).Magnitude();
			if (dist <= rad)
			{
				num_pts++;
				int ix = particle_manager_.AddParticle();
				
				pos_arr[ix] = new_pos;
				vel_arr[ix] = vel;
			}
		}	
	}

	void RasterizeParticleDensityToGrid()
	{
		if (particle_manager_.num_of_pts_ == 0) return;

		array_view<T, 1> density_field_view(grid_.ijk_res_, density_field_);

		array_view<GRID_UNIFORM_3D, 1> grid_view(1, &grid_);

		array_view<T, 1> pts_mass(1, &mass_);
		array_view<T, 1> rest_density(1, &rest_density_);

		array_view<int, 1> num_pts_view(grid_.ijk_res_, particle_manager_.num_pts_cell_);
		array_view<int, 1> start_idx_view(grid_.ijk_res_, particle_manager_.start_idx_cell_);

		array_view<Vec3T, 1> position_arr_view(particle_manager_.num_of_pts_, particle_manager_.position_array_);

		concurrency::extent<1> ext(grid_.ijk_res_);
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

			int start_l, start_m, start_n, end_l, end_m, end_n;
			grid_view[0].StartEndIndices(i, j, k, start_l, start_m, start_n, end_l, end_m, end_n);

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

			density_field_view[ix] = mass_weighted;
		});

		density_field_view.synchronize();
	}


	void ComputeParticleDenistyFromGrid()
	{
		if (particle_manager_.num_of_pts_ == 0) return;

		array_view<Vec3T, 1> pos_arr_view(particle_manager_.num_of_pts_, particle_manager_.position_array_);
		array_view<T, 1> den_field_view(grid_.ijk_res_, density_field_);
		array_view<T, 1> den_arr_view(particle_manager_.num_of_pts_, particle_manager_.density_array_);

		array_view<GRID_UNIFORM_3D, 1> grid_view(1, &grid_);

		parallel_for_each(pos_arr_view.extent, [=](index<1> idx) restrict(amp)
		{
			int i, j, k;
			grid_view[0].CellCenterIndex(pos_arr_view[idx[0]], i, j, k);

			int start_l, start_m, start_n, end_l, end_m, end_n;
			grid_view[0].StartEndIndices(i, j, k, start_l, start_m, start_n, end_l, end_m, end_n);

			T density_weighted = (T)0;
			T weight = (T)0;

			for(int n = start_n; n <= end_n; n++) for (int m = start_m; m <= end_m; m++) for (int l = start_l; l <= end_l; l++)
			{
				int n_ix = grid_view[0].Index3Dto1D(l, m, n);

				Vec3T cell_center = grid_view[0].CellCenterPosition(l, m, n);
				T density_pts = den_field_view[n_ix];

				Vec3T deviation = cell_center - pos_arr_view[idx[0]];
				
				T w = MPMSplineKernel(deviation, grid_view[0].one_over_dx_, grid_view[0].one_over_dy_, grid_view[0].one_over_dz_);

				density_weighted += density_pts * w;
				weight += w;
			}

			if (weight > FLT_EPSILON)
			{
				den_arr_view[idx[0]] = density_weighted / weight;
			}
			else
			{
				den_arr_view[idx[0]] = (T)0;
			}
		});
	}


	void RasterizeParticForceToGrid()
	{
		if (particle_manager_.num_of_pts_ == 0) return;

		array_view<T, 1> density_field_view(grid_.ijk_res_, density_field_);
		array_view<Vec3T, 1> force_field_view(grid_.ijk_res_, force_field_);

		array_view<GRID_UNIFORM_3D, 1> grid_view(1, &grid_);

		array_view<T, 1> mass_view(1, &mass_);
		array_view<T, 1> rest_density_view(1, &rest_density_);
		array_view<T, 1> stiffness_view(1, &stiffness_);

		array_view<int, 1> num_pts_view(grid_.ijk_res_, particle_manager_.num_pts_cell_);
		array_view<int, 1> start_idx_view(grid_.ijk_res_, particle_manager_.start_idx_cell_);

		array_view<Vec3T, 1> position_arr_view(particle_manager_.num_of_pts_, particle_manager_.position_array_);
		array_view<T, 1> density_arr_view(particle_manager_.num_of_pts_, particle_manager_.density_array_);

		concurrency::extent<1> ext(grid_.ijk_res_);
		parallel_for_each(ext, [=](index<1> idx) restrict(amp)
		{
			int i, j, k;
			grid_view[0].Index1Dto3D(idx[0], i, j, k);
			
			int ix = grid_view[0].Index3Dto1D(i, j, k);

			if (grid_view[0].IsGhostCell(i, j, k) == true) return;

			Vec3T cell_center = grid_view[0].CellCenterPosition(i, j, k);

			T rho_i = density_field_view[ix];
			T p_i = stiffness_view[0] * (rho_i - rest_density_view[0]);

			Vec3T pressure(0, 0, 0);
			
			int start_l, start_m, start_n, end_l, end_m, end_n;
			grid_view[0].StartEndIndices(i, j, k, start_l, start_m, start_n, end_l, end_m, end_n);

			

			for (int n = start_n; n <= end_n; n++) for (int m = start_m; m <= end_m; m++) for (int l = start_l; l <= end_l; l++)
			{
				int n_ix = grid_view[0].Index3Dto1D(l, m, n);

				int num = num_pts_view[n_ix];
				int b_ix = start_idx_view[n_ix];

				for (int p = 0; p < num; p++)
				{
					const Vec3T& pos = position_arr_view[b_ix + p];

					Vec3T w = MPMSplineKernelGradient(pos - cell_center, grid_view[0].one_over_dx_, grid_view[0].one_over_dy_, grid_view[0].one_over_dz_);

					T rho_j = density_arr_view[b_ix + p];
					T p_j = stiffness_view[0] * (rho_j - rest_density_view[0]);

					pressure += mass_view[0] * (p_i + p_j) / (rho_j + rho_j) * w;
				}
			}

	//		force_field_view[ix] = -pressure / rho_i;

			
			
		});

//		force_field_view.synchronize();
	}

	void ApplyForce(const T dt)
	{
		if (particle_manager_.num_of_pts_ == 0) return;

		array_view<Vec3T, 1> pos_arr_view(particle_manager_.num_of_pts_, particle_manager_.position_array_);
		array_view<Vec3T, 1> vel_arr_view(particle_manager_.num_of_pts_, particle_manager_.velocity_array_);

		array_view<T, 1> den_arr_view(particle_manager_.num_of_pts_, particle_manager_.density_array_);
		array_view<T, 1> den_field_view(grid_.ijk_res_, density_field_);

		array_view<T, 1> stiffness_view(1, &stiffness_);
		array_view<T, 1> mass_view(1, &mass_);
		array_view<T, 1> rest_density_view(1, &rest_density_);
		array_view<GRID_UNIFORM_3D, 1> grid_view(1, &grid_);

		parallel_for_each(pos_arr_view.extent, [=](index<1> idx) restrict(amp)
		{
			int i, j, k;
			grid_view[0].CellCenterIndex(pos_arr_view[idx[0]], i, j, k);

			int start_l, start_m, start_n, end_l, end_m, end_n;
			grid_view[0].StartEndIndices(i, j, k, start_l, start_m, start_n, end_l, end_m, end_n);

			Vec3T pressure((T)0, (T)0, (T)0);
			
			T d_i = den_arr_view[idx[0]];
			T p_i = stiffness_view[0] * (den_arr_view[idx[0]] - rest_density_view[0]);

			for (int n = start_n; n <= end_n; n++) for (int m = start_m; m <= end_m; m++) for (int l = start_l; l <= end_l; l++)
			{
				int n_ix = grid_view[0].Index3Dto1D(l, m, n);

				Vec3T cell_center = grid_view[0].CellCenterPosition(l, m, n);

				T d_j = den_field_view[n_ix];

				if (d_j > FLT_EPSILON)
				{
					T p_j = stiffness_view[0] * (d_j - rest_density_view[0]);

//					Vec3T r_ij = cell_center - pos_arr_view[idx[0]];
					Vec3T r_ij = pos_arr_view[idx[0]] - cell_center;

					Vec3T w = MPMSplineKernelGradient(r_ij, grid_view[0].one_over_dx_, grid_view[0].one_over_dy_, grid_view[0].one_over_dz_);

					pressure += (w * mass_view[0] * (p_i + p_j) / (d_j + d_j)) *(T)0.01;
//					pressure += w * mass_view[0] * (p_j) / (d_j);
				}
			}

			if (d_i > FLT_EPSILON)
			{
				Vec3T acc = pressure / d_i;
				vel_arr_view[idx[0]] += acc * dt;
			}
		});

		vel_arr_view.synchronize();
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