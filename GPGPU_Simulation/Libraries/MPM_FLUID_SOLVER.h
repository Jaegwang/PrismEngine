
#pragma once

#include "KERNEL_FUNCTIONS.h"
#include "PARTICLE_MANAGER_3D.h"
#include "WALL_CONDITIONS.h"

#include "GL\glut.h"

class MPM_FLUID_SOLVER
{
public:

	GRID_UNIFORM_3D grid_;

	PARTICLE_MANAGER_3D particle_manager_;

	T mass_;
	T rest_density_;

	T stiffness_;
	T smoothing_;

	T*     density_field_;
	Vec3T* velocity_field_;
	Vec3T* force_field_;

	Vec3T gravity_;

	WALL_CONDITIONS wall_conditions_;


public:

	MPM_FLUID_SOLVER() : density_field_(0), velocity_field_(0), force_field_(0), mass_(1), rest_density_(1), smoothing_((T)0.2)
	{}

	~MPM_FLUID_SOLVER()
	{
		if (density_field_) delete[] density_field_;
		if (velocity_field_) delete[] velocity_field_;
		if (force_field_) delete[] force_field_;
	}

	void Initialize(const Vec3T min, const Vec3T max, const int i_res, const int j_res, const int k_res, const int ghost_width, const int num_pts_res)
	{
		grid_.Initialize(min, max, i_res, j_res, k_res, ghost_width);

		particle_manager_.Initialize(grid_, num_pts_res);		

		wall_conditions_.Initialize(grid_);

		mass_ = 1;
		rest_density_ = 1;
		stiffness_ = 3;

		density_field_  = new T[grid_.ijk_res_];
		velocity_field_ = new Vec3T[grid_.ijk_res_];
		force_field_    = new Vec3T[grid_.ijk_res_];

		gravity_ = Vec3T(0, -5, 0);

		memset((void*)density_field_, 0, sizeof(T)*grid_.ijk_res_);
		memset((void*)force_field_, 0, sizeof(Vec3T)*grid_.ijk_res_);
		memset((void*)velocity_field_, 0, sizeof(Vec3T)*grid_.ijk_res_);
	}

	void AdvanceTimeStep(const T spf, const int steps)
	{
		const T dt = spf / (T)steps;

		for (int i = 0; i < steps; i++)
		{			
			ApplyExternalForce(dt);

			RasterizeParticlesDensityAndVelocityToGrid();
			ComputeParticleDenistyFromGrid();

			RasterizeParticlesForceToGrid();
			UpdateParticleAndGridVelocity(dt);

			AdvectParticles(dt);

			SourceFormSphere(Vec3T(0.2, 0.7, 0.5), Vec3T( 2, 0.0, 0.0), 0.05, 500);
			SourceFormSphere(Vec3T(0.8, 0.7, 0.5), Vec3T(-2, 0.0, 0.0), 0.05, 500);

			particle_manager_.RebuildParticleDataStructure();
		}
	}

	void ApplyExternalForce(const T dt)
	{
		if(particle_manager_.num_of_pts_ == 0) return;

		BEGIN_CPU_THREADS_1D(particle_manager_.num_of_pts_)
		{
			for (int p = ix_begin; p <= ix_end; p++)
			{
				Vec3T& pts_pos = particle_manager_.position_array_[p];
				Vec3T& pts_vel = particle_manager_.velocity_array_[p];

				pts_vel += gravity_ * dt;				
			}
		}
		END_CPU_THREADS_1D;
	}

	void AdvectParticles(const T dt)
	{
		if (particle_manager_.num_of_pts_ == 0) return;

		BEGIN_CPU_THREADS_1D(particle_manager_.num_of_pts_)
		{
			for (int p = ix_begin; p <= ix_end; p++)
			{
				Vec3T& pts_pos = particle_manager_.position_array_[p];
				Vec3T& adv_vel = particle_manager_.adv_vel_array_[p];

				Vec3T& pts_vel = particle_manager_.velocity_array_[p];

				pts_pos += adv_vel * dt;

				wall_conditions_.ClampPositionAndVelocity(pts_pos, pts_vel);
			}
		}
		END_CPU_THREADS_1D;
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

	void RasterizeParticlesDensityAndVelocityToGrid()
	{
		if (particle_manager_.num_of_pts_ == 0) return;

		BEGIN_CPU_THREADS_1D(grid_.ijk_res_)
		{
			for (int p = ix_begin; p <= ix_end; p++)
			{
				int i, j, k;
				grid_.Index1Dto3D(p, i, j, k);

				int g_ix = grid_.Index3Dto1D(i, j, k);

				if (grid_.IsGhostCell(i, j, k) == true)
				{
					density_field_[g_ix] = rest_density_;
					velocity_field_[g_ix] = Vec3T();
					continue;
				}

				Vec3T cell_center = grid_.CellCenterPosition(i, j, k);

				Vec3T vel_weighted = Vec3T();
				T mass_weighted = (T)0;

				int start_l, start_m, start_n, end_l, end_m, end_n;
				grid_.StartEndIndices(i, j, k, start_l, start_m, start_n, end_l, end_m, end_n);

				for (int n = start_n; n <= end_n; n++) for (int m = start_m; m <= end_m; m++) for (int l = start_l; l <= end_l; l++)
				{
					int n_ix = grid_.Index3Dto1D(l, m, n);

					int num = particle_manager_.num_pts_cell_[n_ix];
					int b_ix = particle_manager_.start_idx_cell_[n_ix];

					for (int p = 0; p < num; p++)
					{
						const Vec3T& pos = particle_manager_.position_array_[b_ix + p];
						const Vec3T& vel = particle_manager_.velocity_array_[b_ix + p];

						T w = MPMSplineKernel(pos - cell_center, grid_.one_over_dx_, grid_.one_over_dy_, grid_.one_over_dz_);

						mass_weighted += mass_ * w;
						vel_weighted += vel * mass_ * w;
					}
				}

				density_field_[g_ix] = mass_weighted;

				if (mass_weighted > FLT_EPSILON)
				{
					
					velocity_field_[g_ix] = vel_weighted / mass_weighted;
				}
				else
				{
					velocity_field_[g_ix] = Vec3T();
				}
			}
		}
		END_CPU_THREADS_1D;
	}

	void ComputeParticleDenistyFromGrid()
	{
		if (particle_manager_.num_of_pts_ == 0) return;

		BEGIN_CPU_THREADS_1D(particle_manager_.num_of_pts_)
		{
			for (int p = ix_begin; p <= ix_end; p++)
			{
				const Vec3T& pts_pos = particle_manager_.position_array_[p];
				int i, j, k;
				grid_.CellCenterIndex(pts_pos, i, j, k);

				int start_l, start_m, start_n, end_l, end_m, end_n;
				grid_.StartEndIndices(i, j, k, start_l, start_m, start_n, end_l, end_m, end_n);

				T density_weighted = (T)0;

				for (int n = start_n; n <= end_n; n++) for (int m = start_m; m <= end_m; m++) for (int l = start_l; l <= end_l; l++)
				{
					int n_ix = grid_.Index3Dto1D(l, m, n);

					Vec3T cell_center = grid_.CellCenterPosition(l, m, n);
					T density_g = density_field_[n_ix];

					Vec3T deviation = cell_center - pts_pos;

					T w = MPMSplineKernel(deviation, grid_.one_over_dx_, grid_.one_over_dy_, grid_.one_over_dz_);

					density_weighted += density_g * w;
				}

				particle_manager_.density_array_[p] = density_weighted;
			}
		}
		END_CPU_THREADS_1D;
	}

	T ComputePressure(const T density)
	{
//		return (stiffness_*rest_density_) * ((density / rest_density_) - (T)1);
		return (stiffness_*rest_density_ / (T)3) * (POW3(density / rest_density_) - (T)1);
	}

	void RasterizeParticlesForceToGrid()
	{
		if (particle_manager_.num_of_pts_ == 0) return;

		BEGIN_CPU_THREADS_1D(grid_.ijk_res_)
		{
			for (int p = ix_begin; p <= ix_end; p++)
			{
				int i, j, k;
				grid_.Index1Dto3D(p, i, j, k);

				int g_ix = grid_.Index3Dto1D(i, j, k);

				if (grid_.IsGhostCell(i, j, k) == true) continue;

				Vec3T cell_center = grid_.CellCenterPosition(i, j, k);

				T rho_i = density_field_[g_ix];
				T p_i = ComputePressure(rho_i);

				Vec3T pressure((T)0,(T)0,(T)0);

				int start_l, start_m, start_n, end_l, end_m, end_n;
				grid_.StartEndIndices(i, j, k, start_l, start_m, start_n, end_l, end_m, end_n);

				for (int n = start_n; n <= end_n; n++) for (int m = start_m; m <= end_m; m++) for (int l = start_l; l <= end_l; l++)
				{
					int n_ix = grid_.Index3Dto1D(l, m, n);

					int num = particle_manager_.num_pts_cell_[n_ix];
					int b_ix = particle_manager_.start_idx_cell_[n_ix];

					for (int p = 0; p < num; p++)
					{
						const Vec3T& pos = particle_manager_.position_array_[b_ix + p];
						const Vec3T& vel = particle_manager_.velocity_array_[b_ix + p];

						T rho_j = particle_manager_.density_array_[b_ix + p];
						T p_j = ComputePressure(rho_j);

						Vec3T w_grad = MPMSplineKernelGradient(cell_center-pos, grid_.one_over_dx_, grid_.one_over_dy_, grid_.one_over_dz_);

						pressure += w_grad * mass_ * (p_i + p_j) / (rho_j + rho_j);
					}
				}

				force_field_[g_ix] = -pressure;
			}
		}
		END_CPU_THREADS_1D;
	}

	void UpdateParticleAndGridVelocity(const T dt)
	{
		if (particle_manager_.num_of_pts_ == 0) return;

		BEGIN_CPU_THREADS_1D(grid_.ijk_res_)
		{
			for (int p = ix_begin; p <= ix_end; p++)
			{
				int i, j, k;
				grid_.Index1Dto3D(p, i, j, k);

				int g_ix = grid_.Index3Dto1D(i, j, k);

				const T density_g = density_field_[g_ix];
				const Vec3T force_g = force_field_[g_ix];

				if (density_g > FLT_EPSILON) velocity_field_[g_ix] += force_g / density_g*dt;
			}
		}
		END_CPU_THREADS_1D;

		BEGIN_CPU_THREADS_1D(particle_manager_.num_of_pts_)
		{
			for (int p = ix_begin; p <= ix_end; p++)
			{
				const Vec3T& pts_pos = particle_manager_.position_array_[p];
				int i, j, k;
				grid_.CellCenterIndex(pts_pos, i, j, k);

				int start_l, start_m, start_n, end_l, end_m, end_n;
				grid_.StartEndIndices(i, j, k, start_l, start_m, start_n, end_l, end_m, end_n);

				Vec3T force_weighted = Vec3T();
				Vec3T vel_weighted = Vec3T();

				for (int n = start_n; n <= end_n; n++) for (int m = start_m; m <= end_m; m++) for (int l = start_l; l <= end_l; l++)
				{
					int n_ix = grid_.Index3Dto1D(l, m, n);

					Vec3T cell_center = grid_.CellCenterPosition(l, m, n);

					Vec3T velocity_g = velocity_field_[n_ix];
					Vec3T force_g = force_field_[n_ix];
					T density_g = density_field_[n_ix];

					if (density_g > FLT_EPSILON)
					{
						Vec3T deviation = cell_center - pts_pos;
						T w = MPMSplineKernel(deviation, grid_.one_over_dx_, grid_.one_over_dy_, grid_.one_over_dz_);

						force_weighted += force_g*w;
						vel_weighted += velocity_g*w;
					}
				}

				Vec3T& pts_vel = particle_manager_.velocity_array_[p];
				Vec3T& adv_vel = particle_manager_.adv_vel_array_[p];

				const T density_p = particle_manager_.density_array_[p];

				adv_vel = vel_weighted;

				pts_vel += force_weighted / density_p * dt;
				pts_vel = pts_vel*((T)1-smoothing_) + vel_weighted*smoothing_;
			}
		}
		END_CPU_THREADS_1D;
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
