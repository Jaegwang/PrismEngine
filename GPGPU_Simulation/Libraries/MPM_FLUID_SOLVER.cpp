
#include "MPM_FLUID_SOLVER.h"


void MPM_FLUID_SOLVER::Initialize(const Vec3T min, const Vec3T max, const int i_res, const int j_res, const int k_res, const int ghost_width, const int num_pts_res)
{
	grid_.Initialize(min, max, i_res, j_res, k_res, ghost_width);

	vector_temp_arr_  = new Vec3T[num_pts_res]; 

	pts_position_arr_ = new Vec3T[num_pts_res];
	pts_velocity_arr_ = new Vec3T[num_pts_res];
	pts_force_arr_    = new Vec3T[num_pts_res];
	pts_grid_vel_arr_ = new Vec3T[num_pts_res];
	pts_density_arr_  = new T    [num_pts_res];
	pts_tensor_arr_   = new Mat3T[num_pts_res];

	particle_manager_.Initialize(grid_, &pts_position_arr_, &pts_velocity_arr_, &vector_temp_arr_, num_pts_res);

	wall_conditions_.Initialize(grid_);

	mass_ = 1;
	rest_density_ = 2;
	stiffness_ = 0.2;

	normal_stress_coef_ = (T)0;
	shear_stress_coef_ = (T)0;

	density_field_  = new T[grid_.ijk_res_];
	velocity_field_ = new Vec3T[grid_.ijk_res_];
	force_field_    = new Vec3T[grid_.ijk_res_];

	gravity_ = Vec3T(0, -5, 0);

//	memset((void*)pts_density_arr_, 0, sizeof(T)*num_pts_res);
	memset((void*)density_field_, 0, sizeof(T)*grid_.ijk_res_);
	memset((void*)force_field_, 0, sizeof(Vec3T)*grid_.ijk_res_);
	memset((void*)velocity_field_, 0, sizeof(Vec3T)*grid_.ijk_res_);

	SourceParticles();
	particle_manager_.RebuildParticleDataStructure();
}


void MPM_FLUID_SOLVER::RasterizeParticlesDensityAndVelocityToGrid()
{
	if(particle_manager_.num_of_pts_ == 0) return;

	BEGIN_CPU_THREADS(grid_.ijk_res_, p)
	{
		int i, j, k;
		grid_.Index1Dto3D(p, i, j, k);
				
		if (grid_.IsGhostCell(i, j, k) == true)
		{
			density_field_[p] = rest_density_;
			velocity_field_[p] = Vec3T();
			continue;
		}

		Vec3T cell_center = grid_.CellCenterPosition(i, j, k);

		Vec3T vel_weighted = Vec3T();
		T mass_weighted = (T)0;

		BEGIN_STENCIL_LOOP(grid_, i,j,k, l,m,n)
		{
			int s_ix = grid_.Index3Dto1D(l, m, n);

			int num = particle_manager_.num_pts_cell_[s_ix];
			int b_ix = particle_manager_.start_idx_cell_[s_ix];

			for (int x = 0; x < num; x++)
			{
				const Vec3T& pos = pts_position_arr_[b_ix + x];
				const Vec3T& vel = pts_velocity_arr_[b_ix + x];

				T w = QuadBSplineKernel(pos-cell_center, grid_.one_over_dx_, grid_.one_over_dy_, grid_.one_over_dz_);

				mass_weighted += mass_ * w;
				vel_weighted += vel * mass_ * w;
			}
		}
		END_STENCIL_LOOP

		density_field_[p] = mass_weighted;
		velocity_field_[p] = vel_weighted / (mass_weighted+FLT_EPSILON);
	}
	END_CPU_THREADS;
}

void MPM_FLUID_SOLVER::ComputeParticleDenistyFromGrid()
{
	if(particle_manager_.num_of_pts_ == 0) return;

	BEGIN_CPU_THREADS(particle_manager_.num_of_pts_, p)
	{
		const Vec3T& pts_pos = pts_position_arr_[p];

		int i, j, k;
		grid_.CellCenterIndex(pts_pos, i, j, k);

		T density_weighted = (T)0;

		BEGIN_STENCIL_LOOP(grid_, i, j, k, l, m, n)
		{
			int s_ix = grid_.Index3Dto1D(l, m, n);

			Vec3T cell_center = grid_.CellCenterPosition(l, m, n);
			T density_g = density_field_[s_ix];

			T w = QuadBSplineKernel(cell_center-pts_pos, grid_.one_over_dx_, grid_.one_over_dy_, grid_.one_over_dz_);

			density_weighted += density_g*w;
		}
		END_STENCIL_LOOP;

		pts_density_arr_[p] = density_weighted;
	}
	END_CPU_THREADS;
}

void MPM_FLUID_SOLVER::ComputeStressTensors()
{
	if (particle_manager_.num_of_pts_ == 0) return;

	BEGIN_CPU_THREADS(particle_manager_.num_of_pts_, p)
	{
		const Vec3T& pts_pos = pts_position_arr_[p];
		const Vec3T& pts_force = pts_force_arr_[p];
		const T& pts_density = pts_density_arr_[p];

		Mat3T& pts_tensor = pts_tensor_arr_[p];

		T pressure = ComputePressure(pts_density);

		ComputeStrainRate(pts_pos, pts_tensor);

		pts_tensor(0,0) = pts_tensor(0,0)*normal_stress_coef_ + pressure;
		pts_tensor(1,1) = pts_tensor(1,1)*normal_stress_coef_ + pressure;
		pts_tensor(2,2) = pts_tensor(2,2)*normal_stress_coef_ + pressure;

		pts_tensor(0,1) = pts_tensor(0,1)*shear_stress_coef_;
		pts_tensor(0,2) = pts_tensor(0,2)*shear_stress_coef_;
		pts_tensor(1,0) = pts_tensor(1,0)*shear_stress_coef_;
		pts_tensor(1,2) = pts_tensor(1,2)*shear_stress_coef_;
		pts_tensor(2,0) = pts_tensor(2,0)*shear_stress_coef_;
		pts_tensor(2,1) = pts_tensor(2,1)*shear_stress_coef_;
	}
	END_CPU_THREADS;
}

void MPM_FLUID_SOLVER::ComputeGridForces()
{
	if (particle_manager_.num_of_pts_ == 0) return;

	BEGIN_CPU_THREADS(grid_.ijk_res_, p)
	{
		int i, j, k;
		grid_.Index1Dto3D(p, i, j, k);

		if(grid_.IsGhostCell(i, j, k) == true)
		{
			force_field_[p] = Vec3T();
			continue;
		}

		Vec3T cell_center = grid_.CellCenterPosition(i, j, k);

		Vec3T int_force((T)0, (T)0, (T)0);
		Vec3T ext_force((T)0, (T)0, (T)0);

		BEGIN_STENCIL_LOOP(grid_, i,j,k, l,m,n)
		{
			int s_ix = grid_.Index3Dto1D(l, m, n);
			int num = particle_manager_.num_pts_cell_[s_ix];
			int b_ix = particle_manager_.start_idx_cell_[s_ix];

			for (int p = 0; p < num; p++)
			{
				const Vec3T& pos    = pts_position_arr_[b_ix + p];
				const Vec3T& force  = pts_force_arr_[b_ix + p];
				const Mat3T& tensor = pts_tensor_arr_[b_ix + p];
	
				T w = QuadBSplineKernel(cell_center-pos, grid_.one_over_dx_, grid_.one_over_dy_, grid_.one_over_dz_);
				Vec3T grad = QuadBSplineKernelGradient(cell_center-pos, grid_.one_over_dx_, grid_.one_over_dy_, grid_.one_over_dz_);

				int_force += tensor*grad;
				ext_force += w*force;
			}
		}
		END_STENCIL_LOOP;

		force_field_[p] = -int_force + ext_force;
	}
	END_CPU_THREADS;

	Vec3T gravity_force = gravity_ * mass_;

	BEGIN_CPU_THREADS(particle_manager_.num_of_pts_, p)
	{
		pts_force_arr_[p] = gravity_force;
	}
	END_CPU_THREADS;
}

void MPM_FLUID_SOLVER::RebuildParticleDataStructure()
{
	particle_manager_.RebuildParticleDataStructure();
	
	particle_manager_.RearrangeParticleData(&pts_force_arr_, &vector_temp_arr_);
}

void MPM_FLUID_SOLVER::UpdateParticleAndGridVelocity(const T dt)
{
	if (particle_manager_.num_of_pts_ == 0) return;

	BEGIN_CPU_THREADS(grid_.ijk_res_, p)
	{
		int i, j, k;
		grid_.Index1Dto3D(p, i, j, k);

		const T density_g = density_field_[p];
		const Vec3T force_g = force_field_[p];

		velocity_field_[p] += force_g / (density_g+(T)FLT_EPSILON) * dt;
	}
	END_CPU_THREADS;

	BEGIN_CPU_THREADS(particle_manager_.num_of_pts_, p)
	{
		const Vec3T& pts_pos = pts_position_arr_[p];

		int i, j, k;
		grid_.CellCenterIndex(pts_pos, i, j, k);

		Vec3T force_weighted = Vec3T();
		Vec3T vel_weighted = Vec3T();

		BEGIN_STENCIL_LOOP(grid_, i, j, k, l, m, n)
		{
			int s_ix = grid_.Index3Dto1D(l, m, n);
			Vec3T cell_center = grid_.CellCenterPosition(l, m, n);

			Vec3T velocity_g = velocity_field_[s_ix];
			Vec3T force_g = force_field_[s_ix];

			T w = QuadBSplineKernel(cell_center-pts_pos, grid_.one_over_dx_, grid_.one_over_dy_, grid_.one_over_dz_);

			force_weighted += force_g*w;
			vel_weighted += velocity_g*w;
		}
		END_STENCIL_LOOP;

		Vec3T& pts_vel = pts_velocity_arr_[p];
		Vec3T& grid_vel = pts_grid_vel_arr_[p];

		const T density_p = pts_density_arr_[p];

		grid_vel = vel_weighted;

		pts_vel += force_weighted / density_p * dt;
		pts_vel = pts_vel*((T)1-smoothing_) + grid_vel*smoothing_;
	}
	END_CPU_THREADS;
}


void MPM_FLUID_SOLVER::AdvectParticles(const T dt)
{
	if (particle_manager_.num_of_pts_ == 0) return;

	BEGIN_CPU_THREADS(particle_manager_.num_of_pts_, p)
	{
		Vec3T& pts_pos = pts_position_arr_[p];
		Vec3T& grid_vel = pts_grid_vel_arr_[p];

		Vec3T& pts_vel = pts_velocity_arr_[p];
		Vec3T& pts_force = pts_force_arr_[p];

		pts_pos += grid_vel * dt;

		wall_conditions_.ClampPositionAndVelocity(pts_pos, pts_vel, pts_force);
	}
	END_CPU_THREADS;
}