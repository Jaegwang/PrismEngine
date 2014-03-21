
#include "MPM_FLUID_SOLVER.h"


void MPM_FLUID_SOLVER::Initialize(const Vec3 min, const Vec3 max, const int i_res, const int j_res, const int k_res, const int ghost_width, const int num_pts_res)
{
	grid_.Initialize(min, max, i_res, j_res, k_res, ghost_width);

	vector_temp_arr_  = new Vec3[num_pts_res]; 

	pts_position_arr_ = new Vec3[num_pts_res];
	pts_velocity_arr_ = new Vec3[num_pts_res];
	pts_force_arr_    = new Vec3[num_pts_res];
	pts_grid_vel_arr_ = new Vec3[num_pts_res];
	pts_density_arr_  = new FLT [num_pts_res];
	pts_tensor_arr_   = new Mat3[num_pts_res];

	particle_manager_.Initialize(grid_, &pts_position_arr_, &pts_velocity_arr_, num_pts_res);

	wall_conditions_.Initialize(grid_);

	cube_object_.InitializeCube(Vec3(0.4, 0.3, 0.4),Vec3(0.6, 0.5, 0.6), grid_.dx_, grid_.dy_, grid_.dz_);
	particle_world_.Initialize(grid_);

	particle_world_.RasterizeParticles(cube_object_.position_array_, cube_object_.velocity_array_, (FLT)1, cube_object_.pts_num_);

	mass_ = 1;
	rest_density_ = 3;
	stiffness_ = 0.03;

	normal_stress_coef_ = (FLT)0;
	shear_stress_coef_  = (FLT)0;

	density_field_  = new FLT[grid_.ijk_res_];
	velocity_field_ = new Vec3[grid_.ijk_res_];
	force_field_    = new Vec3[grid_.ijk_res_];

	gravity_ = Vec3(0, -3, 0);

//	memset((void*)pts_density_arr_, 0, sizeof(FLT)*num_pts_res);
	memset((void*)density_field_, 0, sizeof(FLT)*grid_.ijk_res_);
	memset((void*)force_field_, 0, sizeof(Vec3)*grid_.ijk_res_);
	memset((void*)velocity_field_, 0, sizeof(Vec3)*grid_.ijk_res_);

	SourceParticles();
	particle_manager_.RebuildParticleDataStructure();
}


void MPM_FLUID_SOLVER::RasterizeParticlesDensityAndVelocityToGrid()
{
	if(particle_manager_.num_of_pts_ == 0) return;

	#pragma omp parallel for
	for (int p = 0 ; p < grid_.ijk_res_ ; p++) 
	{
		int i, j, k;
		grid_.Index1Dto3D(p, i, j, k);
				
		if (grid_.IsGhostCell(i, j, k) == true)
		{
			density_field_[p] = rest_density_;
			velocity_field_[p] = Vec3();
			continue;
		}

		Vec3 cell_center = grid_.CellCenterPosition(i, j, k);

		Vec3 vel_weighted = Vec3();
		FLT mass_weighted = (FLT)0;

		BEGIN_STENCIL_LOOP(grid_, i,j,k, l,m,n)
		{
			int s_ix = grid_.Index3Dto1D(l, m, n);

			int num = particle_manager_.num_pts_cell_[s_ix];
			int b_ix = particle_manager_.start_idx_cell_[s_ix];

			for (int x = 0; x < num; x++)
			{
				const Vec3& pos = pts_position_arr_[b_ix + x];
				const Vec3& vel = pts_velocity_arr_[b_ix + x];

				FLT w = QuadBSplineKernel(pos-cell_center, grid_.one_over_dx_, grid_.one_over_dy_, grid_.one_over_dz_);

				mass_weighted += mass_ * w;
				vel_weighted += vel * mass_ * w;
			}
		}
		END_STENCIL_LOOP

		density_field_[p] = mass_weighted;
		velocity_field_[p] = vel_weighted / (mass_weighted+FLT_EPSILON);
	}
}

void MPM_FLUID_SOLVER::ComputeParticleDenistyFromGrid()
{
	if(particle_manager_.num_of_pts_ == 0) return;

	#pragma omp parallel for
	for(int p=0; p<particle_manager_.num_of_pts_; p++)
	{
		const Vec3& pts_pos = pts_position_arr_[p];

		int i, j, k;
		grid_.CellCenterIndex(pts_pos, i, j, k);

		FLT density_weighted = (FLT)0;

		BEGIN_STENCIL_LOOP(grid_, i, j, k, l, m, n)
		{
			int s_ix = grid_.Index3Dto1D(l, m, n);

			Vec3 cell_center = grid_.CellCenterPosition(l, m, n);
			FLT density_g = density_field_[s_ix];

			FLT w = QuadBSplineKernel(cell_center-pts_pos, grid_.one_over_dx_, grid_.one_over_dy_, grid_.one_over_dz_);

			density_weighted += density_g*w;
		}
		END_STENCIL_LOOP;

		pts_density_arr_[p] = density_weighted;
	}
}

void MPM_FLUID_SOLVER::ComputeStressTensors()
{
	if (particle_manager_.num_of_pts_ == 0) return;

	#pragma omp parallel for
	for(int p=0; p<particle_manager_.num_of_pts_; p++)
	{
		const Vec3& pts_pos = pts_position_arr_[p];
		const Vec3& pts_force = pts_force_arr_[p];
		const FLT& pts_density = pts_density_arr_[p];

		Mat3& pts_tensor = pts_tensor_arr_[p];

		FLT pressure = ComputePressure(pts_density);

		ComputeStrainTensor(pts_pos, pts_tensor);

		pts_tensor[0][0] = pts_tensor[0][0]*normal_stress_coef_ + pressure;
		pts_tensor[1][1] = pts_tensor[1][1]*normal_stress_coef_ + pressure;
		pts_tensor[2][2] = pts_tensor[2][2]*normal_stress_coef_ + pressure;

		pts_tensor[0][1] = pts_tensor[0][1]*shear_stress_coef_;
		pts_tensor[0][2] = pts_tensor[0][2]*shear_stress_coef_;
		pts_tensor[1][0] = pts_tensor[1][0]*shear_stress_coef_;
		pts_tensor[1][2] = pts_tensor[1][2]*shear_stress_coef_;
		pts_tensor[2][0] = pts_tensor[2][0]*shear_stress_coef_;
		pts_tensor[2][1] = pts_tensor[2][1]*shear_stress_coef_;
	}
}

void MPM_FLUID_SOLVER::ComputeGridForces()
{
	if (particle_manager_.num_of_pts_ == 0) return;

	#pragma omp parallel for
	for(int p=0; p<grid_.ijk_res_; p++)
	{
		int i, j, k;
		grid_.Index1Dto3D(p, i, j, k);

		if(grid_.IsGhostCell(i, j, k) == true)
		{
			force_field_[p] = Vec3();
			continue;
		}

		Vec3 cell_center = grid_.CellCenterPosition(i, j, k);

		Vec3 int_force((FLT)0, (FLT)0, (FLT)0);
		Vec3 ext_force((FLT)0, (FLT)0, (FLT)0);

		BEGIN_STENCIL_LOOP(grid_, i,j,k, l,m,n)
		{
			int s_ix = grid_.Index3Dto1D(l, m, n);
			int num = particle_manager_.num_pts_cell_[s_ix];
			int b_ix = particle_manager_.start_idx_cell_[s_ix];

			for (int p = 0; p < num; p++)
			{
				const Vec3& pos    = pts_position_arr_[b_ix + p];
				const Vec3& force  = pts_force_arr_[b_ix + p];
				const Mat3& tensor = pts_tensor_arr_[b_ix + p];
	
				FLT w = QuadBSplineKernel(cell_center-pos, grid_.one_over_dx_, grid_.one_over_dy_, grid_.one_over_dz_);
				Vec3 grad = QuadBSplineKernelGradient(cell_center-pos, grid_.one_over_dx_, grid_.one_over_dy_, grid_.one_over_dz_);

				int_force += tensor*grad;
				ext_force += w*force;
			}
		}
		END_STENCIL_LOOP;

		force_field_[p] = -int_force + ext_force;
	}

	Vec3 gravity_force = gravity_ * mass_;

	#pragma omp parallel for
	for(int p=0; p<particle_manager_.num_of_pts_; p++)
	{
		pts_force_arr_[p] = gravity_force;
	}
}

void MPM_FLUID_SOLVER::ComputeGridForces2()
{
	if (particle_manager_.num_of_pts_ == 0) return;

	#pragma omp parallel for
	for(int p=0; p<grid_.ijk_res_; p++)
	{
		int i, j, k;
		grid_.Index1Dto3D(p, i, j, k);

		if(grid_.IsGhostCell(i, j, k) == true)
		{
			force_field_[p] = Vec3();
			continue;
		}

		Vec3 cell_center = grid_.CellCenterPosition(i, j, k);

		Vec3 int_force((FLT)0, (FLT)0, (FLT)0);
		Vec3 ext_force((FLT)0, (FLT)0, (FLT)0);

		BEGIN_STENCIL_LOOP(grid_, i,j,k, l,m,n)
		{
			int s_ix = grid_.Index3Dto1D(l, m, n);
			int num = particle_manager_.num_pts_cell_[s_ix];
			int b_ix = particle_manager_.start_idx_cell_[s_ix];

			for (int p = 0; p < num; p++)
			{
				const Vec3& position = pts_position_arr_[b_ix + p];
				const Vec3& force = pts_force_arr_[b_ix + p]; 
				const FLT pressure = pts_density_arr_[b_ix + p];
	
				FLT w = QuadBSplineKernel(cell_center-position, grid_.one_over_dx_, grid_.one_over_dy_, grid_.one_over_dz_);
				Vec3 grad = QuadBSplineKernelGradient(cell_center-position, grid_.one_over_dx_, grid_.one_over_dy_, grid_.one_over_dz_);

				int_force += pressure*grad;
				ext_force += w*force;
			}
		}
		END_STENCIL_LOOP;

		force_field_[p] = -int_force + ext_force;
	}

	Vec3 gravity_force = gravity_ * mass_;

	#pragma omp parallel for
	for(int p=0; p<particle_manager_.num_of_pts_; p++)
	{
		pts_force_arr_[p] = gravity_force;
	}

}

void MPM_FLUID_SOLVER::RebuildParticleDataStructure()
{
	particle_manager_.RebuildParticleDataStructure();
	
	particle_manager_.RearrangeParticleData(&pts_position_arr_, &vector_temp_arr_);
	particle_manager_.RearrangeParticleData(&pts_velocity_arr_, &vector_temp_arr_);
	particle_manager_.RearrangeParticleData(&pts_force_arr_   , &vector_temp_arr_);
}

void MPM_FLUID_SOLVER::UpdateParticleAndGridVelocity(const FLT dt)
{
	if (particle_manager_.num_of_pts_ == 0) return;

	#pragma omp parallel for
	for(int p=0; p<grid_.ijk_res_; p++)
	{
		int i, j, k;
		grid_.Index1Dto3D(p, i, j, k);

		const FLT density_g = density_field_[p];
		const Vec3 force_g = force_field_[p];

		velocity_field_[p] += force_g / (density_g+(FLT)FLT_EPSILON) * dt;
	}

	#pragma omp parallel for
	for(int p=0; p<particle_manager_.num_of_pts_; p++)
	{
		const Vec3& pts_pos = pts_position_arr_[p];

		int i, j, k;
		grid_.CellCenterIndex(pts_pos, i, j, k);

		Vec3 force_weighted = Vec3();
		Vec3 vel_weighted = Vec3();

		BEGIN_STENCIL_LOOP(grid_, i, j, k, l, m, n)
		{
			int s_ix = grid_.Index3Dto1D(l, m, n);
			Vec3 cell_center = grid_.CellCenterPosition(l, m, n);

			Vec3 velocity_g = velocity_field_[s_ix];
			Vec3 force_g = force_field_[s_ix];

			FLT w = QuadBSplineKernel(cell_center-pts_pos, grid_.one_over_dx_, grid_.one_over_dy_, grid_.one_over_dz_);

			force_weighted += force_g*w;
			vel_weighted += velocity_g*w;
		}
		END_STENCIL_LOOP;

		Vec3& pts_vel = pts_velocity_arr_[p];
		Vec3& grid_vel = pts_grid_vel_arr_[p];

		const FLT density_p = pts_density_arr_[p];

		grid_vel = vel_weighted;

		pts_vel += force_weighted / density_p * dt;
		pts_vel = pts_vel*((FLT)1-smoothing_) + grid_vel*smoothing_;
	}
}


void MPM_FLUID_SOLVER::AdvectParticles(const FLT dt)
{
	if (particle_manager_.num_of_pts_ == 0) return;

	#pragma omp parallel for
	for(int p=0; p<particle_manager_.num_of_pts_; p++)
	{
		Vec3& pts_pos = pts_position_arr_[p];
		Vec3& grid_vel = pts_grid_vel_arr_[p];

		Vec3& pts_vel = pts_velocity_arr_[p];
		Vec3& pts_force = pts_force_arr_[p];

		pts_pos += grid_vel * dt;
	}
}

void MPM_FLUID_SOLVER::CouplingWithObjects(const FLT dt)
{
	if (particle_manager_.num_of_pts_ == 0) return;

	#pragma omp parallel for
	for(int p=0; p<particle_manager_.num_of_pts_; p++)
	{
		Vec3& pts_pos   = pts_position_arr_[p];
		Vec3& pts_vel   = pts_velocity_arr_[p];
		Vec3& pts_force = pts_force_arr_[p];

		wall_conditions_.ClampPositionAndVelocity(pts_pos, pts_vel, pts_force);

		FLT object_density = particle_world_.world_grid_.TriLinearInterpolate(pts_pos, particle_world_.density_field_);
		Vec3 object_normal = particle_world_.world_grid_.TriLinearInterpolate(pts_pos, particle_world_.normal_field_); 

		if(object_density > FLT_EPSILON && glm::dot(pts_vel, object_normal) < (FLT)0)
		{
			pts_force += glm::normalize(object_normal) * object_density * mass_ * (FLT)500;
		}		
	}
}