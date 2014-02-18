
#pragma once

#include "KERNEL_FUNCTIONS.h"
#include "PARTICLE_MANAGER_3D.h"
#include "WALL_CONDITIONS.h"
#include "MATRIX3_T.h"

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

	T normal_stress_coef_;
	T shear_stress_coef_;

	Vec3T* pts_position_arr_;
	Vec3T* pts_velocity_arr_;
	Vec3T* pts_force_arr_;
	Vec3T* pts_grid_vel_arr_;
	T*     pts_density_arr_;
	Mat3T* pts_tensor_arr_;

	T*     density_field_;
	Vec3T* velocity_field_;
	Vec3T* force_field_;

	Vec3T gravity_;

	WALL_CONDITIONS wall_conditions_;

public:

	MPM_FLUID_SOLVER() : density_field_(0), velocity_field_(0), force_field_(0), mass_(1), rest_density_(1), normal_stress_coef_((T)0), shear_stress_coef_((T)0), smoothing_((T)0.3)
	{}

	~MPM_FLUID_SOLVER()
	{
		if (density_field_)  delete[] density_field_;
		if (velocity_field_) delete[] velocity_field_;
		if (force_field_)    delete[] force_field_;
	}

	void Initialize(const Vec3T min, const Vec3T max, const int i_res, const int j_res, const int k_res, const int ghost_width, const int num_pts_res);

	void AdvanceTimeStep(const T spf, const int steps)
	{
		const T dt = spf / (T)steps;

		for (int i = 0; i < steps; i++)
		{			
			RasterizeParticlesDensityAndVelocityToGrid();

			ComputeParticleDenistyFromGrid();

			ComputeStressTensors();

			ComputeGridForces();

			UpdateParticleAndGridVelocity(dt);

			AdvectParticles(dt);			

			SourceParticles();

			particle_manager_.RebuildParticleDataStructure();
		}
	}

	void RasterizeParticlesDensityAndVelocityToGrid();

	void ComputeParticleDenistyFromGrid();

	void ComputeStressTensors();

	void ComputeGridForces();

	void UpdateParticleAndGridVelocity(const T dt);

	void AdvectParticles(const T dt);

	void SourceParticles()
	{
		SourceFormSphere(Vec3T(0.2, 0.7, 0.5), Vec3T( 2, 0.0, 0.0), 0.05, 500);
		SourceFormSphere(Vec3T(0.8, 0.7, 0.5), Vec3T(-2, 0.0, 0.0), 0.05, 500);
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

	T ComputePressure(const T density)
	{
//		return (stiffness_*rest_density_) * ((density / rest_density_) - (T)1);
		return (stiffness_*rest_density_ / (T)3) * (POW3(density / rest_density_) - (T)1);
	}

	void ComputeStrainRate(const Vec3T& pos, Mat3T& tensor)
	{
		int i, j, k;
		grid_.CellCenterIndex(pos, i, j, k);

		T dudx((T)0), dudy((T)0), dudz((T)0), dvdx((T)0), dvdy((T)0), dvdz((T)0), dwdx((T)0), dwdy((T)0), dwdz((T)0);

		BEGIN_STENCIL_LOOP(grid_, i, j, k, l, m, n)
		{
			int s_ix = grid_.Index3Dto1D(l, m, n);

			Vec3T cell_pos = grid_.CellCenterPosition(l, m, n);
			Vec3T cell_vel = velocity_field_[s_ix];

			Vec3T grad = QuadBSplineKernelGradient(cell_pos-pos, grid_.one_over_dx_, grid_.one_over_dy_, grid_.one_over_dz_);

			dudx += cell_vel.x * grad.x;
			dudy += cell_vel.x * grad.y;
			dudz += cell_vel.x * grad.z;

			dvdx += cell_vel.y * grad.x;
			dvdy += cell_vel.y * grad.y;
			dvdz += cell_vel.y * grad.z;

			dwdx += cell_vel.z * grad.x;
			dwdy += cell_vel.z * grad.y;
			dwdz += cell_vel.z * grad.z;
		}
		END_STENCIL_LOOP;

		T D00 = dudx;
		T D11 = dvdy;
		T D22 = dwdz;
		T D01 = (dudy+dvdx)*(T)0.5;
		T D02 = (dudz+dwdx)*(T)0.5;
		T D12 = (dvdz+dwdy)*(T)0.5;

		tensor.SetRow(0, Vec3T(D00, D01, D02));
		tensor.SetRow(1, Vec3T(D01, D11, D12));
		tensor.SetRow(2, Vec3T(D02, D12, D22));		
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
