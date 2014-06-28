
#pragma once

#include "KERNEL_FUNCTIONS.h"
#include "PARTICLE_MANAGER_3D.h"
#include "WALL_CONDITIONS.h"
#include "MATH_CORE.h"

#include "PARTICLE_OBJECT.h"
#include "PARTICLE_WORLD_BUILDER.h"

#include "GL\glut.h"

class MPM_FLUID_SOLVER
{
public:

	GRID grid_;

	PARTICLE_MANAGER_3D particle_manager_;

	TS mass_;
	TS rest_density_;

	TS stiffness_;
	TS smoothing_;

	TS normal_stress_coef_;
	TS shear_stress_coef_;

	TV3* pts_position_arr_;
	TV3* pts_velocity_arr_;
	TV3* pts_force_arr_;
	TV3* pts_grid_vel_arr_;
	TS*  pts_density_arr_;
	TM3* pts_tensor_arr_;

	TV3* vector_temp_arr_;
	TV3* scalar_temp_arr_;

	TS*  density_field_;
	TV3* velocity_field_;
	TV3* force_field_;

	TV3 gravity_;

	WALL_CONDITIONS wall_conditions_;

	// Object properties
	PARTICLE_OBJECT cube_object_;
	PARTICLE_WORLD_BUILDER particle_world_;

public:

	MPM_FLUID_SOLVER() : density_field_(0), velocity_field_(0), force_field_(0), mass_(1), rest_density_(1), normal_stress_coef_((TS)0), shear_stress_coef_((TS)0), smoothing_((TS)0.3)
	{}

	~MPM_FLUID_SOLVER()
	{
		if (density_field_)  delete[] density_field_;
		if (velocity_field_) delete[] velocity_field_;
		if (force_field_)    delete[] force_field_;
	}

	void Initialize(const TV3 min, const TV3 max, const int i_res, const int j_res, const int k_res, const int ghost_width, const int num_pts_res);

	void AdvanceTimeStep(const TS spf, const int steps)
	{
		SourceParticles();
		RebuildParticleDataStructure();

		const TS dt = spf / (TS)steps;
		for (int i = 0; i < steps; i++)
		{			
			RasterizeParticlesDensityAndVelocityToGrid();

			ComputeParticleDenistyFromGrid();

			ComputeStressTensors();

			ComputeGridForces();

			UpdateParticleAndGridVelocity(dt);

			AdvectParticles(dt);			

			CouplingWithObjects(dt);

			RebuildParticleDataStructure();
		}
	}

	void RasterizeParticlesDensityAndVelocityToGrid();

	void RebuildParticleDataStructure();

	void ComputeParticleDenistyFromGrid();

	void ComputeStressTensors();

	void ComputeGridForces();

	void ComputeGridForces2();

	void UpdateParticleAndGridVelocity(const TS dt);

	void AdvectParticles(const TS dt);

	void CouplingWithObjects(const TS dt);

	void SourceParticles()
	{
		SourceFormSphere(TV3(0.2, 0.7, 0.5), TV3( 0.7, 0.0, 0.0), 0.05, 500);
		SourceFormSphere(TV3(0.8, 0.7, 0.5), TV3(-0.7, 0.0, 0.0), 0.05, 500);
	}

	void SourceFormSphere(const TV3& pos, const TV3& vel, const TS rad, const int num)
	{	
		TV3* pos_arr = pts_position_arr_;
		TV3* vel_arr = pts_velocity_arr_;

		int num_pts = 0;
		while (num_pts < num)
		{
			TV3 new_pos = pos + RandomVector()*rad;

			TS dist = glm::length((new_pos - pos));
			if (dist <= rad)
			{
				num_pts++;
				int ix = particle_manager_.AddParticle();
				
				pos_arr[ix] = new_pos;
				vel_arr[ix] = vel;
			}
		}
	}

	TS ComputePressure(const TS density)
	{
//		return (stiffness_*rest_density_) * ((density / rest_density_) - (TS)1);
		return (stiffness_*rest_density_ / (TS)3) * (POW3(density / rest_density_) - (TS)1);
	}

	void ComputeWeightStencil(const TV3& pos, const TV3& cell, TS wx[3],  TS wy[3],  TS wz[3])
	{
		wx[0] = QuadBSplineKernel(cell.x-grid_.dx_ - pos.x, grid_.one_over_dx_);	
		wx[1] = QuadBSplineKernel(cell.x           - pos.x, grid_.one_over_dx_);
		wx[2] = QuadBSplineKernel(cell.x+grid_.dx_ - pos.x, grid_.one_over_dx_);

		wy[0] = QuadBSplineKernel(cell.y-grid_.dx_ - pos.y, grid_.one_over_dx_);	
		wy[1] = QuadBSplineKernel(cell.y           - pos.y, grid_.one_over_dx_);
		wy[2] = QuadBSplineKernel(cell.y+grid_.dx_ - pos.y, grid_.one_over_dx_);

		wz[0] = QuadBSplineKernel(cell.z-grid_.dx_ - pos.z, grid_.one_over_dx_);	
		wz[1] = QuadBSplineKernel(cell.z           - pos.z, grid_.one_over_dx_);
		wz[2] = QuadBSplineKernel(cell.z+grid_.dx_ - pos.z, grid_.one_over_dx_);
	}

	void ComputeGradientStencil(const TV3& pos, const TV3& cell, TS gx[3],  TS gy[3],  TS gz[3])
	{
		gx[0] = QuadBSplineKernelGradient(cell.x-grid_.dx_ - pos.x, grid_.one_over_dx_);	
		gx[1] = QuadBSplineKernelGradient(cell.x           - pos.x, grid_.one_over_dx_);
		gx[2] = QuadBSplineKernelGradient(cell.x+grid_.dx_ - pos.x, grid_.one_over_dx_);

		gy[0] = QuadBSplineKernelGradient(cell.y-grid_.dx_ - pos.y, grid_.one_over_dx_);	
		gy[1] = QuadBSplineKernelGradient(cell.y           - pos.y, grid_.one_over_dx_);
		gy[2] = QuadBSplineKernelGradient(cell.y+grid_.dx_ - pos.y, grid_.one_over_dx_);

		gz[0] = QuadBSplineKernelGradient(cell.z-grid_.dx_ - pos.z, grid_.one_over_dx_);	
		gz[1] = QuadBSplineKernelGradient(cell.z           - pos.z, grid_.one_over_dx_);
		gz[2] = QuadBSplineKernelGradient(cell.z+grid_.dx_ - pos.z, grid_.one_over_dx_);	
	}

	void ComputeStrainTensor(const TV3& pos, TM3& tensor)
	{
		int ci, cj, ck;
		grid_.CellCenterIndex(pos, ci, cj, ck);

		TV3 cell_pos = grid_.CellCenterPosition(ci,cj,ck);
		
		ci--; cj--; ck--;

		TS wx[3],wy[3],wz[3];
		TS gx[3],gy[3],gz[3];

		ComputeWeightStencil(pos, cell_pos, wx, wy, wz);
		ComputeGradientStencil(pos, cell_pos, gx, gy, gz);

		TS dudx((TS)0), dudy((TS)0), dudz((TS)0);
		TS dvdx((TS)0), dvdy((TS)0), dvdz((TS)0);
		TS dwdx((TS)0), dwdy((TS)0), dwdz((TS)0);

		for(int k=0; k<3; k++) for(int j=0; j<3; j++) for(int i=0; i<3; i++)
		{
			int s_ix = grid_.Index3Dto1D(ci+i, cj+j, ck+k);
			TV3 cell_vel = velocity_field_[s_ix];

			dudx += gx[i] * cell_vel.x;
			dudy += gy[j] * cell_vel.x;
			dudz += gz[k] * cell_vel.x;

			dvdx += gx[i] * cell_vel.y;
			dvdy += gy[j] * cell_vel.y;
			dvdz += gz[k] * cell_vel.y;

			dwdx += gx[i] * cell_vel.z;
			dwdy += gy[j] * cell_vel.z;
			dwdz += gz[k] * cell_vel.z;		
		}

		TS D00 = dudx;
		TS D11 = dvdy;
		TS D22 = dwdz;
		TS D01 = (dudy+dvdx)*(TS)0.5;
		TS D02 = (dudz+dwdx)*(TS)0.5;
		TS D12 = (dvdz+dwdy)*(TS)0.5;

		tensor[0] = TV3(D00, D01, D02);
		tensor[1] = TV3(D01, D11, D12);
		tensor[2] = TV3(D02, D12, D22);		
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

			TV3 cell_center = grid_.CellCenterPosition(i, j, k);

			TS density = density_field_[c];

			if (density > FLT_EPSILON)
			{
				glVertex3f(cell_center.x, cell_center.y, cell_center.z);
			}
		}
		glEnd();

		glPopMatrix();
	}
};
