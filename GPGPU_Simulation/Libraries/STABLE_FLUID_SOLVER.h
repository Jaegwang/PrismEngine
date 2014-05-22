
#pragma once 

#include "FIELD.h"
#include "FIELD_UNIFORM.h"
#include "FIELD_DYNAMIC.h"
#include "ADVECTION_METHOD.h"
#include "RASTERIZER.h"
#include "PROJECTION_METHOD.h"
#include "ARRAY.h"

class STABLE_FLUID_SOLVER
{
public:

	FIELD<FLT>*  density_field_;
	FIELD<Vec3>* velocity_field_;

	FIELD<FLT>* divergence_field_;
	FIELD<FLT>* pressure_field_;
	FIELD<int>* boundary_field_;

	FIELD<FLT>*  scalar_ghost_field_;
	FIELD<Vec3>* vector_ghost_field_;

	// particle data
	ARRAY<Vec3>* particle_position_;
	ARRAY<Vec3>* particle_velocity_;

	ARRAY<Vec3>* particle_vector_ghost_;

	// method
	ADVECTION_METHOD  advection_method_;
	PROJECTION_METHOD projection_method_;

public:

	STABLE_FLUID_SOLVER() : density_field_(0), velocity_field_(0)
	{}
	~STABLE_FLUID_SOLVER()
	{}

	void Initialize(const GRID& grid)
	{
		ARRAY_VECTOR<Vec3>* p_p = new ARRAY_VECTOR<Vec3>;
		ARRAY_VECTOR<Vec3>* p_v = new ARRAY_VECTOR<Vec3>;

		p_p->Initialize(1000000);
		p_v->Initialize(1000000);

		particle_position_ = p_p;
		particle_velocity_ = p_v;			


		FIELD_UNIFORM<FLT>*  density_uniform = new FIELD_UNIFORM<FLT>;
		FIELD_UNIFORM<Vec3>* velocity_uniform = new FIELD_UNIFORM<Vec3>;

		FIELD_UNIFORM<FLT>*  divergence_uniform_ = new FIELD_UNIFORM<FLT>;
		FIELD_UNIFORM<FLT>*  pressure_uniform_ =  new FIELD_UNIFORM<FLT>;
		FIELD_UNIFORM<int>*  boundary_uniform_ = new FIELD_UNIFORM<int>;

		FIELD_UNIFORM<FLT>*  scalar_ghost_uniform_ = new FIELD_UNIFORM<FLT>;
		FIELD_UNIFORM<Vec3>* vector_ghost_uniform_ = new FIELD_UNIFORM<Vec3>;

		density_uniform->Initialize(grid);
		velocity_uniform->Initialize(grid);

		divergence_uniform_->Initialize(grid);
		pressure_uniform_->Initialize(grid);
		boundary_uniform_->Initialize(grid);

		scalar_ghost_uniform_->Initialize(grid);
		vector_ghost_uniform_->Initialize(grid);

		density_field_ = density_uniform;
		velocity_field_ = velocity_uniform;

		divergence_field_ = divergence_uniform_;
		pressure_field_ = pressure_uniform_;
		boundary_field_ = boundary_uniform_;

		scalar_ghost_field_ = scalar_ghost_uniform_;
		vector_ghost_field_ = vector_ghost_uniform_;

		#pragma omp parallel for
		for(int k=0; k<grid.k_res_; k++)
		for(int j=0; j<grid.j_res_; j++)
		for(int i=0; i<grid.i_res_; i++)
		{
			density_field_->Set(i,j,k,(FLT)0);
			velocity_field_->Set(i,j,k,Vec3(0,0,0));
		}
	}

	void SetBoundaryCondition()
	{
		GRID grid = boundary_field_->Grid();

		#pragma omp parallel for
		for(int k=0; k<grid.k_res_; k++)
		for(int j=0; j<grid.j_res_; j++)
		for(int i=0; i<grid.i_res_; i++)
		{
			if(grid.IsGhostCell(i,j,k))
			{
				boundary_field_->Set(i,j,k, BND_WALL);
				continue;
			}

			if(density_field_->Get(i,j,k) > (FLT)0)
				boundary_field_->Set(i,j,k, BND_FULL);
			else
				boundary_field_->Set(i,j,k, BND_NULL);
		}		
	
	}

	void SourceBoundaryFromSphere(const Vec3 pos, const FLT rad)
	{
		GRID grid = density_field_->Grid();

		#pragma omp parallel for
		for(int k=0; k<grid.k_res_; k++)
		for(int j=0; j<grid.j_res_; j++)
		for(int i=0; i<grid.i_res_; i++)
		{
			Vec3 cell_center = grid.CellCenterPosition(i,j,k);
			FLT  dist = glm::length(cell_center-pos);

			if(dist < rad)
			{
				boundary_field_->Set(i,j,k, BND_WALL);
			}
		}	
	}

	void SourceDensityFromSphere(const Vec3 pos, const FLT rad, const FLT den, const Vec3 vel)
	{
		GRID grid = density_field_->Grid();

		#pragma omp parallel for
		for(int k=0; k<grid.k_res_; k++)
		for(int j=0; j<grid.j_res_; j++)
		for(int i=0; i<grid.i_res_; i++)
		{
			Vec3 cell_center = grid.CellCenterPosition(i,j,k);
			FLT  dist = glm::length(cell_center-pos);

			if(dist < rad)
			{
				density_field_->Set(i,j,k,den);
				velocity_field_->Set(i,j,k,vel);
			}
		}	
	}
	
	void SeedParticlesFromSphere(const Vec3 pos, const FLT rad, const Vec3 vel, const int num)
	{
		for(int n=0; n<num; n++)
		{
			Vec3 new_pos = pos + Vec3((FLT)rand()/(FLT)RAND_MAX-0.5, (FLT)rand()/(FLT)RAND_MAX-0.5, (FLT)rand()/(FLT)RAND_MAX-0.5)*2.0f*rad;
			Vec3 new_vel = vel;

			if(glm::length(new_pos-pos) <= rad)
			{
				particle_position_->Push(new_pos);
				particle_velocity_->Push(new_vel);
			}		
		}
	}

	void AdvanceOneTimeStepFLIP(const FLT dt)
	{
		// Rasterization
		{
			RasterizeParticleToField(*velocity_field_, *density_field_, *particle_position_, *particle_velocity_);		
		}

		// Projection
		{
			GRID grid = velocity_field_->Grid();

			FLT flip_coeff = (FLT)0.95;

			#pragma omp parallel for
			for(int i=0; i<grid.ijk_res_; i++)
			{
				vector_ghost_field_->Set(i, velocity_field_->Get(i));			
			}

			SetBoundaryCondition();
			projection_method_.Jacobi(boundary_field_, velocity_field_, divergence_field_, pressure_field_, scalar_ghost_field_, 200);

			#pragma omp parallel for
			for(int i=0; i<grid.ijk_res_; i++)
			{
				vector_ghost_field_->Set(i, velocity_field_->Get(i) - vector_ghost_field_->Get(i));
			}

			#pragma omp parallel for
			for(int i=0; i<particle_position_->Size(); i++)
			{
				Vec3 pos = particle_position_->Get(i);
				Vec3 vel = particle_velocity_->Get(i);

				Vec3 vel_g = velocity_field_->Get(pos);
				Vec3 vel_p = vel + vector_ghost_field_->Get(pos);

				particle_velocity_->Set(i, vel_g*((FLT)1-flip_coeff) + vel_p*flip_coeff);
			}
		}		

		// Advection
		{
			GRID grid = velocity_field_->Grid();

			#pragma omp parallel for
			for(int i=0; i<particle_position_->Size(); i++)
			{
				Vec3 pos = particle_position_->Get(i);
				Vec3 vel = particle_velocity_->Get(i);

				pos += vel*dt;

				particle_position_->Set(i, grid.Clamp(pos));
			}
		}

		// forcing
		{
			#pragma omp parallel for
			for(int i=0; i<particle_position_->Size(); i++)
			{
				Vec3 gravity = Vec3(0,-10,0);
				Vec3 vel = particle_velocity_->Get(i) + gravity * dt;
				
				particle_velocity_->Set(i, vel);
			}
		}
	}

	void AdvanceOneTimeStep(const FLT dt)
	{
		// Projection
		{
			SetBoundaryCondition();
			projection_method_.Jacobi(boundary_field_, velocity_field_, divergence_field_, pressure_field_, scalar_ghost_field_, 20);
		}

		// Advection
		{
			advection_method_.SemiLagrangian(velocity_field_, dt, density_field_, scalar_ghost_field_);
			advection_method_.SemiLagrangian(velocity_field_, dt, velocity_field_, vector_ghost_field_);
			FIELD<FLT>*  scalar_temp;
			FIELD<Vec3>* vector_temp;

			SWAP(density_field_, scalar_ghost_field_, scalar_temp);
			SWAP(velocity_field_, vector_ghost_field_, vector_temp);
		}

		// Rebulid
		{
			density_field_->Rebuild();
			velocity_field_->Rebuild();

			divergence_field_->Rebuild();
			pressure_field_->Rebuild();
			boundary_field_->Rebuild();

			scalar_ghost_field_->Rebuild();
			vector_ghost_field_->Rebuild();
		}
	}

	void Render()
	{
		GRID grid = density_field_->Grid();

		glDisable(GL_LIGHTING);
		glBegin(GL_POINTS);

		for(int k=0; k<grid.k_res_; k++)
		for(int j=0; j<grid.j_res_; j++)
		for(int i=0; i<grid.i_res_; i++)
		{
			if(density_field_->Get(i,j,k) == 0.0) continue;

			if(ABS(divergence_field_->Get(i,j,k)) < 0.05)
			{
				Vec3 cell_center = grid.CellCenterPosition(i,j,k);
				glColor3f(1,0,0);
				glVertex3fv(&cell_center.x);			
			}
		}	

		glEnd();
		glEnable(GL_LIGHTING);
	}

	void RenderParticles()
	{
		glDisable(GL_LIGHTING);

		const int size = particle_position_->Size();

		glBegin(GL_POINTS);
		for(int i=0; i<size; i++)
		{
			Vec3 pos = particle_position_->Get(i);
			glColor3f(0,0,0);
			glVertex3fv(&pos.x);
		}

		glEnd();
		glEnable(GL_LIGHTING);
	}

	void RenderVelocity()
	{
		GRID grid = density_field_->Grid();

		glDisable(GL_LIGHTING);
		glBegin(GL_LINES);

		glColor3f(1,0,0);
		for(int k=0; k<grid.k_res_; k++)
		for(int j=0; j<grid.j_res_; j++)
		for(int i=0; i<grid.i_res_; i++)
		{
			Vec3 velocity = velocity_field_->Get(i,j,k);
			Vec3 cell_center = grid.CellCenterPosition(i,j,k);

			if(glm::length(velocity) > 0)
			{
				Vec3 v = cell_center + (velocity * 0.01f);

				glVertex3fv(&cell_center.x);			
				glVertex3fv(&v.x);
			}
		}	

		glEnd();
		glEnable(GL_LIGHTING);
	}
};