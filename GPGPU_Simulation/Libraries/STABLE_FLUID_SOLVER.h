
#pragma once 

#include "FIELD.h"
#include "FIELD_UNIFORM.h"
#include "FIELD_ENCODED.h"
#include "ADVECTION_METHOD.h"
#include "RASTERIZER.h"
#include "PROJECTION_METHOD.h"
#include "ARRAY.h"

class STABLE_FLUID_SOLVER
{
public:

	FIELD<TS>*  density_field_;
	FIELD<TV3>* velocity_field_;

	FIELD<TS>*  divergence_field_;
	FIELD<TS>*  pressure_field_;
	FIELD<int>* boundary_field_;

	FIELD<TS>*  scalar_ghost_field_;
	FIELD<TV3>* vector_ghost_field_;

	// particle data
	ARRAY<TV3>*  particle_position_;
	ARRAY<TV3>*  particle_velocity_;
	ARRAY<bool>* particle_active_;

	ARRAY<TV3>*  particle_vector_ghost_;

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
		ARRAY_VECTOR<TV3>* p_p = new ARRAY_VECTOR<TV3>;
		ARRAY_VECTOR<TV3>* p_v = new ARRAY_VECTOR<TV3>;
		ARRAY_VECTOR<bool>* p_a = new ARRAY_VECTOR<bool>;

		p_p->Initialize(2000000);
		p_v->Initialize(2000000);
		p_a->Initialize(2000000);

		particle_position_ = p_p;
		particle_velocity_ = p_v;			
		particle_active_   = p_a;

		FIELD_ENCODED<TS>*  density_uniform = new FIELD_ENCODED<TS>;
		FIELD_ENCODED<TV3>* velocity_uniform = new FIELD_ENCODED<TV3>;

		FIELD_ENCODED<TS>*  divergence_uniform_ = new FIELD_ENCODED<TS>;
		FIELD_ENCODED<TS>*  pressure_uniform_ =  new FIELD_ENCODED<TS>;
		FIELD_ENCODED<int>*  boundary_uniform_ = new FIELD_ENCODED<int>;

		FIELD_ENCODED<TS>*  scalar_ghost_uniform_ = new FIELD_ENCODED<TS>;
		FIELD_ENCODED<TV3>* vector_ghost_uniform_ = new FIELD_ENCODED<TV3>;

		density_uniform->Initialize(grid, (TS)0);
		velocity_uniform->Initialize(grid, TV3());

		divergence_uniform_->Initialize(grid, (TS)0);
		pressure_uniform_->Initialize(grid, (TS)0);
		boundary_uniform_->Initialize(grid, BND_NULL);

		scalar_ghost_uniform_->Initialize(grid, (TS)0);
		vector_ghost_uniform_->Initialize(grid, TV3());

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
			density_field_->Set(i,j,k,(TS)0);
			velocity_field_->Set(i,j,k,TV3(0,0,0));
		}
	}

	GRID Grid()
	{
		return velocity_field_->Grid();
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

			if(density_field_->Get(i,j,k) > (TS)0.0)
				boundary_field_->Set(i,j,k, BND_FULL);
			else
				boundary_field_->Set(i,j,k, BND_NULL);
		}		
	
	}

	void SourceBoundaryFromSphere(const TV3 pos, const TS rad)
	{
		GRID grid = density_field_->Grid();

		#pragma omp parallel for
		for(int k=0; k<grid.k_res_; k++)
		for(int j=0; j<grid.j_res_; j++)
		for(int i=0; i<grid.i_res_; i++)
		{
			TV3 cell_center = grid.CellCenterPosition(i,j,k);
			TS  dist = glm::length(cell_center-pos);

			if(dist < rad)
			{
				boundary_field_->Set(i,j,k, BND_WALL);
			}
		}	
	}

	void SourceDensityFromSphere(const TV3 pos, const TS rad, const TS den, const TV3 vel)
	{
		GRID grid = density_field_->Grid();

		#pragma omp parallel for
		for(int k=0; k<grid.k_res_; k++)
		for(int j=0; j<grid.j_res_; j++)
		for(int i=0; i<grid.i_res_; i++)
		{
			TV3 cell_center = grid.CellCenterPosition(i,j,k);
			TS  dist = glm::length(cell_center-pos);

			if(dist < rad)
			{
				density_field_->Set(i,j,k,den);
				velocity_field_->Set(i,j,k,vel);
			}
		}	
	}
	
	void SeedParticlesFromSphere(const TV3 pos, const TS rad, const TV3 vel, const int num)
	{
		for(int n=0; n<num; n++)
		{
			TV3 new_pos = pos + TV3((TS)rand()/(TS)RAND_MAX-0.5, (TS)rand()/(TS)RAND_MAX-0.5, (TS)rand()/(TS)RAND_MAX-0.5)*2.0f*rad;
			TV3 new_vel = vel;

			if(glm::length(new_pos-pos) <= rad)
			{
				particle_position_->Push(new_pos);
				particle_velocity_->Push(new_vel);
				particle_active_->Push(true);
			}		
		}
	}

	void SeedParticlesFromHeight(const TS height, const int num)
	{
		GRID grid = velocity_field_->Grid();

		for(int i=0; i<grid.i_res_; i++)
		for(int j=0; j<grid.j_res_; j++)
		for(int k=0; k<grid.k_res_; k++)
		{
			if(grid.IsGhostCell(i,j,k) == true) continue;

			TV3 cell_center = grid.CellCenterPosition(i,j,k);

			TS rad = grid.dx_*(TS)0.5;

			if(cell_center.y < height)
			{
				for(int n=0; n<num; n++)
				{
					TV3 new_pos = cell_center + TV3((TS)rand()/(TS)RAND_MAX-0.5, (TS)rand()/(TS)RAND_MAX-0.5, (TS)rand()/(TS)RAND_MAX-0.5)*2.0f*rad;
					TV3 new_vel = TV3();

					particle_position_->Push(new_pos);
					particle_velocity_->Push(new_vel);
					particle_active_->Push(false);
				}
			}		
		}	
	}

	void AdvanceOneTimeStepFLIP(const TS dt)
	{
		// Advection
		{
			GRID grid = velocity_field_->Grid();

			FOR_EACH_PARALLER(i, 0, particle_position_->Size()-1)		
			{
				TV3 pos = particle_position_->Get(i);
				TV3 vel_g = velocity_field_->Get(pos);

				if(glm::length(vel_g) != (TS)0)
				{
					pos += vel_g*dt;
					particle_active_->Set(i, true);
				}

				particle_position_->Set(i, grid.Clamp(pos));				
			}
		}

		// Rasterization
		{
			RasterizeParticleToField(*velocity_field_, *density_field_, *particle_position_, *particle_velocity_);		
		}

		// Projection
		{
			GRID grid = velocity_field_->Grid();

			TS flip_coeff = (TS)0.95;

			FOR_EACH_PARALLER(i, 0, grid.ijk_res_-1)
			{
				vector_ghost_field_->Set(i, velocity_field_->Get(i));			
			}

			SetBoundaryCondition();
			projection_method_.Jacobi(boundary_field_, velocity_field_, divergence_field_, pressure_field_, scalar_ghost_field_, dt, 200);

			FOR_EACH_PARALLER(i, 0, grid.ijk_res_-1)
			{
				vector_ghost_field_->Set(i, velocity_field_->Get(i) - vector_ghost_field_->Get(i));
			}

			FOR_EACH_PARALLER(i, 0, particle_position_->Size()-1)
			{
				TV3 pos = particle_position_->Get(i);
				TV3 vel = particle_velocity_->Get(i);

				TV3 vel_g = velocity_field_->Get(pos);
				TV3 vel_p = vel + vector_ghost_field_->Get(pos);

				particle_velocity_->Set(i, vel_g*((TS)1-flip_coeff) + vel_p*flip_coeff);
			}
		}		

		// forcing
		{
			FOR_EACH_PARALLER(i, 0, particle_position_->Size()-1)
			{
				if(particle_active_->Get(i) == true)
				{
					TV3 gravity = TV3(0,-10,0);
					TV3 vel = particle_velocity_->Get(i) + gravity * dt;
				
					particle_velocity_->Set(i, vel);
				}
			}
		}
	}

	void AdvanceOneTimeStep(const TS dt)
	{
		// Projection
		{
			SetBoundaryCondition();
			projection_method_.Jacobi(boundary_field_, velocity_field_, divergence_field_, pressure_field_, scalar_ghost_field_, dt, 20);
		}

		// Advection
		{
			advection_method_.SemiLagrangian(velocity_field_, dt, density_field_, scalar_ghost_field_);
			advection_method_.SemiLagrangian(velocity_field_, dt, velocity_field_, vector_ghost_field_);
			FIELD<TS>*  scalar_temp;
			FIELD<TV3>* vector_temp;

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
				TV3 cell_center = grid.CellCenterPosition(i,j,k);
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

		const TS dt = 0.01;
		const TS cfl = velocity_field_->Grid().dx_*(TS)2;

		glBegin(GL_POINTS);
		for(int i=0; i<size; i++)
		{
			TV3 pos = particle_position_->Get(i);
			TV3 vel = particle_velocity_->Get(i);

			TS mag_vel = glm::length(vel);

			TS c = mag_vel*dt/cfl;

			glColor3f(c,c,1);

			glVertex3fv(&pos.x);
		}

		glEnd();
		glEnable(GL_LIGHTING);
	}

	void RenderVelocity()
	{
		GRID grid = density_field_->Grid();

		glDisable(GL_LIGHTING);
		glBegin(GL_POINTS);

		glColor3f(1,0,0);
		for(int k=0; k<grid.k_res_; k++)
		for(int j=0; j<grid.j_res_; j++)
		for(int i=0; i<grid.i_res_; i++)
		{
			TV3 velocity = velocity_field_->Get(i,j,k);

			TV3 cell_center = grid.CellCenterPosition(i,j,k);

			TS den = density_field_->Get(i,j,k);

			if(den > (TS)0)
			{
//				glVertex3fv(&cell_center.x);							
			}

			//if(glm::length(velocity) > 0)
			//{
			//	TV3 v = cell_center + (velocity * 0.01f);

			//	glVertex3fv(&cell_center.x);			
			//	glVertex3fv(&v.x);
			//}
		}	

		glEnd();
		glEnable(GL_LIGHTING);
	}
};