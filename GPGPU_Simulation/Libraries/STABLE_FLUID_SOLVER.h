
#pragma once 

#include "FIELD.h"
#include "FIELD_UNIFORM.h"
#include "ADVECTION_METHOD.h"
#include "PROJECTION_METHOD.h"

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

	ADVECTION_METHOD  advection_method_;
	PROJECTION_METHOD projection_method_;

public:

	STABLE_FLUID_SOLVER() : density_field_(0), velocity_field_(0)
	{}
	~STABLE_FLUID_SOLVER()
	{}

	void Initialize(const GRID& grid)
	{
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
				boundary_field_->Set(i,j,k, BND_WALL);
			else
				boundary_field_->Set(i,j,k, BND_FULL);
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
			if(density_field_->Get(i,j,k) > 0.1)
			{
				Vec3 cell_center = grid.CellCenterPosition(i,j,k);
				glVertex3fv(&cell_center.x);			
			}
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
				Vec3 v = cell_center + (velocity * 0.1f);

				glVertex3fv(&cell_center.x);			
				glVertex3fv(&v.x);
			}
		}	

		glEnd();
		glEnable(GL_LIGHTING);
	}
};