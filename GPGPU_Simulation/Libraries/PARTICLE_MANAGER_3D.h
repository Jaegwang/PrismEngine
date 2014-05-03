
#pragma once

#include <amp.h>
#include "GRID_UNIFORM_3D.h"
#include "PARTICLE.h"
#include "GL\glut.h"
#include "iostream"
#include "GENERAL_MACRO.h"
#include "MATH_CORE.h"

using namespace std;
using namespace concurrency;

class PARTICLE_MANAGER_3D
{
public:
	GRID_UNIFORM_3D grid_;

	// particle properties
	Vec3 **position_array_pointer_;
	Vec3 **velocity_array_pointer_;

	int* particle_id_array_;
	int* particle_index_array_;

	atomic<int> num_of_pts_;
	atomic<int> max_of_pts_;

	atomic<int>* num_pts_cell_;
	atomic<int>* start_idx_cell_;

public:

	PARTICLE_MANAGER_3D() : 
		position_array_pointer_(0), velocity_array_pointer_(0),
		num_pts_cell_(0), start_idx_cell_(0), particle_id_array_(0), particle_index_array_(0), max_of_pts_(0)
	{}

	~PARTICLE_MANAGER_3D()
	{
		Finalize();
	}

	void Initialize(const GRID_UNIFORM_3D& grid_input, Vec3** pos_arr_input, Vec3** vel_arr_input, const int num_pts)
	{
		grid_ = grid_input;
		max_of_pts_ = num_pts;

		position_array_pointer_ = pos_arr_input;
		velocity_array_pointer_ = vel_arr_input;

		particle_id_array_ = new int[num_pts];
		particle_index_array_ = new int[num_pts];

		num_pts_cell_ = new atomic<int>[grid_.ijk_res_];
		start_idx_cell_ = new atomic<int>[grid_.ijk_res_];
	}

	void Finalize()
	{
		SAFE_DELETE_ARRAY(particle_id_array_);
		SAFE_DELETE_ARRAY(particle_index_array_);

		SAFE_DELETE_ARRAY(num_pts_cell_);
		SAFE_DELETE_ARRAY(start_idx_cell_);
	}

	int AddParticle()
	{
		return (int)std::atomic_fetch_add(&num_of_pts_, 1);
	}

	void DelParticle(int ix)
	{
		(*position_array_pointer_)[ix] = Vec3((FLT)FLT_MAX, (FLT)FLT_MAX, (FLT)FLT_MAX);
	}

	template<class TT>
	void RearrangeParticleData(TT** pts_arr, TT** pts_arr_temp)
	{
		#pragma omp parallel for
		for(int p=0; p<grid_.ijk_res_; p++)
		{
			const int& num = num_pts_cell_[p];
			const int& b_ix = start_idx_cell_[p];

			for(int n = 0; n < num; n++)
			{
				int v_ix = particle_index_array_[b_ix + n];
				(*pts_arr_temp)[b_ix + n] = (*pts_arr)[v_ix];
			}
		}

		TT* temp = (*pts_arr_temp);
		(*pts_arr_temp) = (*pts_arr);
		(*pts_arr) = temp;
	}

	void RebuildParticleDataStructure()
	{
		if (num_of_pts_ == 0) return;

		Vec3 *pos_arr = *position_array_pointer_;

		atomic<int> start_idx = 0;
		atomic<int> count_pts = 0;

		#pragma omp for 
		for (int p = 0 ; p < grid_.ijk_res_ ; p++) 
		{
			num_pts_cell_[p] = 0;
			start_idx_cell_[p] = -1;
		}

		#pragma omp for 
		for (int p = 0 ; p < num_of_pts_; p++) 
		{			
			const Vec3 pos = pos_arr[p];
			if (grid_.IsInsideValid(pos) == false) continue;

			int i, j, k;
			grid_.CellCenterIndex(pos, i, j, k);

			int g_ix = grid_.Index3Dto1D(i, j, k);

			particle_id_array_[p] = std::atomic_fetch_add(&num_pts_cell_[g_ix], 1);

			std::atomic_fetch_add(&count_pts, 1);
		}

		#pragma omp for 
		for (int p = 0 ; p < grid_.ijk_res_; p++) 
		{
			const atomic<int>& num = num_pts_cell_[p];
			start_idx_cell_[p] = std::atomic_fetch_add(&start_idx, num);
		}

		#pragma omp for 
		for (int p = 0 ; p < num_of_pts_; p++) 
		{
			const Vec3 pos = pos_arr[p];
			const int id = particle_id_array_[p];

			if (grid_.IsInsideValid(pos) == false) continue;

			int i, j, k;
			grid_.CellCenterIndex(pos, i, j, k);

			int g_ix = grid_.Index3Dto1D(i, j, k);
			int b_ix = start_idx_cell_[g_ix];

			particle_index_array_[b_ix + id] = p;
		}

		num_of_pts_ = (int)count_pts;
	}

	void Rendering()
	{
		glDisable(GL_LIGHTING);

		glPushMatrix();
		glColor3f(0.3, 0.3, 1);

		glPointSize(1.3);

		FLT dt = (FLT)0.01;
		FLT dist = MAX3(grid_.dx_*(FLT)3, grid_.dx_*(FLT)3, grid_.dx_*(FLT)3);
		FLT max_vel = dist / dt;

		glBegin(GL_POINTS);

		for (int i = 0; i < num_of_pts_; i++)
		{
			const Vec3& pos =(*position_array_pointer_)[i];
			const FLT v = glm::length((*velocity_array_pointer_)[i]);
			const FLT g = v / max_vel;

			glColor3f(g, g, (FLT)1);

			glVertex3f(pos.x, pos.y, pos.z);
		}
		glEnd();

		glPopMatrix();
	}
};