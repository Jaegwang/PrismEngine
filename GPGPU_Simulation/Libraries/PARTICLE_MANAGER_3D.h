
#pragma once

#include <amp.h>
#include "GRID_UNIFORM_3D.h"
#include "PARTICLE.h"
#include "VECTOR3_T.h"
#include "GL\glut.h"
#include "iostream"
#include "PARALLEL_MACRO.h"
#include "GENERAL_MACRO.h"

using namespace std;
using namespace concurrency;

class PARTICLE_MANAGER_3D
{
public:
	GRID_UNIFORM_3D grid_;

	// particle properties
	Vec3T *position_array_;
	Vec3T *velocity_array_;
	Vec3T *temp_array_;

	Vec3T **position_array_pointer_;
	Vec3T **velocity_array_pointer_;
	Vec3T **temp_array_pointer_;

	int* particle_id_array_;
	int* particle_index_array_;

	atomic<int> num_of_pts_;
	atomic<int> max_of_pts_;

	atomic<int>* num_pts_cell_;
	atomic<int>* start_idx_cell_;

public:

	PARTICLE_MANAGER_3D() : 
		position_array_(0), velocity_array_(0), temp_array_(0), position_array_pointer_(0), velocity_array_pointer_(0), temp_array_pointer_(0),
		num_pts_cell_(0), start_idx_cell_(0), particle_id_array_(0), particle_index_array_(0), max_of_pts_(0)
	{}

	~PARTICLE_MANAGER_3D()
	{
		Finalize();
	}

	void Initialize(const GRID_UNIFORM_3D& grid_input, Vec3T** pos_arr_input, Vec3T** vel_arr_input, Vec3T** temp_arr_input, const int num_pts)
	{
		grid_ = grid_input;
		max_of_pts_ = num_pts;

		position_array_ = *pos_arr_input;
		velocity_array_ = *vel_arr_input;
		temp_array_     = *temp_arr_input;

		position_array_pointer_ = pos_arr_input;
		velocity_array_pointer_ = vel_arr_input;
		temp_array_pointer_     = temp_arr_input;

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
		position_array_[ix] = Vec3T((T)FLT_MAX, (T)FLT_MAX, (T)FLT_MAX);
	}

	template<class TT>
	void RearrangeParticleData(TT** pts_arr, TT** pts_arr_temp)
	{
		BEGIN_CPU_THREADS(grid_.ijk_res_, p)
		{
			const int& num = num_pts_cell_[p];
			const int& b_ix = start_idx_cell_[p];

			for(int n = 0; n < num; n++)
			{
				int v_ix = particle_index_array_[b_ix + n];
				(*pts_arr_temp)[b_ix + n] = (*pts_arr)[v_ix];
			}
		}
		END_CPU_THREADS;	

		TT* temp = (*pts_arr_temp);
		(*pts_arr_temp) = (*pts_arr);
		(*pts_arr) = temp;
	}

	void RebuildParticleDataStructure()
	{
		if (num_of_pts_ == 0) return;

		atomic<int> start_idx = 0;
		atomic<int> count_pts = 0;

		BEGIN_CPU_THREADS(grid_.ijk_res_, p)
		{
			num_pts_cell_[p] = 0;
			start_idx_cell_[p] = -1;
		}
		END_CPU_THREADS;

		BEGIN_CPU_THREADS(num_of_pts_, p)
		{			
			const Vec3T pos = position_array_[p];
			if (grid_.IsInsideValid(pos) == false) continue;

			int i, j, k;
			grid_.CellCenterIndex(pos, i, j, k);

			int g_ix = grid_.Index3Dto1D(i, j, k);

			particle_id_array_[p] = std::atomic_fetch_add(&num_pts_cell_[g_ix], 1);

			std::atomic_fetch_add(&count_pts, 1);
		}
		END_CPU_THREADS;

		BEGIN_CPU_THREADS(grid_.ijk_res_, p)
		{
			const atomic<int>& num = num_pts_cell_[p];
			start_idx_cell_[p] = std::atomic_fetch_add(&start_idx, num);
		}
		END_CPU_THREADS;

		BEGIN_CPU_THREADS(num_of_pts_, p)
		{
			const Vec3T pos = position_array_[p];
			const int id = particle_id_array_[p];

			if (grid_.IsInsideValid(pos) == false) continue;

			int i, j, k;
			grid_.CellCenterIndex(pos, i, j, k);

			int g_ix = grid_.Index3Dto1D(i, j, k);
			int b_ix = start_idx_cell_[g_ix];

			particle_index_array_[b_ix + id] = p;
		}
		END_CPU_THREADS;

		// Rearrange Paritlce array for position and velocity
		RearrangeParticleData(position_array_pointer_, temp_array_pointer_);
		RearrangeParticleData(velocity_array_pointer_, temp_array_pointer_);

		position_array_ = *position_array_pointer_;
		velocity_array_ = *velocity_array_pointer_;

		num_of_pts_ = (int)count_pts;
	}

	void Rendering()
	{
		glDisable(GL_LIGHTING);

		glPushMatrix();
		glColor3f(0.3, 0.3, 1);

		glPointSize(1.2);

		T dt = (T)0.01;
		T dist = MAX(grid_.dx_*(T)3, grid_.dy_*(T)3, grid_.dz_*(T)3);
		T max_vel = dist / dt;

		glBegin(GL_POINTS);

		for (int i = 0; i < num_of_pts_; i++)
		{
			const Vec3T& pos = position_array_[i];
			const T v = velocity_array_[i].Magnitude();
			const T g = v / max_vel;

			glColor3f(g, g, (T)1);

			glVertex3f(pos.x, pos.y, pos.z);
		}
		glEnd();

		glPopMatrix();
	}
};