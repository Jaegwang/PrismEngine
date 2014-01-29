
#pragma once

#include <amp.h>
#include "GRID_UNIFORM_3D.h"
#include "PARTICLE.h"
#include "VECTOR3_T.h"
#include "GL\glut.h"
#include "iostream"
#include "PARALLEL_MACRO.h"

using namespace std;
using namespace concurrency;

class PARTICLE_MANAGER_3D
{
public:

	GRID_UNIFORM_3D grid_;

	// particle properties
	Vec3T *position_array_, *velocity_array_, *force_array_;
	Vec3T *position_array_old_, *velocity_array_old_, *force_array_old_;

	T *density_array_;
	T *density_array_old_;

	int* particle_id_array_;
	int* particle_index_array_;

	atomic<int> num_of_pts_;
	atomic<int> max_of_pts_;

	atomic<int>* num_pts_cell_;
	atomic<int>* start_idx_cell_;

public:

	PARTICLE_MANAGER_3D() : 
		position_array_(0), position_array_old_(0), velocity_array_(0), velocity_array_old_(0), density_array_(0), density_array_old_(0),
		force_array_(0), force_array_old_(0),
		num_pts_cell_(0), start_idx_cell_(0), particle_id_array_(0), particle_index_array_(0), max_of_pts_(0)
	{}

	~PARTICLE_MANAGER_3D()
	{
		if (position_array_) delete[] position_array_;
		if (position_array_old_) delete[] position_array_old_;

		if (velocity_array_) delete[] velocity_array_;
		if (velocity_array_old_) delete[] velocity_array_old_;

		if (force_array_) delete[] force_array_;
		if (force_array_old_) delete[] force_array_old_;

		if (density_array_) delete[] density_array_;
		if (density_array_old_) delete[] density_array_old_;

		if(particle_id_array_) delete[] particle_id_array_;
		if(particle_index_array_) delete[] particle_index_array_;

		if(num_pts_cell_) delete[] num_pts_cell_;
		if(start_idx_cell_) delete[] start_idx_cell_;
	}

	void Initialize(const GRID_UNIFORM_3D& grid_input, const int num_pts)
	{
		grid_ = grid_input;
		max_of_pts_ = num_pts;

		position_array_ = new Vec3T[num_pts];
		position_array_old_ = new Vec3T[num_pts];

		velocity_array_ = new Vec3T[num_pts];
		velocity_array_old_ = new Vec3T[num_pts];

		force_array_ = new Vec3T[num_pts];
		force_array_old_ = new Vec3T[num_pts];

		density_array_ = new T[num_pts];
		density_array_old_ = new T[num_pts];


		particle_id_array_ = new int[num_pts];
		particle_index_array_ = new int[num_pts];

		num_pts_cell_ = new atomic<int>[grid_.ijk_res_];
		start_idx_cell_ = new atomic<int>[grid_.ijk_res_];

		memset((void*)density_array_, 0, sizeof(T)*max_of_pts_);
		memset((void*)density_array_old_, 0, sizeof(T)*max_of_pts_);
	}

	int AddParticle()
	{
		return (int)std::atomic_fetch_add(&num_of_pts_, 1);
	}

	void DelParticle(int ix)
	{
		T max_pos = FLT_MAX;
		position_array_[ix] = Vec3T(max_pos, max_pos, max_pos);
	}

	void RebuildParticleDataStructure()
	{
		if (num_of_pts_ == 0) return;

		atomic<int> start_idx = 0;
		atomic<int> count_pts = 0;

		BEGIN_CPU_THREADS_1D(grid_.ijk_res_)
		{
			for (int p = ix_begin; p <= ix_end; p++)
			{
				num_pts_cell_[p] = 0;
				start_idx_cell_[p] = -1;
			}
		}
		END_CPU_THREADS_1D;

		BEGIN_CPU_THREADS_1D(num_of_pts_)
		{
			for (int p = ix_begin; p <= ix_end; p++)
			{
				const Vec3T pos = position_array_[p];

				if (grid_.IsInsideValid(pos) == false) continue;

				int i, j, k;
				grid_.CellCenterIndex(pos, i, j, k);

				int g_ix = grid_.Index3Dto1D(i, j, k);

				particle_id_array_[p] = std::atomic_fetch_add(&num_pts_cell_[g_ix], 1);

				std::atomic_fetch_add(&count_pts, 1);
			}
		}
		END_CPU_THREADS_1D;

		BEGIN_CPU_THREADS_1D(grid_.ijk_res_)
		{
			for (int p = ix_begin; p <= ix_end; p++)
			{
				const atomic<int>& num = num_pts_cell_[p];
				start_idx_cell_[p] = std::atomic_fetch_add(&start_idx, num);
			}
		}
		END_CPU_THREADS_1D;

		BEGIN_CPU_THREADS_1D(num_of_pts_)
		{
			for (int p = ix_begin; p <= ix_end; p++)
			{
				const Vec3T pos = position_array_[p];
				const int pts_id = particle_id_array_[p];

				if (grid_.IsInsideValid(pos) == false) continue;

				int i, j, k;
				grid_.CellCenterIndex(pos, i, j, k);

				int g_ix = grid_.Index3Dto1D(i, j, k);

				int b_ix = start_idx_cell_[g_ix];

				particle_index_array_[b_ix + pts_id] = p;
			}
		}
		END_CPU_THREADS_1D;

		BEGIN_CPU_THREADS_1D(grid_.ijk_res_)
		{
			for (int p = ix_begin; p <= ix_end; p++)
			{
				const int& num = num_pts_cell_[p];
				const int& b_ix = start_idx_cell_[p];

				for (int n = 0; n < num; n++)
				{
					int v_ix = particle_index_array_[b_ix + n];

					force_array_old_[b_ix + n] = force_array_[v_ix];
					position_array_old_[b_ix + n] = position_array_[v_ix];
					velocity_array_old_[b_ix + n] = velocity_array_[v_ix];
					density_array_old_[b_ix + n]  = density_array_[v_ix];
				}
			}
		}
		END_CPU_THREADS_1D;

		Vec3T* v_temp = 0; T* t_temp = 0;

		SWAP(position_array_, position_array_old_, v_temp);
		SWAP(velocity_array_, velocity_array_old_, v_temp);
		SWAP(force_array_, force_array_old_, v_temp);
		SWAP(density_array_, density_array_old_, t_temp);

		num_of_pts_ = (int)count_pts;
	}

	void Rendering()
	{
		glDisable(GL_LIGHTING);

		glPushMatrix();
		glColor3f(0.3, 0.3, 1);

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