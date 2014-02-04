
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
	Vec3T *position_array_, *velocity_array_;

	std::vector<Vec3T**> vector_data_pointer_;
	std::vector<T**>     scalar_data_pointer_;

	std::vector<Vec3T*>  vector_data_;
	std::vector<T*>      scalar_data_;

	std::vector<Vec3T*>  vector_data_temp_;
	std::vector<T*>      scalar_data_temp_;

	int* particle_id_array_;
	int* particle_index_array_;

	atomic<int> num_of_pts_;
	atomic<int> max_of_pts_;

	atomic<int>* num_pts_cell_;
	atomic<int>* start_idx_cell_;

	// array for cell including particles
	int* pts_cell_index_array_;
	atomic<int> num_of_pts_cell_;

public:

	PARTICLE_MANAGER_3D() : 
		position_array_(0), velocity_array_(0), num_pts_cell_(0), start_idx_cell_(0), particle_id_array_(0), particle_index_array_(0), max_of_pts_(0), num_of_pts_cell_(0)
	{}

	~PARTICLE_MANAGER_3D()
	{
		Finalize();
	}

	void Initialize(const GRID_UNIFORM_3D& grid_input, Vec3T** pos_arr_input, Vec3T** vel_arr_input, const int num_pts)
	{
		grid_ = grid_input;
		max_of_pts_ = num_pts;

		position_array_ = *pos_arr_input;
		vector_data_pointer_.push_back(pos_arr_input);
		vector_data_.push_back(*pos_arr_input);
		vector_data_temp_.push_back(new Vec3T[num_pts]);		

		velocity_array_ = *vel_arr_input;
		vector_data_pointer_.push_back(vel_arr_input);
		vector_data_.push_back(*vel_arr_input);
		vector_data_temp_.push_back(new Vec3T[num_pts]);

		particle_id_array_ = new int[num_pts];
		particle_index_array_ = new int[num_pts];

		num_pts_cell_ = new atomic<int>[grid_.ijk_res_];
		start_idx_cell_ = new atomic<int>[grid_.ijk_res_];

		pts_cell_index_array_ = new int[grid_.ijk_res_];
	}

	void Finalize()
	{
		SAFE_DELETE_ARRAY(particle_id_array_);
		SAFE_DELETE_ARRAY(particle_index_array_);
		SAFE_DELETE_ARRAY(num_pts_cell_);
		SAFE_DELETE_ARRAY(start_idx_cell_);
		SAFE_DELETE_ARRAY(pts_cell_index_array_);

		for(int i = 0; i < (int)vector_data_temp_.size(); i++)
			SAFE_DELETE_ARRAY(vector_data_temp_[i]);

		for(int i = 0; i < (int)scalar_data_temp_.size(); i++)
			SAFE_DELETE_ARRAY(scalar_data_temp_[i]);
	}

	void AddVectorData(Vec3T** data_arr_input, const int num_pts)
	{
		vector_data_pointer_.push_back(data_arr_input);
		vector_data_.push_back(*data_arr_input);
		vector_data_temp_.push_back(new Vec3T[num_pts]);
	}

	void AddScalarData(T** data_arr_input, const int num_pts)
	{
		scalar_data_pointer_.push_back(data_arr_input);
		scalar_data_.push_back(*data_arr_input);
		scalar_data_temp_.push_back(new T[num_pts]);
	}

	int AddParticle()
	{
		return (int)std::atomic_fetch_add(&num_of_pts_, 1);
	}

	void DelParticle(int ix)
	{
		position_array_[ix] = Vec3T((T)FLT_MAX, (T)FLT_MAX, (T)FLT_MAX);
	}

	void RebuildParticleDataStructure()
	{
		if (num_of_pts_ == 0) return;

		atomic<int> start_idx = 0;
		atomic<int> count_pts = 0;
		num_of_pts_cell_ = 0;

		array_view<int, 1> num_pts_cell_view(grid_.ijk_res_, (int*)num_pts_cell_);
		array_view<int, 1> start_idx_cell_view(grid_.ijk_res_, (int*)start_idx_cell_);

		BEGIN_PARALLEL_FOR_EACH_1D(grid_.ijk_res_, p)
		{
			num_pts_cell_view[p] = 0;
			start_idx_cell_view[p] = -1;
		}
		END_PARALLEL_FOR_EACH_1D

		num_pts_cell_view.synchronize();
		start_idx_cell_view.synchronize();

		BEGIN_CPU_THREADS_1D(num_of_pts_)
		{
			THREAD_LOOPS_1D(p)
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
			THREAD_LOOPS_1D(p)
			{
				const atomic<int>& num = num_pts_cell_[p];
				start_idx_cell_[p] = std::atomic_fetch_add(&start_idx, num);

				if (num > 0)
				{
					int ix = std::atomic_fetch_add(&num_of_pts_cell_, 1);
					pts_cell_index_array_[ix] = p;
				}
			}
		}
		END_CPU_THREADS_1D;

		BEGIN_CPU_THREADS_1D(num_of_pts_)
		{
			THREAD_LOOPS_1D(p)
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
			THREAD_LOOPS_1D(p)
			{
				const int& num = num_pts_cell_[p];
				const int& b_ix = start_idx_cell_[p];

				for (int n = 0; n < num; n++)
				{
					int v_ix = particle_index_array_[b_ix + n];
					
					for (int m = 0; m < (int)vector_data_.size(); m++)
					{
						Vec3T* data = vector_data_[m];
						Vec3T* data_temp = vector_data_temp_[m];
						data_temp[b_ix + n] = data[v_ix];
					}
					for (int m = 0; m < (int)scalar_data_.size(); m++)
					{
						T* data = scalar_data_[m];
						T* data_temp = scalar_data_temp_[m];
						data_temp[b_ix + n] = data[v_ix];
					}
				}
			}
		}
		END_CPU_THREADS_1D;

		for (int m = 0; m < (int)vector_data_.size(); m++)
		{
			Vec3T* temp = vector_data_temp_[m];
			vector_data_temp_[m] = vector_data_[m];
			vector_data_[m] = temp;

			(*vector_data_pointer_[m]) = temp;
		}
		for (int m = 0; m < (int)scalar_data_.size(); m++)
		{
			T* temp = scalar_data_temp_[m];
			scalar_data_temp_[m] = scalar_data_[m];
			scalar_data_[m] = temp;

			(*scalar_data_pointer_[m]) = temp;
		}

		position_array_ = vector_data_[0];
		velocity_array_ = vector_data_[1];
		
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