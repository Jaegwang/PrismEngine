
#pragma once

#include <amp.h>
#include "GRID_UNIFORM_3D.h"
#include "PARTICLE.h"
#include "VECTOR3_T.h"
#include "GL\glut.h"
#include "iostream"

using namespace std;
using namespace concurrency;

class PARTICLE_MANAGER_3D
{
public:

	GRID_UNIFORM_3D grid_;

	// particle properties
	Vec3T *position_array_, *velocity_array_;
	Vec3T *position_array_old_, *velocity_array_old_;

	T *density_array_;
	T *density_array_old_;

	int* particle_id_array_;
	int* particle_index_array_;

	int num_of_pts_;
	int max_of_pts_;

	int* num_pts_cell_;
	int* start_idx_cell_;

public:

	PARTICLE_MANAGER_3D() : 
		position_array_(0), position_array_old_(0), velocity_array_(0), velocity_array_old_(0), density_array_(0), density_array_old_(0),
		num_pts_cell_(0), start_idx_cell_(0), particle_id_array_(0), particle_index_array_(0), max_of_pts_(0)
	{}

	~PARTICLE_MANAGER_3D()
	{
		if (position_array_) delete[] position_array_;
		if (position_array_old_) delete[] position_array_old_;

		if (velocity_array_) delete[] velocity_array_;
		if (velocity_array_old_) delete[] velocity_array_old_;

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

		density_array_ = new T[num_pts];
		density_array_old_ = new T[num_pts];


		particle_id_array_ = new int[num_pts];
		particle_index_array_ = new int[num_pts];

		num_pts_cell_ = new int[grid_.ijk_res_];
		start_idx_cell_ = new int[grid_.ijk_res_];



		memset((void*)density_array_, 0, sizeof(T)*max_of_pts_);
		memset((void*)density_array_old_, 0, sizeof(T)*max_of_pts_);
	}

	int AddParticle()
	{
		int ix = num_of_pts_++;
		return ix;
	}

	void DelParticle(int ix)
	{
		T max_pos = (T)100000;
		position_array_[ix] = Vec3T(max_pos, max_pos, max_pos);
	}

	void RebuildParticleDataStructure()
	{
		if (num_of_pts_ == 0) return;
		
		array_view<Vec3T, 1> pos_arr_view(num_of_pts_, position_array_);
		array_view<Vec3T, 1> pos_arr_old_view(num_of_pts_, position_array_old_);

		array_view<Vec3T, 1> vel_arr_view(num_of_pts_, velocity_array_);
		array_view<Vec3T, 1> vel_arr_old_view(num_of_pts_, velocity_array_old_);

		array_view<T, 1> den_arr_view(num_of_pts_, density_array_);
		array_view<T, 1> den_arr_old_view(num_of_pts_, density_array_old_);



		array_view<int, 1> pts_id_view(num_of_pts_, particle_id_array_);
		array_view<int, 1> pts_index_view(num_of_pts_, particle_index_array_);
		
		array_view<int, 1> num_pts_cell_view(grid_.ijk_res_, num_pts_cell_);
		array_view<int, 1> start_idx_cell_view(grid_.ijk_res_, start_idx_cell_);		
		
		int start_idx = 0; array_view<int, 1> start_idx_view(1, &start_idx);

		int valid_num_pts = 0; array_view<int, 1> valid_num_pts_view(1, &valid_num_pts);

		array_view<GRID_UNIFORM_3D, 1> grid_view(1, &grid_);

		parallel_for_each(num_pts_cell_view.extent, [=](index<1> idx) restrict(amp)
		{
			num_pts_cell_view[idx] = 0;
			start_idx_cell_view[idx] = -1;
		});
		
		parallel_for_each(pos_arr_view.extent, [=](index<1> idx) restrict(amp)
		{
			Vec3T& position = (Vec3T)pos_arr_view[idx[0]];
			int& pts_id = pts_id_view[idx[0]];

			if (grid_view[0].IsInside(position) == false) return;

			int i, j, k;
			grid_view[0].CellCenterIndex(position, i, j, k);
			int ix = grid_view[0].Index3Dto1D(i, j, k);
			
			pts_id = atomic_fetch_add(&num_pts_cell_view[ix], 1);

			atomic_fetch_add(&valid_num_pts_view[0], 1);
		});

		parallel_for_each(num_pts_cell_view.extent, [=](index<1> idx) restrict(amp)
		{
			int num = num_pts_cell_view[idx];
			start_idx_cell_view[idx] = atomic_fetch_add(&start_idx_view[0], num);
		});
	
		parallel_for_each(pos_arr_view.extent, [=](index<1> idx) restrict(amp)
		{
			Vec3T& position = pos_arr_view[idx[0]];
			int& pts_id = pts_id_view[idx[0]];

			if (grid_view[0].IsInside(position) == false) return;

			int i, j, k;
			grid_view[0].CellCenterIndex(position, i, j, k);

			int ix = grid_view[0].Index3Dto1D(i, j, k);

			int b_ix = start_idx_cell_view[ix];

			pts_index_view[b_ix + pts_id] = idx[0];
		});

		parallel_for_each(num_pts_cell_view.extent, [=](index<1> idx) restrict(amp)
		{
			int num = num_pts_cell_view[idx];
			int b_ix = start_idx_cell_view[idx];

			for (int n = 0; n < num; n++)
			{
				int v_ix = pts_index_view[b_ix + n];
				pos_arr_old_view[b_ix + n] = pos_arr_view[v_ix];

			}
		});
		
		parallel_for_each(num_pts_cell_view.extent, [=](index<1> idx) restrict(amp)
		{
			int num = num_pts_cell_view[idx];
			int b_ix = start_idx_cell_view[idx];

			for (int n = 0; n < num; n++)
			{
				int v_ix = pts_index_view[b_ix + n];
				vel_arr_old_view[b_ix + n] = vel_arr_view[v_ix];
				den_arr_old_view[b_ix + n] = den_arr_view[v_ix];
			}
		});
		

		num_pts_cell_view.synchronize();
		start_idx_cell_view.synchronize();

		valid_num_pts_view.synchronize();
		pos_arr_view.synchronize();
		vel_arr_view.synchronize();
		den_arr_view.synchronize();

		Vec3T* v_temp = 0; T* t_temp = 0;

		SWAP(position_array_, position_array_old_, v_temp);
		SWAP(velocity_array_, velocity_array_old_, v_temp);
		SWAP(density_array_, density_array_old_, t_temp);

		num_of_pts_ = valid_num_pts;
	}

	void Rendering()
	{
		glDisable(GL_LIGHTING);

		glPushMatrix();
		glColor3f(0, 0, 0);

		glBegin(GL_POINTS);

		int count = 0;

		for (int n = 0; n < grid_.ijk_res_; n++)
		{
			int num = num_pts_cell_[n];
			int b_ix = start_idx_cell_[n];


			for (int x = 0; x < num; x++)
			{
				Vec3T pos = position_array_[b_ix + x];

				glVertex3f(pos.x, pos.y, pos.z);

				count++;
			}
		}

		glEnd();

		glPopMatrix();

		if (count != num_of_pts_)
		{
			std::cout << "Error" << std::endl;
		}

		std::cout << num_of_pts_ << std::endl;
	}
};