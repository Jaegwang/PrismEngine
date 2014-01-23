
#pragma once

#include "GRID_UNIFORM_3D.h"

template<class TT>
class PARTICLE_MANAGER_3D
{
public:

	GRID_UNIFORM_3D grid_;

	TT* particle_array_;
	TT* particle_array_old_;

	int* num_pts_cell_;
	int* start_idx_cell_;

public:

	PARTICLE_MANAGER_3D() : particle_array_(0), particle_array_old_(0), num_pts_cell_(0), start_idx_cell_(0)
	{}

	~PARTICLE_MANAGER_3D()
	{
		if(particle_array_) delete[] particle_array_;
		if(particle_array_old_) delete[] particle_array_old_;

		if(num_pts_cell_) delete[] num_pts_cell_;
		if(start_idx_cell_) delete[] start_idx_cell_;
	}

	void Initialize(const GRID_UNIFORM_3D& grid_input, const int num_pts)
	{
		grid_ = grid_input;

		particle_array_ = new TT[num_pts];
		particle_array_old_ = new TT[num_pts];

		num_pts_cell_ = new int[grid_.ijk_res_];
		start_idx_cell_ = new int[grid_.ijk_res_];
	}

	void RebuildParticleDataStructure()
	{

	}


};