

#pragma once

#include "GRID_UNIFORM_3D.h"

class WALL_CONDITIONS
{
public:

	GRID_UNIFORM_3D grid_;	

public:

	void Initialize(GRID_UNIFORM_3D& grid_input)
	{
		grid_ = grid_input;
	}

	void ClampPositionAndVelocity(Vec3T& pos, Vec3T& vel, Vec3T& force)
	{
		T penalty_coef = 0;
		
		if (grid_.min_.x + grid_.gx_ + grid_.dx_*(T)0.5 >= pos.x)
		{
			pos.x = grid_.min_.x + grid_.gx_ + grid_.dx_*(T)0.5 + (T)FLT_EPSILON;
			vel.x = (T)0;
		}
		if (grid_.max_.x - grid_.gx_ - grid_.dx_*(T)0.5 <= pos.x)
		{
			pos.x = grid_.max_.x - grid_.gx_ - grid_.dx_*(T)0.5 - (T)FLT_EPSILON;
			vel.x = (T)0;
		}


		if (grid_.min_.y + grid_.gy_ + grid_.dy_*(T)0.5 >= pos.y)
		{
			pos.y = grid_.min_.y + grid_.gy_ + grid_.dy_*(T)0.5 + (T)FLT_EPSILON;
			vel.y = (T)0;
		}
		if (grid_.max_.y - grid_.gy_ - grid_.dy_*(T)0.5 <= pos.y)
		{
			pos.y = grid_.max_.y - grid_.gy_ - grid_.dy_*(T)0.5 - (T)FLT_EPSILON;
			vel.y = (T)0;
		}


		if (grid_.min_.z + grid_.gz_ + grid_.dz_*(T)0.5 >= pos.z)
		{
			pos.z = grid_.min_.z + grid_.gz_ + grid_.dz_*(T)0.5 + (T)FLT_EPSILON;
			vel.z = (T)0;
		}
		if (grid_.max_.z - grid_.gz_ - grid_.dz_*(T)0.5 <= pos.z)
		{
			pos.z = grid_.max_.z - grid_.gz_ - grid_.dz_*(T)0.5 - (T)FLT_EPSILON;
			vel.z = (T)0;
		}
	}
};