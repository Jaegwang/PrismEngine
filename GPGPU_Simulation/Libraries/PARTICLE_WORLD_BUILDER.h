
#include "GRID_UNIFORM_3D.h"
#include "MATRIX3_T.h"

class PARTICLE_WORLD_BUILDER
{
public:

	GRID_UNIFORM_3D world_grid_;

	T*     density_field_;
	Vec3T* velocity_field_;

public:

	PARTICLE_WORLD_BUILDER()
	{}
	~PARTICLE_WORLD_BUILDER()
	{}

	void Initialize(const GRID_UNIFORM_3D& grid)
	{
		world_grid_ = grid;

		density_field_  = new T[world_grid_.ijk_res_];
		velocity_field_ = new Vec3T[world_grid_.ijk_res_];
	}

	void RasterizeObject();
};