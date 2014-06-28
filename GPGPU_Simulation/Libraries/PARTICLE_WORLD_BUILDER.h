
#include "GRID.h"
#include "MATH_CORE.h"
#include "KERNEL_FUNCTIONS.h"

using namespace std;
using namespace concurrency;

class PARTICLE_WORLD_BUILDER
{
public:

	GRID world_grid_;

	TS*  density_field_;
	TV3* velocity_field_;
	TV3* normal_field_;

	atomic<int>*  pts_num_buffer_;
	atomic<int>*  start_idx_buffer_;

	atomic<int>*  index_array_;
	atomic<int>*  id_array_;

public:

	PARTICLE_WORLD_BUILDER() : 
		density_field_(0), velocity_field_(0), normal_field_(0), pts_num_buffer_(0), start_idx_buffer_(0), index_array_(0), id_array_(0)
	{}
	~PARTICLE_WORLD_BUILDER()
	{}

	void Initialize(const GRID& grid)
	{
		world_grid_ = grid;

		density_field_  = new TS[world_grid_.ijk_res_];
		velocity_field_ = new TV3[world_grid_.ijk_res_];
		normal_field_   = new TV3[world_grid_.ijk_res_];

		pts_num_buffer_   = new atomic<int>[world_grid_.ijk_res_];
		start_idx_buffer_ = new atomic<int>[world_grid_.ijk_res_];
	}

	void RasterizeParticles(const TV3* pos_arr, const TV3* vel_arr, const TS pts_mass, const int num_pts)
	{
		index_array_ = new atomic<int>[num_pts];
		id_array_ = new atomic<int>[num_pts];

		BuildParticleDataStructure(pos_arr, num_pts);

		#pragma omp parallel for
		for(int p=0; p<world_grid_.ijk_res_; p++)
		{
			int i, j, k;
			world_grid_.Index1Dto3D(p, i, j, k);

			TV3 cell_center = world_grid_.CellCenterPosition(i,j,k);

			TV3 velocity_weighted = TV3();
			TV3 normal_weighted = TV3();
			TS mass_weighted = (TS)0;

			BEGIN_STENCIL_LOOP(world_grid_, i, j, k, l, m, n)
			{
				int s_ix = world_grid_.Index3Dto1D(l, m, n);

				int b_ix = start_idx_buffer_[s_ix];
				int num = pts_num_buffer_[s_ix];

				for(int x=0; x<num; x++)
				{
					int v_ix = index_array_[b_ix+x];

					const TV3& pos = pos_arr[v_ix];
					const TV3& vel = vel_arr[v_ix];

					TS w = QuadBSplineKernel(pos-cell_center, world_grid_.one_over_dx_, world_grid_.one_over_dx_, world_grid_.one_over_dx_);		
					TV3 grad = QuadBSplineKernelGradient(pos-cell_center, world_grid_.one_over_dx_, world_grid_.one_over_dx_, world_grid_.one_over_dx_);	

					mass_weighted += pts_mass * w;
					velocity_weighted += vel * pts_mass * w;
					normal_weighted += grad * pts_mass;
				}
			}
			END_STENCIL_LOOP;

			density_field_[p] = mass_weighted;
			velocity_field_[p] = velocity_weighted / (mass_weighted + FLT_EPSILON);
			normal_field_[p] = normal_weighted;
		}
	}

	void BuildParticleDataStructure(const TV3* pos_arr, const int num_pts)
	{
		#pragma omp parallel for
		for(int p=0; p<world_grid_.ijk_res_; p++)
		{
			pts_num_buffer_[p] = 0;
			start_idx_buffer_[p] = -1;
		}

		#pragma omp parallel for
		for(int p=0; p<num_pts; p++)
		{
			int i, j, k, ix;
			if(world_grid_.IsInsideGhost(pos_arr[p]) == false) continue;

			world_grid_.CellCenterIndex(pos_arr[p], i, j, k);
			ix = world_grid_.Index3Dto1D(i, j, k);

			id_array_[p] = atomic_fetch_add(&pts_num_buffer_[ix], 1);
		}

		atomic<int> start_idx = 0;		

		#pragma omp parallel for
		for(int p=0; p<world_grid_.ijk_res_; p++)
		{
			const atomic<int>& num = pts_num_buffer_[p];
			start_idx_buffer_[p] = std::atomic_fetch_add(&start_idx, num);
		}

		#pragma omp parallel for
		for(int p=0; p<num_pts; p++)
		{
			int i, j, k, ix;
			if(world_grid_.IsInsideGhost(pos_arr[p]) == false) continue;

			world_grid_.CellCenterIndex(pos_arr[p], i, j, k);
			ix = world_grid_.Index3Dto1D(i, j, k);

			int s_ix = start_idx_buffer_[ix];
			int id = id_array_[p];

			index_array_[s_ix + id] = p;
		}
	}

	void Render()
	{
		glDisable(GL_LIGHTING);

		glPushMatrix();
		glColor3f(0, 1, 0);

		glPointSize(3);

		glBegin(GL_POINTS);
		for (int c = 0; c < world_grid_.ijk_res_; c++)
		{
			int i, j, k;
			world_grid_.Index1Dto3D(c, i, j, k);

			TV3 cell_center = world_grid_.CellCenterPosition(i, j, k);

			const TS& density = density_field_[c];

			if (density > FLT_EPSILON)
			{
				glVertex3f(cell_center.x, cell_center.y, cell_center.z);
			}
		}
		glEnd();
/*
		glColor3f(1, 0, 0);
		glBegin(GL_LINES);
		for (int c = 0; c < world_grid_.ijk_res_; c++)
		{
			int i, j, k;
			world_grid_.Index1Dto3D(c, i, j, k);

			TV3 cell_center = world_grid_.CellCenterPosition(i, j, k);

			const TS& density = density_field_[c];

			const TV3 n = normal_field_[c];
			const TS d = (TS)0.01;

			if(density > FLT_EPSILON)
			{
				glVertex3f(cell_center.x, cell_center.y, cell_center.z);
				glVertex3f(cell_center.x+n.x*d, cell_center.y+n.y*d, cell_center.z+n.z*d);
			}
		}
		glEnd();
*/
		glPopMatrix();
	}
};