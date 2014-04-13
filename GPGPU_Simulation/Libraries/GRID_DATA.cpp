
#include <GL\glut.h>
#include "GRID_DATA.h"

template<class TT>
void GRID_DYNAMIC<TT>::Initialize(const Vec3& min, const Vec3& max, int i, int j, int k, int g, int block_i)
{
	ghost_width_ = g;

	i_res_ = i + g*2;
	j_res_ = j + g*2;
	k_res_ = k + g*2;

	block_i_res_ = MIN(block_i, i_res_);
	block_size_  = i_res_/block_i_res_;	

	i_res_ = block_i_res_ * block_size_;
	i = i_res_ - g*2;

	ij_res_  = i_res_ * j_res_;
	ijk_res_ = i_res_ * j_res_ * k_res_;

	//

	dx_ = (max.x - min.x)/(FLT)i;
	dy_ = (max.y - min.y)/(FLT)j;
	dz_ = (max.z - min.z)/(FLT)k;

	one_over_dx_ = (FLT)1/dx_;
	one_over_dy_ = (FLT)1/dy_;
	one_over_dz_ = (FLT)1/dz_;

	gx_ = dx_*(FLT)g;
	gy_ = dy_*(FLT)g;
	gz_ = dz_*(FLT)g;

	min_ = min - Vec3(gx_, gy_, gz_);
	max_ = max + Vec3(gx_, gy_, gz_);

	compressed_arr_ = new TT*[block_i_res_*j_res_*k_res_];
	memset(compressed_arr_, 0, sizeof(TT*)*block_i_res_*j_res_*k_res_);

	num_blocks_ = ij_res_;

	for(int n=0; n<num_blocks_; n++)
	{
		TT* new_block = new TT[block_size_];
		free_blocks_.push(new_block);
	}
}

template<class TT>
void GRID_DYNAMIC<TT>::Finalize()
{

}

template<class TT>
void GRID_DYNAMIC<TT>::Render()
{
	glDisable(GL_LIGHTING);
	glPointSize(2.0f);
	glBegin(GL_POINTS);
	glColor3f(1.0f, 0.0f, 0.0f);
	for(int k = 0; k < k_res_; k++) for (int j = 0; j < j_res_; j++) for(int b=0; b<block_i_res_; b++)
	{
		int b_ix = k*(block_i_res_*j_res_) + j*block_i_res_ + b;
		TT* block = compressed_arr_[b_ix];

		if(block == 0) continue;

		int n = b*block_size_;

		for(int p=0; p<block_size_; p++)
		{
			const int i = n+p;

			Vec3 pos = CellCenterPosition(i,j,k);

			glVertex3f(pos.x, pos.y, pos.z);
		}				
	}
	glEnd();

	glPointSize(1.0f);
	glEnable(GL_LIGHTING);
}

template class GRID_DYNAMIC<FLT>;