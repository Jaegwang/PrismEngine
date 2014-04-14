
#include <GL\glut.h>
#include "GRID_DATA.h"

void GRID_DYNAMIC::Initialize(const Vec3& min, const Vec3& max, const FLT dw, const int g, const int block_size)
{
	ghost_width_ = g;
	block_size_  = block_size;

	dw_ = dw;
	gw_ = dw*(FLT)g;

	i_res_ = (max.x - min.x + dw)/dw + g*2;
	j_res_ = (max.y - min.y + dw)/dw + g*2;
	k_res_ = (max.z - min.z + dw)/dw + g*2;

	b_res_ = (i_res_+block_size_-1)/block_size_;
	i_res_ = b_res_*block_size_;

	min_ = min - Vec3(gw_, gw_, gw_);
	max_ = min + Vec3((FLT)i_res_*dw_, (FLT)j_res_*dw_, (FLT)k_res_*dw_);

	compressed_arr_ = new int[b_res_*j_res_*k_res_];

	for(int n=0; n<b_res_*j_res_*k_res_; n++) compressed_arr_[n] = -1;
}

void GRID_DYNAMIC::Initialize(const Vec3& min, const Vec3& max, int i, int j, int k, int g, int block_i)
{
	ghost_width_ = g;

	i_res_ = i + g*2;
	j_res_ = j + g*2;
	k_res_ = k + g*2;
	b_res_ = MIN(block_i, i_res_);

	block_size_  = i_res_/b_res_;	

	i_res_ = b_res_ * block_size_;
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

	compressed_arr_ = new int[b_res_*j_res_*k_res_];

	for(int n=0; n<b_res_*j_res_*k_res_; n++)
		compressed_arr_[n] = -1;
}

void GRID_DYNAMIC::Finalize()
{

}

void GRID_DYNAMIC::Render()
{
	glDisable(GL_LIGHTING);
	glPointSize(2.0f);
	glBegin(GL_POINTS);
	glColor3f(1.0f, 0.0f, 0.0f);

	for(int k = 0; k < k_res_; k++) for (int j = 0; j < j_res_; j++) for(int b=0; b<b_res_; b++)
	{
		int b_ix = k*(b_res_*j_res_) + j*b_res_ + b;
		int start_idx = compressed_arr_[b_ix];

		if(start_idx < 0) continue;

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