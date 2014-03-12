

#pragma once
#include <GL\glut.h>
#include <amp.h>
#include "MATH_CORE.h"

using namespace concurrency;

class GRID_UNIFORM_3D
{
public:
	
	Vec3 min_, max_;

	FLT dx_, dy_, dz_;

	FLT gx_, gy_, gz_;

	FLT one_over_dx_, one_over_dy_, one_over_dz_;

	int i_res_, j_res_, k_res_;
	int ij_res_, ijk_res_;

	int ghost_width_;
	
public:
	
	GRID_UNIFORM_3D()
	{}
	~GRID_UNIFORM_3D()
	{}

	void Initialize(const Vec3& min, const Vec3& max, const int i_res, const int j_res, const int k_res, const int g)
	{
		dx_ = (max.x-min.x)/(FLT)i_res;
		dy_ = (max.y-min.y)/(FLT)j_res;
		dz_ = (max.z-min.z)/(FLT)k_res;

		one_over_dx_ = (FLT)1 / dx_;
		one_over_dy_ = (FLT)1 / dy_;
		one_over_dz_ = (FLT)1 / dz_;

		gx_ = dx_*(FLT)g;
		gy_ = dy_*(FLT)g;
		gz_ = dz_*(FLT)g;
		
		min_ = min - Vec3(gx_, gy_, gz_);
		max_ = max + Vec3(gx_, gy_, gz_);

		ghost_width_ = g;

		i_res_ = i_res + g*2;
		j_res_ = j_res + g*2;
		k_res_ = k_res + g*2;

		ij_res_  = i_res_ * j_res_;
		ijk_res_ = i_res_ * j_res_ * k_res_;
	}

	int Index3Dto1D(const int i, const int j, const int k) 
	{
		return i + j*i_res_ + k*ij_res_;
	}

	void Index1Dto3D(const int idx, int& i, int& j, int& k) 
	{
		i = idx           % i_res_;
		j = (idx/i_res_ ) % j_res_;
		k = (idx/ij_res_) % k_res_;
	}

	bool IsGhostCell(const int i, const int j, const int k) 
	{
		if(i < ghost_width_ || i >= i_res_-ghost_width_) return true;
		if(j < ghost_width_ || j >= j_res_-ghost_width_) return true;
		if(k < ghost_width_ || k >= k_res_-ghost_width_) return true;

		return false;
	}

	bool IsInsideValid(const Vec3& p)
	{
		if (min_.x+gx_ < p.x && max_.x-gx_ > p.x &&
			min_.y+gy_ < p.y && max_.y-gy_ > p.y &&
			min_.z+gz_ < p.z && max_.z-gz_ > p.z)
			return true;

		return false;
	}

	bool IsInsideGhost(const Vec3& p)
	{
		if (min_.x < p.x && max_.x > p.x &&
			min_.y < p.y && max_.y > p.y &&
			min_.z < p.z && max_.z > p.z)
			return true;

		return false;
	}

	Vec3 CellCenterPosition(const int i, const int j, const int k) 
	{
		return min_ + Vec3(((FLT)i+(FLT)0.5)*dx_, ((FLT)j+(FLT)0.5)*dy_, ((FLT)k+(FLT)0.5)*dz_);
	}

	void CellCenterIndex(const Vec3& p, int& i, int& j, int& k) 
	{
		Vec3 v = p-min_;
		i = (int)(v.x/dx_);
		j = (int)(v.y/dy_);
		k = (int)(v.z/dz_);
	}

	void LeftBottomIndex(const Vec3& p, int& i, int& j, int& k)  
	{
		Vec3 v = p-min_;
		i = (int)((v.x/dx_)-(FLT)0.5);
		j = (int)((v.y/dy_)-(FLT)0.5);
		k = (int)((v.z/dz_)-(FLT)0.5);
	}

	void RightUpIndex(const Vec3& p, int& i, int& j, int& k) 
	{
		Vec3 v = p-min_;
		i = (int)((v.x/dx_)+(FLT)0.5); 
		j = (int)((v.y/dy_)+(FLT)0.5);
		k = (int)((v.z/dz_)+(FLT)0.5);
	}

	void StartEndIndices(const int l, const int m, const int n, int& start_l, int& start_m, int& start_n, int& end_l, int& end_m, int& end_n, const int pad=1) 
	{
		start_l = MAX(l-pad, 0);
		start_m = MAX(m-pad, 0);
		start_n = MAX(n-pad, 0);

		end_l = MIN(l+pad, i_res_-1);
		end_m = MIN(m+pad, j_res_-1);
		end_n = MIN(n+pad, k_res_-1);
	}

	void StencilIndexBuffer(const int i, const int j, const int k, int* buf)
	{
		int start_l, start_m, start_n, end_l, end_m, end_n;
		StartEndIndices(i, j, k, start_l, start_m, start_n, end_l, end_m, end_n, 1);

		int ix = 0;

		for (int n = start_n; n <= end_n; n++) for (int m = start_m; m <= end_m; m++) for (int l = start_l; n <= end_l; l++)
		{
			buf[ix++] = Index3Dto1D(l, m, n);
		}
	}

	template<class TT>
	TT TriLinearInterpolate(const Vec3& p, TT* arr) 
	{ //http://en.wikipedia.org/wiki/Trilinear_interpolation

		int b_i, b_j, b_k;

		LeftBottomIndex(p, b_i, b_j, b_k);
		int ix = Index3Dto1D(b_i, b_j, b_k);

		Vec3 p_0 = CellCenterPosition(b_i ,b_j ,b_k);

		FLT x_d = (p.x-p_0.x)/dx_;
		FLT y_d = (p.y-p_0.y)/dy_;
		FLT z_d = (p.z-p_0.z)/dz_;

		TT c_00 = arr[ix]*((FLT)1-x_d) + arr[ix+1]*(x_d);
		TT c_10 = arr[ix+i_res_]*((FLT)1-x_d) + arr[ix+i_res_+1]*(x_d);
		TT c_01 = arr[ix+ij_res_]*((FLT)1-x_d) + arr[ix+ij_res_+1]*(x_d);
		TT c_11 = arr[ix+i_res_+ij_res_]*((FLT)1-x_d) + arr[ix+i_res_+ij_res_+1]*(x_d);

		TT c_0 = c_00*((FLT)1-y_d) + c_10*y_d;
		TT c_1 = c_01*((FLT)1-y_d) + c_11*y_d;

		return c_0*((FLT)1-z_d) + c_1*z_d;
	}


	void RenderGrid()
	{
		glDisable(GL_LIGHTING);

		Vec3 grid_center = (min_+max_)*(FLT)0.5;
		Vec3 deviation = max_-min_;

		glColor3f(0.5, 0.5, 0.5);

		glPushMatrix();

		glTranslatef(grid_center.x, grid_center.y, grid_center.z);
		glScaled(deviation.x, deviation.y, deviation.z);
		glutWireCube(1);
		glPopMatrix();
		
		deviation -= Vec3(gx_*(FLT)2, gy_*(FLT)2, gz_*(FLT)2);

		glColor3f(0, 0, 0);

		glPushMatrix();

		glTranslatef(grid_center.x, grid_center.y, grid_center.z);
		glScaled(deviation.x, deviation.y, deviation.z);
		glutWireCube(1);
		glPopMatrix();
				
		glEnable(GL_LIGHTING);
	}

	void RenderCells()
	{
		glDisable(GL_LIGHTING);

		glPushMatrix();
		glColor3f(0, 0, 0);
		
		glBegin(GL_POINTS);
		for(int x=0; x<ijk_res_; x++)
		{
			int i,j,k;
			Index1Dto3D(x,i,j,k);

			if(!IsGhostCell(i,j,k))
			{
				Vec3 cell_center = CellCenterPosition(i,j,k);

				glVertex3f(cell_center.x, cell_center.y, cell_center.z);
			}
		}
		glEnd();

		glPopMatrix();
		glEnable(GL_LIGHTING);
	}

};


#define BEGIN_STENCIL_LOOP(_grid, _i, _j, _k, _l, _m, _n) { int _start_l, _start_m, _start_n, _end_l, _end_m, _end_n; \
															_grid.StartEndIndices(_i, _j, _k, _start_l, _start_m, _start_n, _end_l, _end_m, _end_n); \
															int _b_ix = _grid.Index3Dto1D(_start_l, _start_m, _start_n); \
															for (int _ns = 0, _n = _start_n; _n <= _end_n; _ns += _grid.ij_res_, _n++)\
															for (int _ms = 0, _m = _start_m; _m <= _end_m; _ms += _grid.i_res_, _m++)\
															for (int _ls = 0, _l = _start_l; _l <= _end_l; _ls += 1, _l++) { int _p = _b_ix + _ns + _ms + _ls;
#define END_STENCIL_LOOP }}