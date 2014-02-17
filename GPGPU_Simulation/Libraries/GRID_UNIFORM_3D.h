

#pragma once
#include <GL\glut.h>
#include <amp.h>
#include "VECTOR3_T.h"

using namespace concurrency;

class GRID_UNIFORM_3D
{
public:
	
	Vec3T min_, max_;

	T dx_, dy_, dz_;

	T gx_, gy_, gz_;

	T one_over_dx_, one_over_dy_, one_over_dz_;

	int i_res_, j_res_, k_res_;
	int ij_res_, ijk_res_;

	int ghost_width_;
	int g_;
	
public:
	
	GRID_UNIFORM_3D()
	{}
	~GRID_UNIFORM_3D()
	{}

	void Initialize(const Vec3T& min, const Vec3T& max, const int i_res, const int j_res, const int k_res, const int g) restrict(cpu,amp)
	{
		dx_ = (max.x-min.x)/(T)i_res;
		dy_ = (max.y-min.y)/(T)j_res;
		dz_ = (max.z-min.z)/(T)k_res;

		one_over_dx_ = (T)1 / dx_;
		one_over_dy_ = (T)1 / dy_;
		one_over_dz_ = (T)1 / dz_;

		gx_ = dx_*(T)g;
		gy_ = dy_*(T)g;
		gz_ = dz_*(T)g;
		
		min_ = min - Vec3T(gx_, gy_, gz_);
		max_ = max + Vec3T(gx_, gy_, gz_);

		ghost_width_ = g;
		g_ = g;

		i_res_ = i_res + g*2;
		j_res_ = j_res + g*2;
		k_res_ = k_res + g*2;

		ij_res_  = i_res_ * j_res_;
		ijk_res_ = i_res_ * j_res_ * k_res_;
	}

	void UseGhostStartIndex(bool b)
	{
		if(b == true) g_ = 0;
		else g_ = ghost_width_;
	}

	int Index3Dto1D(const int i, const int j, const int k) restrict(cpu,amp)
	{
		return (i+g_) + (j+g_)*i_res_ + (k+g_)*ij_res_;
	}

	void Index1Dto3D(const int idx, int& i, int& j, int& k) restrict(cpu,amp)
	{
		i = idx           % i_res_ - g_;
		j = (idx/i_res_ ) % j_res_ - g_;
		k = (idx/ij_res_) % k_res_ - g_;
	}

	bool IsGhostCell(const int i, const int j, const int k) restrict(cpu,amp)
	{
		if((i+g_) < ghost_width_ || (i+g_) >= i_res_-ghost_width_) return true;
		if((j+g_) < ghost_width_ || (j+g_) >= j_res_-ghost_width_) return true;
		if((k+g_) < ghost_width_ || (k+g_) >= k_res_-ghost_width_) return true;

		return false;
	}

	bool IsInsideValid(const Vec3T& p) restrict(cpu, amp)
	{
		if (min_.x+gx_ < p.x && max_.x-gx_ > p.x &&
			min_.y+gy_ < p.y && max_.y-gy_ > p.y &&
			min_.z+gz_ < p.z && max_.z-gz_ > p.z)
			return true;

		return false;
	}

	bool IsInsideGhost(const Vec3T& p) restrict(cpu, amp)
	{
		if (min_.x < p.x && max_.x > p.x &&
			min_.y < p.y && max_.y > p.y &&
			min_.z < p.z && max_.z > p.z)
			return true;

		return false;
	}

	Vec3T CellCenterPosition(const int i, const int j, const int k) restrict(cpu,amp)
	{
		return min_ + Vec3T(((T)(i+g_)+(T)0.5)*dx_, ((T)(j+g_)+(T)0.5)*dy_, ((T)(k+g_)+(T)0.5)*dz_);
	}

	void CellCenterIndex(const Vec3T& p, int& i, int& j, int& k) restrict(cpu,amp)
	{
		Vec3T v = p-min_;
		i = (int)(v.x/dx_) - g_;
		j = (int)(v.y/dy_) - g_;
		k = (int)(v.z/dz_) - g_;
	}

	void LeftBottomIndex(const Vec3T& p, int& i, int& j, int& k) restrict(cpu,amp) 
	{
		Vec3T v = p-min_;
		i = (int)((v.x/dx_)-(T)0.5) - g_;
		j = (int)((v.y/dy_)-(T)0.5) - g_;
		k = (int)((v.z/dz_)-(T)0.5) - g_;
	}

	void RightUpIndex(const Vec3T& p, int& i, int& j, int& k) restrict(cpu,amp)
	{
		Vec3T v = p-min_;
		i = (int)((v.x/dx_)+(T)0.5) - g_; 
		j = (int)((v.y/dy_)+(T)0.5) - g_;
		k = (int)((v.z/dz_)+(T)0.5) - g_;
	}

	void StartEndIndices(const int l, const int m, const int n, int& start_l, int& start_m, int& start_n, int& end_l, int& end_m, int& end_n, const int pad=1) restrict(cpu,amp)
	{
		start_l = MAX(l-pad, -g_);
		start_m = MAX(m-pad, -g_);
		start_n = MAX(n-pad, -g_);

		end_l = MIN(l+pad, i_res_-g_-1);
		end_m = MIN(m+pad, j_res_-g_-1);
		end_n = MIN(n+pad, k_res_-g_-1);
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
	TT TriLinearInterpolate(const Vec3T& p, TT* arr) restrict(cpu,amp)
	{ //http://en.wikipedia.org/wiki/Trilinear_interpolation

		int b_i, b_j, b_k;

		LeftBottomIndex(p, b_i, b_j, b_k);
		int ix = Index3Dto1D(b_i, b_j, b_k);

		Vec3 p_0 = CellCenter(b_i ,b_j ,b_k);

		T_ x_d = (p.x-p_0.x)/dx_;
		T_ y_d = (p.y-p_0.y)/dy_;
		T_ z_d = (p.z-p_0.z)/dz_;

		TT c_00 = arr[ix]*((T_)1-x_d) + arr[ix+1]*(x_d);
		TT c_10 = arr[ix+i_res_]*((T_)1-x_d) + arr[ix+i_res_+1]*(x_d);
		TT c_01 = arr[ix+ij_res_]*((T_)1-x_d) + arr[ix+ij_res_+1]*(x_d);
		TT c_11 = arr[ix+i_res_+ij_res_]*((T_)1-x_d) + arr[ix+i_res_+ij_res_+1]*(x_d);

		TT c_0 = c_00*((T_)1-y_d) + c_10*y_d;
		TT c_1 = c_01*((T_)1-y_d) + c_11*y_d;

		return c_0*((T_)1-z_d) + c_1*z_d;
	}


	void RenderGrid()
	{
		glDisable(GL_LIGHTING);

		Vec3T grid_center = (min_ + max_)*0.5;
		Vec3T deviation = max_-min_;

		glColor3f(0.5, 0.5, 0.5);

		glPushMatrix();

		glTranslatef(grid_center.x, grid_center.y, grid_center.z);
		glScaled(deviation.x, deviation.y, deviation.z);
		glutWireCube(1);
		glPopMatrix();
		
		deviation -= Vec3T(gx_*(T)2, gy_*(T)2, gz_*(T)2);

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
				Vec3T cell_center = CellCenterPosition(i,j,k);

				glVertex3f(cell_center.x, cell_center.y, cell_center.z);
			}
		}
		glEnd();

		glPopMatrix();
		glEnable(GL_LIGHTING);
	}

};