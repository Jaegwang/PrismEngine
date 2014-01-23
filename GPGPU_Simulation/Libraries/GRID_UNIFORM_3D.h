

#pragma once
#include <GL\glut.h>
#include "VECTOR3_T.h"

class GRID_UNIFORM_3D
{
public:

	Vec3T min_, max_;

	T dx_, dy_, dz_;

	T gx_, gy_, gz_;

	int i_res_, j_res_ ,k_res_;
	int ij_res_, ijk_res_;

	int ghost_width_;	
	
public:
	
	GRID_UNIFORM_3D()
	{}
	~GRID_UNIFORM_3D()
	{}

	void Initialize(const Vec3T& min, const Vec3T& max, const int i_res, const int j_res, const int k_res, const int g) restrict(cpu,amp)
	{
		dx_ = (max.x-min.x)/(float)i_res;
		dy_ = (max.y-min.y)/(float)j_res;
		dz_ = (max.z-min.z)/(float)k_res;

		gx_ = dx_*(float)g;
		gy_ = dy_*(float)g;
		gz_ = dz_*(float)g;
		
		min_ = min - Vec3T(gx_, gy_, gz_);
		max_ = max + Vec3T(gx_, gy_, gz_);

		ghost_width_ = g;

		i_res_ = i_res + g*2;
		j_res_ = j_res + g*2;
		k_res_ = k_res + g*2;

		ij_res_ = i_res_ * j_res_;
		ijk_res_ = i_res_ * j_res_ * k_res_;
	}

	int Index3Dto1D(const int i, const int j, const int k) restrict(cpu,amp)
	{
		return (i+ghost_width_) + (j+ghost_width_)*i_res_ + (k+ghost_width_)*ij_res_;
	}

	void Index1Dto3D(const int idx, int& i, int& j, int& k) restrict(cpu,amp)
	{
		i = idx           %i_res_ - ghost_width_;
		j = (idx/i_res_)  %j_res_ - ghost_width_;
		k = (idx/ij_res_) %k_res_ - ghost_width_;
	}

	bool IsGhostCell(const int i, const int j, const int k) restrict(cpu,amp)
	{
		if(i < 0 || i >= i_res_-ghost_width_*2) return true;
		if(j < 0 || j >= j_res_-ghost_width_*2) return true;
		if(k < 0 || k >= k_res_-ghost_width_*2) return true;

		return false;
	}

	Vec3T CellCenterPosition(const int i, const int j, const int k) restrict(cpu,amp)
	{
		return min_ + Vec3T(((T)(i+ghost_width_)+(T)0.5)*dx_, ((T)(j+ghost_width_)+(T)0.5)*dy_, ((T)(k+ghost_width_)+(T)0.5)*dz_);
	}

	void CellCenterIndex(const Vec3T& p, int& i, int& j, int& k) restrict(cpu,amp)
	{
		Vec3T v = p-min_;
		i = (int)(v.x/dx_) - ghost_width_; 
		j = (int)(v.y/dy_) - ghost_width_; 
		k = (int)(v.z/dz_) - ghost_width_; 
	}

	void LeftBottomIndex(const Vec3T& p, int& i, int& j, int& k) restrict(cpu,amp) 
	{
		Vec3T v = p-min_;
		i = (int)((v.x/dx_)-0.5) - ghost_width_;  
		j = (int)((v.y/dy_)-0.5) - ghost_width_;  
		k = (int)((v.z/dz_)-0.5) - ghost_width_; 
	}

	void RightUpIndex(const Vec3T& p, int& i, int& j, int& k) restrict(cpu,amp)
	{
		Vec3T v = p-min_;
		i = (int)((v.x/dx_)+0.5) - ghost_width_; 
		j = (int)((v.y/dy_)+0.5) - ghost_width_; 
		k = (int)((v.z/dz_)+0.5) - ghost_width_; 
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

		glPushMatrix();
		glColor3f(0, 0, 1);
		glTranslatef(grid_center.x, grid_center.y, grid_center.z);
		glScaled(deviation.x, deviation.y, deviation.z);
		glutWireCube(1);
		glPopMatrix();


		glPushMatrix();
		glColor3f(1, 0, 0);
		glTranslatef(grid_center.x, grid_center.y, grid_center.z);
		glScaled(deviation.x-gx_*2.0, deviation.y-gx_*2.0, deviation.z-gx_*2.0);
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

//			if(!IsGhostCell(i,j,k))
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