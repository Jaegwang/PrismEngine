
#pragma once 

#include "VECTOR3_T.h"
#include "GL\glut.h"

class PARTICLE_OBJECT
{
public:

	int pts_num_;
	T   pts_mass_;

	Vec3T* position_array_;
	Vec3T* velocity_array_;

public:

	PARTICLE_OBJECT()
	{}
	~PARTICLE_OBJECT()
	{}

	void InitializeCube(const Vec3T& min, const Vec3T& max, const T dx, const T dy, const T dz)
	{
		int i_res = (max.x-min.x)/dx;
		int j_res = (max.y-min.y)/dy;
		int k_res = (max.z-min.z)/dz;

		position_array_ = new Vec3T[i_res*j_res*k_res];
		velocity_array_ = new Vec3T[i_res*j_res*k_res];

		pts_mass_ = (T)1;

		int ix = 0;
		for(int k=0; k<k_res; k++) for(int j=0; j<j_res; j++) for(int i=0; i<i_res; i++)
		{
			position_array_[ix] = min + Vec3T(dx*(T)i, dy*(T)j, dz*(T)k);
			velocity_array_[ix] = Vec3T();		

			ix++;
			pts_num_++;
		}	
	}

	void RenderObject()
	{
		glDisable(GL_LIGHTING);

		glPushMatrix();
		glColor3f(0, 0, 0);
		
		glBegin(GL_POINTS);
		for(int x=0; x<pts_num_; x++)
		{
			const Vec3T& pos = position_array_[x];
			glVertex3f(pos.x, pos.y, pos.z);		
		}
		glEnd();

		glPopMatrix();
		glEnable(GL_LIGHTING);
	}
};