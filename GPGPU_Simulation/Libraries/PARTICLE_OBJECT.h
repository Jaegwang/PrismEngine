
#pragma once 

#include "MATH_CORE.h"
#include "GL\glut.h"

class PARTICLE_OBJECT
{
public:

	int pts_num_;
	FLT   pts_mass_;

	Vec3* position_array_;
	Vec3* velocity_array_;

public:

	PARTICLE_OBJECT()
	{}
	~PARTICLE_OBJECT()
	{}

	void InitializeCube(const Vec3& min, const Vec3& max, const FLT dx, const FLT dy, const FLT dz)
	{
		int i_res = (max.x-min.x)/dx;
		int j_res = (max.y-min.y)/dy;
		int k_res = (max.z-min.z)/dz;

		position_array_ = new Vec3[i_res*j_res*k_res];
		velocity_array_ = new Vec3[i_res*j_res*k_res];

		pts_mass_ = (FLT)1;

		int ix = 0;
		for(int k=0; k<k_res; k++) for(int j=0; j<j_res; j++) for(int i=0; i<i_res; i++)
		{
			position_array_[ix] = min + Vec3(dx*(FLT)i, dy*(FLT)j, dz*(FLT)k);
			velocity_array_[ix] = Vec3();		

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
			const Vec3& pos = position_array_[x];
			glVertex3f(pos.x, pos.y, pos.z);		
		}
		glEnd();

		glPopMatrix();
		glEnable(GL_LIGHTING);
	}
};