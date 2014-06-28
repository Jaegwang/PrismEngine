
#pragma once 

#include "MATH_CORE.h"
#include "GL\glut.h"

class PARTICLE_OBJECT
{
public:

	int  pts_num_;
	TS   pts_mass_;

	TV3* position_array_;
	TV3* velocity_array_;

public:

	PARTICLE_OBJECT()
	{}
	~PARTICLE_OBJECT()
	{}

	void InitializeCube(const TV3& min, const TV3& max, const TS dx, const TS dy, const TS dz)
	{
		int i_res = (max.x-min.x)/dx;
		int j_res = (max.y-min.y)/dy;
		int k_res = (max.z-min.z)/dz;

		position_array_ = new TV3[i_res*j_res*k_res];
		velocity_array_ = new TV3[i_res*j_res*k_res];

		pts_mass_ = (TS)1;

		int ix = 0;
		for(int k=0; k<k_res; k++) for(int j=0; j<j_res; j++) for(int i=0; i<i_res; i++)
		{
			position_array_[ix] = min + TV3(dx*(TS)i, dy*(TS)j, dz*(TS)k);
			velocity_array_[ix] = TV3();		

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
			const TV3& pos = position_array_[x];
			glVertex3f(pos.x, pos.y, pos.z);		
		}
		glEnd();

		glPopMatrix();
		glEnable(GL_LIGHTING);
	}
};