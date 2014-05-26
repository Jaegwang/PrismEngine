
#include "PROJECTION_METHOD.h"

void PROJECTION_METHOD::DetermineDivergence(const FIELD<int>* bnd, const FIELD<Vec3>* vel, FIELD<FLT>* div)
{
	GRID grid = vel->Grid();

	const FLT coef = (FLT)1.0;

	#pragma omp parallel for
	for(int k=1; k<grid.k_res_-1; k++)
	for(int j=1; j<grid.j_res_-1; j++) 
	for(int i=1; i<grid.i_res_-1; i++)
	{
		if(bnd->Get(i,j,k) == BND_FULL)
		{
			FLT d = (FLT)0;

			Vec3 v = vel->Get(i,j,k);
			Vec3 a;

			a.x = ABS(v.x);
			a.y = ABS(v.y);
			a.z = ABS(v.z);

			a *= coef;

			int bc;

			bc = bnd->Get(i+1,j,k);
			if(bc >= 0) d += vel->Get(i+1,j,k).x;
			else if(bc == BND_NULL) d += v.x;
			else d -= a.x;

			bc = bnd->Get(i-1,j,k);
			if(bc >= 0) d -= vel->Get(i-1,j,k).x;
			else if(bc == BND_NULL) d -= v.x;
			else d -= a.x;

			bc = bnd->Get(i,j+1,k);
			if(bc >= 0) d += vel->Get(i,j+1,k).y;
			else if(bc == BND_NULL) d += v.y;
			else d -= a.y;

			bc = bnd->Get(i,j-1,k);
			if(bc >= 0) d -= vel->Get(i,j-1,k).y;
			else if(bc == BND_NULL) d -= v.y;
			else d -= a.y;


			bc = bnd->Get(i,j,k+1);
			if(bc >= 0) d += vel->Get(i,j,k+1).z;
			else if(bc == BND_NULL) d += v.z;
			else d -= a.z;

			bc = bnd->Get(i,j,k-1);
			if(bc >= 0) d -= vel->Get(i,j,k-1).z;
			else if(bc == BND_NULL) d -= v.z;
			else d -= a.z;
			
			d *= grid.one_over_dx_*(FLT)0.5;

			div->Set(i,j,k,d);
			continue;
		}
		else
		{			
			div->Set(i,j,k,(FLT)0);
			continue;
		}
	}	
}

void PROJECTION_METHOD::DeterminePressure(const FIELD<int>* bnd, const FIELD<FLT>* div, FIELD<FLT>* press, FIELD<FLT>* press_temp, const int itr)
{
	GRID grid = press->Grid();

	#pragma omp parallel for
	for(int p=0; p<grid.ijk_res_; p++)
	{
		press->Set(p,(FLT)0);
		press_temp->Set(p,(FLT)0);
	}

	for(int t=0; t<itr; t++)
	{
		#pragma omp parallel for
		for(int k=1; k<grid.k_res_-1; k++)
		for(int j=1; j<grid.j_res_-1; j++) 
		for(int i=1; i<grid.i_res_-1; i++)
		{
			if(bnd->Get(i,j,k) >= 0)
			{
				FLT cp = press_temp->Get(i,j,k);
				FLT p  = (FLT)0;

				int bc;

				bc = bnd->Get(i+1,j,k);
				if(bc >= 0) p += press_temp->Get(i+1,j,k);
				else if(bc != BND_NULL) p += cp;

				bc = bnd->Get(i-1,j,k);
				if(bc >= 0) p += press_temp->Get(i-1,j,k);
				else if(bc != BND_NULL) p += cp;


				bc = bnd->Get(i,j+1,k);
				if(bc >= 0) p += press_temp->Get(i,j+1,k);
				else if(bc != BND_NULL) p += cp;

				bc = bnd->Get(i,j-1,k);
				if(bc >= 0) p += press_temp->Get(i,j-1,k);
				else if(bc != BND_NULL) p += cp;


				bc = bnd->Get(i,j,k+1);
				if(bc >= 0) p += press_temp->Get(i,j,k+1);
				else if(bc != BND_NULL) p += cp;

				bc = bnd->Get(i,j,k-1);
				if(bc >= 0) p += press_temp->Get(i,j,k-1);
				else if(bc != BND_NULL) p += cp;


				p -= div->Get(i,j,k)*grid.dx_*grid.dx_;
				p /= (FLT)6;

				press->Set(i,j,k,p);
			}
		}

		FIELD<FLT>* temp;
		SWAP(press, press_temp, temp);
	}
}

void PROJECTION_METHOD::DetermineVelocity(const FIELD<int>* bnd, const FIELD<FLT>* press, FIELD<Vec3>* vel)
{
	GRID grid = press->Grid();

	#pragma omp parallel for
	for(int k=1; k<grid.k_res_-1; k++)
	for(int j=1; j<grid.j_res_-1; j++) 
	for(int i=1; i<grid.i_res_-1; i++)
	{
		if(bnd->Get(i,j,k) >= 0)
		{
			FLT  cp = press->Get(i,j,k);
			Vec3 p  = Vec3();
			Vec3 v  = vel->Get(i,j,k); 

			int bc;

			bc = bnd->Get(i+1,j,k);
			if(bc >= 0) p.x += press->Get(i+1,j,k); 
			else if(bc != BND_NULL) p.x += cp;

			bc = bnd->Get(i-1,j,k);
			if(bc >= 0) p.x -= press->Get(i-1,j,k);
			else if(bc != BND_NULL) p.x -= cp;

			bc = bnd->Get(i,j+1,k);
			if(bc >= 0) p.y += press->Get(i,j+1,k); 
			else if(bc != BND_NULL) p.y += cp;

			bc = bnd->Get(i,j-1,k);
			if(bc >= 0) p.y -= press->Get(i,j-1,k);
			else if(bc != BND_NULL) p.y -= cp;

			bc = bnd->Get(i,j,k+1);
			if(bc >= 0) p.z += press->Get(i,j,k+1); 
			else if(bc != BND_NULL) p.z += cp;

			bc = bnd->Get(i,j,k-1);
			if(bc >= 0) p.z -= press->Get(i,j,k-1);
			else if(bc != BND_NULL) p.z -= cp;

			p *= grid.one_over_dx_*(FLT)0.5;

			vel->Set(i,j,k,v-p);			
		}
	}
}

void PROJECTION_METHOD::Jacobi(FIELD<int>* bnd, FIELD<Vec3>* vel, FIELD<FLT>* div, FIELD<FLT>* press, FIELD<FLT>* press_temp, const int itr)
{
	DetermineDivergence(bnd, vel, div);

	{
		InitializeProjectionLinearSystem(bnd, div, matrix_a_, vector_x_, vector_b_);
		ConjugateGradient_(matrix_a_, vector_b_, vector_x_);

		DeterminePressureField(bnd, press);
	}

//	DeterminePressure(bnd, div, press, press_temp, itr);

	DetermineVelocity(bnd, press, vel);
}


void PROJECTION_METHOD::Diffuse(const FIELD<int>* bnd, FIELD<Vec3>* vel, FIELD<Vec3>* temp, const int itr)
{


}

