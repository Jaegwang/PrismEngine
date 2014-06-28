
#include "PROJECTION_METHOD.h"
#include "GENERAL_MACRO.h"

void PROJECTION_METHOD::DetermineDivergence(const FIELD<int>* bnd, const FIELD<TV3>* vel, FIELD<TS>* div)
{
	GRID grid = vel->Grid();

	const TS coef = (TS)1.0;

	FOR_EACH_PARALLER(k, 1, grid.k_res_-2)
	{
		for(int j=1; j<grid.j_res_-1; j++) for(int i=1; i<grid.i_res_-1; i++)
		{
			if(bnd->Get(i,j,k) == BND_FULL)
			{
				TS d = (TS)0;

				TV3 v = vel->Get(i,j,k);
				TV3 a;

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
			
				d *= grid.one_over_dx_*(TS)0.5;

				div->Set(i,j,k,d);
				continue;
			}
			else
			{			
				div->Set(i,j,k,(TS)0);
				continue;
			}
		}	
	}
}

void PROJECTION_METHOD::DeterminePressure(const FIELD<int>* bnd, const FIELD<TS>* div, FIELD<TS>* press, FIELD<TS>* press_temp, const int itr)
{
	GRID grid = press->Grid();

	FOR_EACH_PARALLER(p, 0, grid.ijk_res_-1)
	{
		press->Set(p,(TS)0);
		press_temp->Set(p,(TS)0);
	}

	for(int t=0; t<itr; t++)
	{
		FOR_EACH_PARALLER(k, 1, grid.k_res_-2)
		{
			for(int j=1; j<grid.j_res_-1; j++) for(int i=1; i<grid.i_res_-1; i++)
			{
				if(bnd->Get(i,j,k) >= 0)
				{
					TS cp = press_temp->Get(i,j,k);
					TS p  = (TS)0;

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
					p /= (TS)6;

					press->Set(i,j,k,p);
				}
			}
		
		}

		FIELD<TS>* temp;
		SWAP(press, press_temp, temp);
	}
}

void PROJECTION_METHOD::DetermineVelocity(const FIELD<int>* bnd, const FIELD<TS>* press, FIELD<TV3>* vel)
{
	GRID grid = press->Grid();

	FOR_EACH_PARALLER(k, 1, grid.k_res_-2)
	{
		for(int j=1; j<grid.j_res_-1; j++) for(int i=1; i<grid.i_res_-1; i++)
		{
			if(bnd->Get(i,j,k) >= 0)
			{
				TS  cp = press->Get(i,j,k);
				TV3 p  = TV3();
				TV3 v  = vel->Get(i,j,k); 

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

				p *= grid.one_over_dx_*(TS)0.5;

				vel->Set(i,j,k,v-p);			
			}
		}
	}
}

void PROJECTION_METHOD::Jacobi(FIELD<int>* bnd, FIELD<TV3>* vel, FIELD<TS>* div, FIELD<TS>* press, FIELD<TS>* press_temp, const TS dt, const int itr)
{
	DetermineDivergence(bnd, vel, div);

	{
		InitializeProjectionLinearSystem(bnd, div, matrix_a_, vector_x_, vector_b_);
		ConjugateGradient_(matrix_a_, vector_b_, vector_x_);

		DeterminePressureField(bnd, press, dt);
	}

//	DeterminePressure(bnd, div, press, press_temp, itr);

	DetermineVelocity(bnd, press, vel);
}


void PROJECTION_METHOD::Diffuse(const FIELD<int>* bnd, FIELD<TV3>* vel, FIELD<TV3>* temp, const int itr)
{


}

