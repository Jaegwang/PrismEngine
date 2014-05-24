#pragma once
#include "FIELD.h"
#include "SPARSE_MATRIX.h"
#include "FIELD_UNIFORM.h"

enum BOUNDARY_CELL
{
	BND_FRAG= -2,
	BND_WALL= -1,
	BND_FULL=  0,
	BND_NULL=  1,
};

class PROJECTION_METHOD
{
public:

	SPARSE_MATRIX matrix_a_;

	ARRAY_VECTOR<FLT> vector_x_;
	ARRAY_VECTOR<FLT> vector_b_;

	ARRAY_VECTOR<FLT> vector_r_;
	ARRAY_VECTOR<FLT> vector_p_;

	ARRAY_VECTOR<FLT> vector_ap_;
	ARRAY_VECTOR<FLT> vector_temp_;

public:

	PROJECTION_METHOD()
	{

	}

public:

	void InitializeLinearSystem(FIELD<int>* bnd, FIELD<int>* ptr, FIELD<FLT>* div)
	{
		int full_ix(0);
		int start_full_ix(0);

		int nnz(0);

		GRID grid = bnd->Grid();

		const int i_res = grid.i_res_;
		const int ij_res = grid.ij_res_;

		for(int i=0; i<grid.ijk_res_; i++)
		{
			if(bnd->Get(i) == 0)
			{
				full_ix++;

				nnz++;
				ptr->Set(i, start_full_ix++);

				if(bnd->Get(i+1)      >= 0){ nnz++; ptr->Set(i+1     , start_full_ix++);}
				if(bnd->Get(i-1)      >= 0){ nnz++; ptr->Set(i-1     , start_full_ix++);}
				if(bnd->Get(i+i_res)  >= 0){ nnz++; ptr->Set(i+i_res , start_full_ix++);}
				if(bnd->Get(i-i_res)  >= 0){ nnz++; ptr->Set(i-i_res , start_full_ix++);}
				if(bnd->Get(i+ij_res) >= 0){ nnz++;	ptr->Set(i+ij_res, start_full_ix++);}
				if(bnd->Get(i-ij_res) >= 0){ nnz++;	ptr->Set(i-ij_res, start_full_ix++);}
			}
		}

		// Initialize Matrix and Vectors
		{
			matrix_a_.Initialize(nnz, full_ix);

			vector_b_.Initialize(full_ix, true);
			vector_x_.Initialize(full_ix, true);

			vector_r_.Initialize(full_ix, true);
			vector_p_.Initialize(full_ix, true);

			vector_ap_.Initialize(full_ix, true);
			vector_temp_.Initialize(full_ix, true);
		}

		ARRAY_VECTOR<FLT> val;
		ARRAY_VECTOR<int> col;

		val.Initialize(10);
		col.Initialize(10);

		for(int p=0; p<grid.ijk_res_; p++)
		{
			if(bnd->Get(p) != 0) continue;

			const int bc_ix = bnd->Get(p);
			int p_ijk = 0;

			// TODO : assign matrix A
			if(bnd->Get(p+1) >= 0)
			{
				p_ijk++;
				val.Push(-1.0f); col.Push(ptr->Get(p+1));
			}
			else p_ijk +=2;			

			if(bnd->Get(p-1) >= 0)
			{
				p_ijk++;
				val.Push(-1.0f); col.Push(ptr->Get(p-1));
			}
			else p_ijk += 2;

			if(bnd->Get(p+i_res) >= 0)
			{
				p_ijk++;
				val.Push(-1.0f); col.Push(ptr->Get(p+i_res));
			}
			else p_ijk += 2;

			if(bnd->Get(p-i_res) >= 0)
			{
				p_ijk++;
				val.Push(-1.0f); col.Push(ptr->Get(p-i_res));
			}
			else p_ijk += 2;

			if(bnd->Get(p+ij_res) >= 0)
			{
				p_ijk++;
				val.Push(-1.0f); col.Push(ptr->Get(p+ij_res));
			}
			else p_ijk += 2;

			if(bnd->Get(p-ij_res) >= 0)
			{
				p_ijk++;
				val.Push(-1.0f); col.Push(ptr->Get(p-ij_res));
			}
			else p_ijk += 2;

			val.Push(p_ijk); col.Push(ptr->Get(p));

			matrix_a_.AddRow(val, col);
						
			vector_b_.Set(bc_ix, div->Get(p));
			vector_x_.Set(bc_ix, (FLT)0);			

			val.Clear();
			col.Clear();
		}
	}

	void ConjugateGradient()
	{
		matrix_a_.Multiply(vector_x_, vector_r_);
		ArrayMinus(vector_b_, vector_r_, vector_r_);

		ArrayEqual(vector_r_, vector_p_);

		FLT rsold = ArrayDot(vector_r_, vector_r_);

		for(int i=0; i<10; i++)
		{
			matrix_a_.Multiply(vector_p_, vector_ap_);

			FLT alpha = rsold / ArrayDot(vector_p_, vector_ap_);

			ArrayMultipy(vector_p_, alpha, vector_temp_);
			ArrayPlus(vector_x_, vector_temp_, vector_x_);

			ArrayMultipy(vector_ap_, alpha, vector_temp_);
			ArrayMinus(vector_r_, vector_temp_, vector_r_);

			FLT rsnew = ArrayDot(vector_r_, vector_r_);

			if(sqrt(rsnew) < 1e-6) break;
			
			FLT beta = rsnew/rsold;
			ArrayMultipy(vector_p_, beta, vector_temp_);
			ArrayPlus(vector_r_, vector_temp_, vector_p_);

			rsold = rsnew;		
		}
	}
	
	void DetermineDivergence(const FIELD<int>* bnd, const FIELD<Vec3>* vel, FIELD<FLT>* div);

	void DeterminePressure(const FIELD<int>* bnd, const FIELD<FLT>* div, FIELD<FLT>* press, FIELD<FLT>* press_temp, const int itr=20);

	void DetermineVelocity(const FIELD<int>* bnd, const FIELD<FLT>* press, FIELD<Vec3>* vel);

	void Jacobi(FIELD<int>* bnd, FIELD<Vec3>* vel, FIELD<FLT>* div, FIELD<FLT>* press, FIELD<FLT>* press_temp, const int itr)
	{
		

		DetermineDivergence(bnd, vel, div);


		{
			FIELD_UNIFORM<int>* d = new FIELD_UNIFORM<int>;
			d->Initialize(bnd->Grid());

			for(int i=0; i<bnd->Grid().ijk_res_; i++)
			{
				d->Set(i,-1);
			}

			FIELD<int>* p = d;

			InitializeLinearSystem(bnd, p, div);
			ConjugateGradient();
		}

		DeterminePressure(bnd, div, press, press_temp, itr);

		DetermineVelocity(bnd, press, vel);
	}

	void Diffuse(const FIELD<int>* bnd, FIELD<Vec3>* vel, FIELD<Vec3>* temp, const int itr);
};
