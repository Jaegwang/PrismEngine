#pragma once
#include "FIELD.h"
#include "SPARSE_MATRIX.h"
#include "FIELD_UNIFORM.h"

enum BOUNDARY_CELL
{
	BND_FRAG= -3,
	BND_WALL= -2,
	BND_NULL= -1,
	BND_FULL=  0
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

	void InitializeLinearSystem(FIELD<int>* bnd, FIELD<FLT>* div)
	{
		int full_ix(0);
		int start_full_ix(0);

		int nnz(0);

		GRID grid = bnd->Grid();

		const int i_res = grid.i_res_;
		const int ij_res = grid.ij_res_;

		for(int i=0; i<grid.ijk_res_; i++)
		{
			if(bnd->Get(i) == BND_FULL)
			{
				nnz++;
				bnd->Set(i, full_ix++);

				if(bnd->Get(i+1)      >= 0) nnz++;
				if(bnd->Get(i-1)      >= 0) nnz++;
				if(bnd->Get(i+i_res)  >= 0) nnz++;
				if(bnd->Get(i-i_res)  >= 0) nnz++;
				if(bnd->Get(i+ij_res) >= 0) nnz++;
				if(bnd->Get(i-ij_res) >= 0) nnz++;
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
			if(bnd->Get(p) < 0) continue;

			const int bc_ix = bnd->Get(p);
			int p_ijk = 0;

			if(bnd->Get(p+1) >= 0)
			{
				p_ijk++;
				val.Push(-1.0f); col.Push(bnd->Get(p+1));
			}
			else if(bnd->Get(p+1) != BND_NULL) p_ijk +=2;			

			if(bnd->Get(p-1) >= 0)
			{
				p_ijk++;
				val.Push(-1.0f); col.Push(bnd->Get(p-1));
			}
			else if(bnd->Get(p-1) != BND_NULL) p_ijk += 2;

			if(bnd->Get(p+i_res) >= 0)
			{
				p_ijk++;
				val.Push(-1.0f); col.Push(bnd->Get(p+i_res));
			}
			else if(bnd->Get(p+i_res) != BND_NULL) p_ijk += 2;

			if(bnd->Get(p-i_res) >= 0)
			{
				p_ijk++;
				val.Push(-1.0f); col.Push(bnd->Get(p-i_res));
			}
			else if(bnd->Get(p-i_res) != BND_NULL) p_ijk += 2;

			if(bnd->Get(p+ij_res) >= 0)
			{
				p_ijk++;
				val.Push(-1.0f); col.Push(bnd->Get(p+ij_res));
			}
			else if(bnd->Get(p+ij_res) != BND_NULL) p_ijk += 2;

			if(bnd->Get(p-ij_res) >= 0)
			{
				p_ijk++;
				val.Push(-1.0f); col.Push(bnd->Get(p-ij_res));
			}
			else if(bnd->Get(p-ij_res) != BND_NULL) p_ijk += 2;

			val.Push(p_ijk); col.Push(bnd->Get(p));

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

		for(int i=0; i<20; i++)
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



	void DeterminePressureField(const FIELD<int>* bnd, FIELD<FLT>* press)
	{
		GRID grid = bnd->Grid();

		for(int i=0; i<grid.ijk_res_; i++)
		{
			int bc = bnd->Get(i);

			if(bc >= 0)
			{
				FLT p = vector_x_.Get(bc);
				press->Set(i, vector_x_.Get(bc));			
			}
			else
			{
				press->Set(i, (FLT)0);
			}
		}			
	}
	
	void DetermineDivergence(const FIELD<int>* bnd, const FIELD<Vec3>* vel, FIELD<FLT>* div);

	void DeterminePressure(const FIELD<int>* bnd, const FIELD<FLT>* div, FIELD<FLT>* press, FIELD<FLT>* press_temp, const int itr=20);

	void DetermineVelocity(const FIELD<int>* bnd, const FIELD<FLT>* press, FIELD<Vec3>* vel);

	void Jacobi(FIELD<int>* bnd, FIELD<Vec3>* vel, FIELD<FLT>* div, FIELD<FLT>* press, FIELD<FLT>* press_temp, const int itr);

	void Diffuse(const FIELD<int>* bnd, FIELD<Vec3>* vel, FIELD<Vec3>* temp, const int itr);
};

	static void InitializeProjectionLinearSystem(FIELD<int>* bnd, FIELD<FLT>* div, SPARSE_MATRIX& matrix_a, ARRAY_VECTOR<FLT>& vector_x, ARRAY_VECTOR<FLT>& vector_b)
	{
		int full_ix(0);
		int nnz(0);

		GRID grid = bnd->Grid();

		FLT h  = 1/(grid.dx_*grid.dx_);
		FLT h2 = h*(FLT)2; 

		const int i_res = grid.i_res_;
		const int ij_res = grid.ij_res_;

		for(int i=0; i<grid.ijk_res_; i++)
		{
			if(bnd->Get(i) == BND_FULL)
			{
				nnz++;
				bnd->Set(i, full_ix++);

				if(bnd->Get(i+1)      >= 0) nnz++;
				if(bnd->Get(i-1)      >= 0) nnz++;
				if(bnd->Get(i+i_res)  >= 0) nnz++;
				if(bnd->Get(i-i_res)  >= 0) nnz++;
				if(bnd->Get(i+ij_res) >= 0) nnz++;
				if(bnd->Get(i-ij_res) >= 0) nnz++;
			}
		}

		// Initialize Matrix and Vectors
		{
			matrix_a.Initialize(nnz, full_ix);

			vector_b.Initialize(full_ix, true);
			vector_x.Initialize(full_ix, true);
		}

		ARRAY_VECTOR<FLT> val; val.Initialize(16);
		ARRAY_VECTOR<int> col; col.Initialize(16);

		for(int p=0; p<grid.ijk_res_; p++)
		{
			if(bnd->Get(p) < 0) continue;

			const int bc_ix = bnd->Get(p);
			FLT p_ijk = 0;

			if(bnd->Get(p+1) >= 0)
			{
				p_ijk -= h;
				val.Push(h); col.Push(bnd->Get(p+1));
			}
			else if(bnd->Get(p+1) == BND_NULL) p_ijk -= h;			

			if(bnd->Get(p-1) >= 0)
			{
				p_ijk -= h;
				val.Push(h); col.Push(bnd->Get(p-1));
			}
			else if(bnd->Get(p-1) == BND_NULL) p_ijk -= h;

			if(bnd->Get(p+i_res) >= 0)
			{
				p_ijk -= h;
				val.Push(h); col.Push(bnd->Get(p+i_res));
			}
			else if(bnd->Get(p+i_res) == BND_NULL) p_ijk -= h;			

			if(bnd->Get(p-i_res) >= 0)
			{
				p_ijk -= h;
				val.Push(h); col.Push(bnd->Get(p-i_res));
			}
			else if(bnd->Get(p-i_res) == BND_NULL) p_ijk -= h;			

			if(bnd->Get(p+ij_res) >= 0)
			{
				p_ijk -= h;
				val.Push(h); col.Push(bnd->Get(p+ij_res));
			}
			else if(bnd->Get(p+ij_res) == BND_NULL) p_ijk -= h;			

			if(bnd->Get(p-ij_res) >= 0)
			{
				p_ijk -= h;
				val.Push(h); col.Push(bnd->Get(p-ij_res));
			}
			else if(bnd->Get(p-ij_res) == BND_NULL) p_ijk -= h;			

			val.Push(p_ijk); col.Push(bnd->Get(p));

			matrix_a.AddRow(val, col);

			vector_b.Set(bc_ix, div->Get(p));
			vector_x.Set(bc_ix, (FLT)0);			

			val.Clear();
			col.Clear();
		}
	}

	static void ConjugateGradient_(SPARSE_MATRIX& matrix_a, ARRAY_VECTOR<FLT>& vector_b, ARRAY_VECTOR<FLT>& vector_x)
	{
		int size = vector_b.Size();
		ARRAY_VECTOR<FLT> vector_r, vector_p, vector_ap, vector_t;

		vector_r.Initialize(size, true);
		vector_p.Initialize(size, true);
		vector_t.Initialize(size, true);
		vector_ap.Initialize(size, true);

		matrix_a.Multiply(vector_x, vector_r);
		ArrayMinus(vector_b, vector_r, vector_r);

		ArrayEqual(vector_r, vector_p);

		FLT rsold = ArrayDot(vector_r, vector_r);
		if(sqrt(rsold) < (FLT)1E-6) return;

		for(int i=0; i<300; i++)
		{
			matrix_a.Multiply(vector_p, vector_ap);

			FLT alpha = rsold / ArrayDot(vector_p, vector_ap);

			ArrayMultipy(vector_p, alpha, vector_t);
			ArrayPlus(vector_x, vector_t, vector_x);

			ArrayMultipy(vector_ap, alpha, vector_t);
			ArrayMinus(vector_r, vector_t, vector_r);

			FLT rsnew = ArrayDot(vector_r, vector_r);

			if(sqrt(rsnew) < (FLT)1E-6) break;
			
			FLT beta = rsnew/rsold;
			ArrayMultipy(vector_p, beta, vector_t);
			ArrayPlus(vector_r, vector_t, vector_p);

			rsold = rsnew;		
		}
	}