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

	SPARSE_VECTOR vector_x_;
	SPARSE_VECTOR vector_b_;

public:

	PROJECTION_METHOD()
	{

	}

public:

	void DeterminePressureField(const FIELD<int>* bnd, FIELD<FLT>* press)
	{
		GRID grid = bnd->Grid();

		for(int i=0; i<grid.ijk_res_; i++)
		{
			int bc = bnd->Get(i);

			if(bc >= 0)
			{
				FLT p = vector_x_(bc);
				press->Set(i, vector_x_(bc));			
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

static void InitializeProjectionLinearSystem(FIELD<int>* bnd, FIELD<FLT>* div, SPARSE_MATRIX& matrix_a, SPARSE_VECTOR& vector_x, SPARSE_VECTOR& vector_b)
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

		vector_b.Initialize(full_ix);
		vector_x.Initialize(full_ix);
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

		vector_b(bc_ix) = div->Get(p);
		vector_x(bc_ix) = (FLT)0;			

		val.Clear();
		col.Clear();
	}
}

static void ConjugateGradient_(SPARSE_MATRIX& matrix_a, SPARSE_VECTOR& vector_b, SPARSE_VECTOR& vector_x)
{
	int size = vector_b.Size();
	SPARSE_VECTOR vector_r, vector_p, vector_ap, vector_t;

	vector_r.Initialize(size);  vector_p.Initialize(size);
	vector_ap.Initialize(size); vector_t.Initialize(size);

	SPARSE_MATRIX::Mul(matrix_a, vector_x, vector_r);
	SPARSE_VECTOR::Sub(vector_b, vector_r, vector_r);

	SPARSE_VECTOR::Equal(vector_r, vector_p);

	FLT rsold = SPARSE_VECTOR::Dot(vector_r, vector_r);
	if(sqrt(rsold) < (FLT)1E-6) return;

	for(int i=0; i<300; i++)
	{
		SPARSE_MATRIX::Mul(matrix_a, vector_p, vector_ap);

		FLT alpha = rsold / SPARSE_VECTOR::Dot(vector_p, vector_ap);

		SPARSE_VECTOR::Mul(vector_p, alpha, vector_t);
		SPARSE_VECTOR::Add(vector_x, vector_t, vector_x);

		SPARSE_VECTOR::Mul(vector_ap, alpha, vector_t);
		SPARSE_VECTOR::Sub(vector_r, vector_t, vector_r);

		FLT rsnew = SPARSE_VECTOR::Dot(vector_r, vector_r);

		if(sqrt(rsnew) < (FLT)1E-6) break;
			
		FLT beta = rsnew/rsold;
		SPARSE_VECTOR::Mul(vector_p, beta, vector_t);
		SPARSE_VECTOR::Add(vector_r, vector_t, vector_p);

		rsold = rsnew;		
	}
}