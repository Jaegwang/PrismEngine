#pragma once
#include "MATH_CORE.h"
#include "SPARSE_VECTOR.h"

class SPARSE_MATRIX
{ //http://netlib.org/linalg/html_templates/node91.html
private:

	FLT* val_arr_;

	int* col_ind_arr_;

	int* row_ptr_arr_;

	int  val_ptr_;

	int  row_ptr_;

	int  nnz_num_;
	int  row_num_;

public:

	SPARSE_MATRIX() : 
		val_arr_(0), 
		col_ind_arr_(0), 
		row_ptr_arr_(0),
		val_ptr_(0),
		row_ptr_(0)
	{}

	void Initialize(const int nnz_num, const int row_num)
	{
		Finalize();

		val_arr_     = new FLT[nnz_num];
		col_ind_arr_ = new int[nnz_num];

		row_ptr_arr_ = new int[row_num];

		row_ptr_ = -1;
		val_ptr_ = -1;

		nnz_num_ = nnz_num;
		row_num_ = row_num;
	}

	void Finalize()
	{
		if(val_arr_)     delete[] val_arr_;
		if(col_ind_arr_) delete[] col_ind_arr_;
		if(row_ptr_arr_) delete[] row_ptr_arr_;
	}

	void Clear()
	{
		row_ptr_ = -1;
		val_ptr_ = -1;
	}

	void AddRow(const ARRAY<FLT>& val_arr, const ARRAY<int>& col_arr)
	{
		int size = val_arr.Size();

		row_ptr_arr_[++row_ptr_] = val_ptr_+1;

		for(int p=0; p<size; p++)
		{
			int ix = ++val_ptr_;

			val_arr_[ix] = val_arr.Get(p);
			col_ind_arr_[ix] = col_arr.Get(p);
		}
	}

	static void Mul(const SPARSE_MATRIX matrix_a, const SPARSE_VECTOR vector_b, SPARSE_VECTOR vector_c)
	{
		const int size = vector_c.Size();

		for(int p=0; p<size; p++)
		{
			int start_ptr = matrix_a.row_ptr_arr_[p];
			int end_ptr;

			if(p==size-1) end_ptr = matrix_a.val_ptr_;
			else end_ptr = matrix_a.row_ptr_arr_[p+1]-1;

			FLT v = (FLT)0;

			for(int ix=start_ptr; ix<=end_ptr; ix++)
			{
				int c = matrix_a.col_ind_arr_[ix];

				v += matrix_a.val_arr_[ix] * vector_b(c);
			}

			vector_c(p) = v;
		}	
	}
};