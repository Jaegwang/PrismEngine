#pragma once
#include "MATH_CORE.h"
#include "ARRAY.h"

class SPARSE_MATRIX
{ //http://netlib.org/linalg/html_templates/node91.html
private:

	FLT* val_arr_;

	int* col_ind_arr_;

	int* row_ptr_arr_;

	int  val_ptr_;

	int  row_ptr_;

	int  total_;

public:

	SPARSE_MATRIX() : 
		val_arr_(0), 
		col_ind_arr_(0), 
		row_ptr_arr_(0),
		val_ptr_(0),
		row_ptr_(0),
		total_(0)
	{}

	void Initialize(const int nnz_num, const int row_num)
	{
		Finalize();

		total_ = nnz_num;

		val_arr_     = new FLT[nnz_num];
		col_ind_arr_ = new int[nnz_num];
		row_ptr_arr_ = new int[row_num];

		row_ptr_ = -1;
		val_ptr_ = -1;
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

	void Multiply(const ARRAY<FLT>& x, ARRAY<FLT>& b)
	{
		const int size = x.Size();

		for(int p=0; p<size; p++)
		{
			int start_ptr = row_ptr_arr_[p];
			int end_ptr;

			if(p==size-1) end_ptr = total_-1;
			else end_ptr = row_ptr_arr_[p+1];

			FLT v = (FLT)0;

			for(int ix=start_ptr; ix<=end_ptr; ix++)
			{
				int c = col_ind_arr_[ix];

				v += val_arr_[ix] * x.Get(c);
			}

			b.Set(p, v);
		}	
	}
};