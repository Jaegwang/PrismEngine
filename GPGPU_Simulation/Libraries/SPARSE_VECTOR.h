#pragma once
#include "MATH_CORE.h"

class SPARSE_VECTOR
{
private:

	FLT* val_arr_;
	int  size_;

public:

	SPARSE_VECTOR() :
		val_arr_(0), size_(0)
	{}

	~SPARSE_VECTOR()
	{}

	void Initialize(const int num)
	{
		Finalize();

		val_arr_ = new FLT[num];
		size_ = num;	
	}

	void Finalize()
	{
		if(val_arr_) delete[] val_arr_;		
	}

	int Size() { return size_; }

	FLT& operator() (const int p) const { return val_arr_[p]; };

	static void Add(const SPARSE_VECTOR& a, const SPARSE_VECTOR& b, SPARSE_VECTOR& c)
	{
		#pragma omp parallel for
		for(int p=0; p<c.size_; p++)
		{
			c.val_arr_[p] = a.val_arr_[p] + b.val_arr_[p];		
		}
	}

	static void Sub(const SPARSE_VECTOR& a, const SPARSE_VECTOR& b, SPARSE_VECTOR& c)
	{
		#pragma omp parallel for
		for(int p=0; p<c.size_; p++)
		{
			c.val_arr_[p] = a.val_arr_[p] - b.val_arr_[p];		
		}
	}

	static void Equal(const SPARSE_VECTOR& a, SPARSE_VECTOR& c)
	{
		#pragma omp parallel for
		for(int p=0; p<c.size_; p++)
		{
			c.val_arr_[p] = a.val_arr_[p];
		}	
	}

	static void Mul(const SPARSE_VECTOR& a, const FLT b, SPARSE_VECTOR& c)
	{
		#pragma omp parallel for
		for(int p=0; p<c.size_; p++)
		{
			c.val_arr_[p] = a.val_arr_[p] * b;
		}		
	}

	static FLT Dot(const SPARSE_VECTOR& a, const SPARSE_VECTOR& b)
	{
		FLT d = FLT(0);

		for(int p=0; p<a.size_; p++)
		{
			d += a.val_arr_[p] * b.val_arr_[p];
		}	

		return d;
	}
};


