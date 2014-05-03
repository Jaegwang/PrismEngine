
#pragma once 

#include <atomic>
#include <iostream>

template<class TT>
class ARRAY_ATOMIC
{
private:

	TT* arr_;

	int num_tot_;

	std::atomic<int> num_cur_;	

public:

	ARRAY_ATOMIC() : num_cur_(0), num_tot_(0), arr_(0)
	{}

	~ARRAY_ATOMIC()
	{
		Finalize();
	}

	void Initialize(const int num)
	{
		Finalize();

		arr_ = new TT[num];

		num_tot_ = num;
		num_cur_ = 0;
	}

	void Finalize()
	{
		if(arr_) delete[] arr_;
	}

	void Reset()
	{
		num_cur_ = 0;
	}

	const int Size()
	{
		return (int)num_cur_;
	}

	const int Total()
	{
		return num_tot_;
	}

	int Push(const int count=1)
	{	
		int ix = atomic_fetch_add(&num_cur_, count);
		Rearray();

		return ix;
	}

	int Pop(const int count=1)
	{
		int ix = atomic_fetch_sub(&num_cur_, count)-count;
		return ix;	
	}

	void Rearray()
	{ 
		#pragma omp critical
		{
			if((int)num_cur_ >= (int)num_tot_) 
			{
				int inc = (int)num_tot_/4 + (num_cur_-(int)num_tot_);

				TT* temp = new TT[inc+num_tot_];

				for(int i=0; i<num_tot_; i++) temp[i] = arr_[i];

				delete[] arr_; arr_ = temp;

				num_tot_ = inc+num_tot_;
			}
		}
	}

	TT& operator() (const int ix) const
	{
		return arr_[ix];
	}
};