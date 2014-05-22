
#pragma once 

#include <atomic>
#include <iostream>

template<class TT>
class ARRAY_DYNAMIC
{
private:

	TT* arr_;

	int num_tot_;

	std::atomic<int>  num_cur_;	
	std::atomic<bool> mem_lock_;

public:

	ARRAY_DYNAMIC() : num_cur_(0), num_tot_(0), arr_(0), mem_lock_(false)
	{}

	~ARRAY_DYNAMIC()
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
		while(mem_lock_){};
		int ix = atomic_fetch_add(&num_cur_, count);		
		Rearray();
		return ix;
	}

	int Pop(const int count=1)
	{
		while(mem_lock_){};
		return atomic_fetch_sub(&num_cur_, count)-count;
	}

	void Rearray()
	{ 
		if(num_cur_ >= num_tot_) 
		{
			if(atomic_exchange(&mem_lock_, true) == false)
			{
				int inc = (int)num_tot_/4 + ((int)num_cur_-(int)num_tot_);

				TT* temp = new TT[inc+num_tot_];

				for(int i=0; i<(int)num_tot_; i++) temp[i] = arr_[i];

				delete[] arr_; arr_ = temp;

				num_tot_ = inc+num_tot_;

				mem_lock_.store(false);
			}
		}
	}

	TT& operator() (const int ix) const
	{
		while(mem_lock_){};
		return arr_[ix];
	}
};