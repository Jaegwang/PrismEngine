
#pragma once

#include <atomic>
#include <iostream>
#include <stack>

template<class TT>
class FACTORY
{
private:

	TT** arr_;

	int num_tot_;

	std::atomic<int> num_cur_;	

public:

	FACTORY() : num_cur_(0), num_tot_(0), arr_(0)
	{}
	~FACTORY()
	{
		Finalize();
	}

	void Initialize(const int num)
	{
		Finalize();

		arr_ = new TT*[num];

		num_tot_ = num;
		num_cur_ = num;

		for(int i=0; i<num; i++) arr_[i] = new TT;
	}

	void Finalize()
	{	
		for(int i=0; i<num_tot_; i++)
		{
			delete arr_[i];
		}
		if(arr_) delete[] arr_;
	}

	void Rearray()
	{ // TODO : run on one thread
		if((int)num_cur_ > 0) return;

		#pragma omp single
		{
			int ppp = (int)num_cur_;

			int inc = num_tot_/4;

			TT** temp = new TT*[inc+num_tot_];

			for(int i=0; i<inc; i++)
				temp[i] = new TT;

			for(int i=inc; i<inc+num_tot_; i++)
			{
				temp[i] = arr_[i-inc];

	//			(*temp[i]) = TT();

				int aaa = 0;
			}

			delete[] arr_;
			arr_ = temp;

			num_cur_ = inc;
			num_tot_ = inc+num_tot_;
		}
	}

	TT* Pop()
	{
		Rearray();
		int ix = atomic_fetch_sub(&num_cur_, 1)-1;

		return arr_[ix];		
	}

	void Push(TT* data)
	{
		int ix = atomic_fetch_add(&num_cur_, 1);
		arr_[ix] = data;
	}
};


template<class TT>
class INSTANCE
{
private:

	std::stack<TT*> free_stack_;
	std::stack<TT*> used_stack_;

public:

	void Push(TT* data)
	{
	
	}

}