
#pragma once

#include <atomic>
#include <iostream>
#include <vector>
#include <stack>

template<class TT>
class FACTORY
{
private:

	std::vector<TT*> arr_;
	std::atomic<int> tail_ptr_;	

public:

	FACTORY() : tail_ptr_(-1)
	{}
	~FACTORY()
	{}

	void Initialize(const int num)
	{
		arr_.clear();
		arr_.reserve(num);

		tail_ptr_ = num-1;

		for(int i=0; i<num; i++) arr_.push_back(new TT());
	}

	void Rearray()
	{ // TODO : run on one thread

	}

	TT* Pop()
	{
		int ix = atomic_fetch_sub(&tail_ptr_, 1);
		Rearray();

		return arr_[ix];		
	}

	void Push(TT* data)
	{
		int ix = atomic_fetch_add(&tail_ptr_, 1)+1;
		Rearray();

		TT** value = arr_.data();
		value[ix] = data;
	}
};