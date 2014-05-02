
#pragma once 

#include <atomic>
#include <iostream>
#include <vector>

template<class TT>
class VECTOR_ATOMIC
{
private:

	std::vector<TT>  arr_;
	std::atomic<int> tail_ptr_;

public:

	VECTOR_ATOMIC() : tail_ptr_(-1)
	{}

	~VECTOR_ATOMIC()
	{}

	void Initialize(const int num)
	{
		arr_.clear();
		arr_.reserve(num);

		for(int i=0; i<num; i++) arr_.push_back(TT());

		tail_ptr_ = -1;
	}

	void Reset()
	{
		tail_ptr_ = -1;
	}

	const int Size()
	{
		return (int)tail_ptr_+1;
	}

	int Push(const int count)
	{	
		int ix = atomic_fetch_add(&tail_ptr_, count)+count;
		Rearray();

		return ix;
	}

	void Rearray()
	{ // TODO : run on one thread
		if((int)tail_ptr_+1 < (int)arr_.size()) return;

		int increase = (int)arr_.size()/4 + (tail_ptr_-(int)arr_.size()+1);

		for(int i=0; i<increase; i++) arr_.push_back(TT());
	}

	TT& operator() (const int ix) const
	{
		TT* value = (TT*)arr_.data();
		return value[ix];
	}
};