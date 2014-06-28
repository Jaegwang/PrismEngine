#pragma once

#include <atomic> 
 
template<class TT>
class STACK_FREE
{
private:

	TT* array_;
	std::atomic<int> ptr_;

	int total_;

public:

	STACK_FREE()
		: array_(0)
		, ptr_(-1)
	{}
	~STACK_FREE()
	{
		Finalize();
	}

	void Initialize(const int total_size)
	{
		array_ = new TT[total_size];	
		total_ = total_size;
	}

	void Finalize()
	{
		if(array_) { delete[] array_; array_ = 0; };		
	}

	void Push(TT data)
	{
		int p = std::atomic_fetch_add(&ptr_, 1);
		array_[p+1] = data;	
	}

	TT Pop()
	{
		int p = std::atomic_fetch_sub(&ptr_, 1);
		return array_[p];
	}	

	bool IsEmpty()
	{
		if(ptr_ < 0) return true;
		else return false;
	}
};