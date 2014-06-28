
#pragma once


template<class TT>
class CACHE_BLOCK
{ 
public:

	TT* array_;

	int size_;
	int i_,j_,k_;

public:

	CACHE_BLOCK(const int size, const TT value)
		: array_(0)
	{
		size_  = size;
		array_ = new TT[size];		
	}

	~CACHE_BLOCK()
	{
		if(array_) { delete[] array_; array_ = 0; };		
	}

	void Fill(const TT value)
	{
		for(int i=0; i<size_; i++) array_[i] = value;
	}
};