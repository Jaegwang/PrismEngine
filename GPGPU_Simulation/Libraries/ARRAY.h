#pragma once 

template<class TT>
class ARRAY
{
public:

	virtual void Set(const int idx, const TT& data)=0;
	virtual TT   Get(const int idx) const=0;

	virtual int  Size() const=0;

	virtual void Rebuild()=0;
};

template class ARRAY<float>;
template class ARRAY<double>;

template<class TT>
class ARRAY_VECTOR
{
public:

	virtual void Set(const int idx, const TT& data)
	{
		arr_[idx] = data;	
	}

	virtual TT Get(const int idx)
	{
		return arr_[idx];
	}

	virtual int Size()
	{
		return size_;
	}

	virtual void Rebuild()
	{
	
	}

	void Initialize(const int num)
	{
		arr_ = new TT[num];
		size_ = num;
	}

private:

	TT* arr_;

	int size_;
};