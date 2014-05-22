#pragma once 

template<class TT>
class ARRAY
{
public:

	virtual void Set(const int idx, const TT& data)=0;
	virtual TT   Get(const int idx) const=0;

	virtual int  Size() const=0;

	virtual void Push(TT& data)=0;
	virtual TT   Pop()=0;

	virtual void Rebuild()=0;
};

template class ARRAY<float>;
template class ARRAY<double>;

template<class TT>
class ARRAY_VECTOR : public ARRAY<TT>
{
public:

	virtual void Set(const int idx, const TT& data)
	{
		arr_[idx] = data;	
	}

	virtual TT Get(const int idx) const
	{
		return arr_[idx];
	}

	virtual int Size() const
	{
		return ptr_+1;
	}

	virtual void Push(TT& data)
	{
		arr_[++ptr_] = data;
	}

	virtual TT Pop()
	{
		return arr_[ptr_--];
	}

	virtual void Rebuild()
	{
	
	}

	void Initialize(const int num)
	{
		arr_ = new TT[num];
		ptr_ = -1;
		size_ = num;
	}

	void Finalize()
	{
		if(arr_) delete[] arr_;	
	}

	ARRAY_VECTOR() : arr_(0), ptr_(-1), size_(0)
	{
	
	}

	~ARRAY_VECTOR()
	{
		Finalize();
	}

private:

	TT* arr_;

	int ptr_;

	int size_;
};