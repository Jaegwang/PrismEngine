#pragma once 

template<class TT>
class ARRAY
{
public:

	virtual void Set(const int idx, const TT& data)=0;
	virtual TT   Get(const int idx) const=0;

	virtual int  Size() const=0;

	virtual void Push(TT data)=0;
	virtual TT   Pop()=0;

	virtual void Clear()=0;
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

	virtual void Push(TT data)
	{
		arr_[++ptr_] = data;
	}

	virtual TT Pop()
	{
		return arr_[ptr_--];
	}

	virtual void Clear()
	{
		ptr_ = -1;
	}

	virtual void Rebuild()
	{
	
	}

	void Initialize(const int num, const bool rev=false)
	{
		Finalize();

		arr_   = new TT[num];
		total_ = num;

		for(int i=0; i<num; i++)
		{
			arr_[i] = TT();
		}

		if(rev == true) ptr_ = num-1;
		else ptr_ = -1;
	}

	void Finalize()
	{
		if(arr_) delete[] arr_;	
	}

	ARRAY_VECTOR() : arr_(0), ptr_(-1), total_(0)
	{
	
	}

	~ARRAY_VECTOR()
	{
		Finalize();
	}

private:

	TT* arr_;

	int ptr_;

	int total_;
};

template<class TT> void ArrayMinus(const ARRAY_VECTOR<TT>& a, const ARRAY_VECTOR<TT>& b, ARRAY_VECTOR<TT>& c);

template<class TT> void ArrayPlus(const ARRAY_VECTOR<TT>& a, const ARRAY_VECTOR<TT>& b, ARRAY_VECTOR<TT>& c);

template<class TT> TT ArrayDot(const ARRAY_VECTOR<TT>& a, const ARRAY_VECTOR<TT>& b);

template<class TT> void ArrayEqual(const ARRAY_VECTOR<TT>& a, ARRAY_VECTOR<TT>& b);

template<class TT> void ArrayMultipy(const ARRAY_VECTOR<TT>& a, const TT b, ARRAY_VECTOR<TT>& c);