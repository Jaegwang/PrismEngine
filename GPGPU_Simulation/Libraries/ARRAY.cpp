
#include "ARRAY.h"
#include "MATH_CORE.h"

template<class TT> 
void ArrayMinus(const ARRAY_VECTOR<TT>& a, const ARRAY_VECTOR<TT>& b, ARRAY_VECTOR<TT>& c)
{
	const int size = c.Size();
	for(int i=0; i<size; i++)
	{
		c.Set(i, a.Get(i)-b.Get(i));	
	}
}

template<class TT>
void ArrayPlus(const ARRAY_VECTOR<TT>& a, const ARRAY_VECTOR<TT>& b, ARRAY_VECTOR<TT>& c)
{
	const int size = c.Size();
	for(int i=0; i<size; i++)
	{
		c.Set(i, a.Get(i)+b.Get(i));	
	}
}

template<class TT>
TT ArrayDot(const ARRAY_VECTOR<TT>& a, const ARRAY_VECTOR<TT>& b)
{
	const int size = b.Size();
	TT acc = TT();

	for(int i=0; i<size; i++)
	{
		acc += a.Get(i)*b.Get(i);
	}
	return acc;
}

template<class TT> 
void ArrayEqual(const ARRAY_VECTOR<TT>& a, ARRAY_VECTOR<TT>& b)
{
	const int size = b.Size();
	for(int i=0; i<size; i++)
	{
		b.Set(i, a.Get(i));
	}
}

template<class TT> 
void ArrayMultipy(const ARRAY_VECTOR<TT>& a, const TT b, ARRAY_VECTOR<TT>& c)
{
	const int size = c.Size();
	for(int i=0; i<size; i++)
	{
		c.Set(i, a.Get(i)*b);
	}
}


template void ArrayMinus(const ARRAY_VECTOR<FLT>& a, const ARRAY_VECTOR<FLT>& b, ARRAY_VECTOR<FLT>& c);
template void ArrayPlus(const ARRAY_VECTOR<FLT>& a, const ARRAY_VECTOR<FLT>& b, ARRAY_VECTOR<FLT>& c);
template FLT ArrayDot(const ARRAY_VECTOR<FLT>& a, const ARRAY_VECTOR<FLT>& b);
template void ArrayEqual(const ARRAY_VECTOR<FLT>& a, ARRAY_VECTOR<FLT>& b);
template void ArrayMultipy(const ARRAY_VECTOR<FLT>& a, const FLT b, ARRAY_VECTOR<FLT>& c);