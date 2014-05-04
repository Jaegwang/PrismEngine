
#include "FIELD_UNIFORM.h"

template<class TT>
void FIELD_UNIFORM<TT>::Initialize(const GRID& grid_input)
{
	Finalize();

	grid_ = grid_input;
	arr_ = new TT[grid_.ijk_res_];
}

template<class TT>
void FIELD_UNIFORM<TT>::Finalize()
{
	if(arr_) delete[] arr_;
}

template<class TT>
void FIELD_UNIFORM<TT>::Set(const int idx, const TT& data)
{
	arr_[idx] = data;
}

template<class TT>
void FIELD_UNIFORM<TT>::Set(const int i, const int j, const int k, const TT& data)
{
	int idx = grid_.Index3Dto1D(i,j,k);
	arr_[idx] = data;	
}

template<class TT>
TT FIELD_UNIFORM<TT>::Get(const int idx) const
{
	return arr_[idx];
}

template<class TT>
TT FIELD_UNIFORM<TT>::Get(const int i, const int j, const int k) const
{
	int idx = grid_.Index3Dto1D(i,j,k);
	return arr_[idx];
}

template<class TT>
TT FIELD_UNIFORM<TT>::Get(const Vec3& p) const
{
	return grid_.TriLinearInterpolate(p, arr_);
}


template class FIELD_UNIFORM<float>;
template class FIELD_UNIFORM<double>;
template class FIELD_UNIFORM<int>;
template class FIELD_UNIFORM<Vec3>;