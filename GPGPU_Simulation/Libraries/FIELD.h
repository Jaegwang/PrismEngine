#pragma once 

#include "GRID.h"

template<class TT>
class FIELD
{
protected:

	GRID grid_;

public:

	GRID Grid() const { return grid_; }

	virtual void Set(const int idx, const TT& data)=0;
	virtual void Set(const int i, const int j, const int k, const TT& data)=0;

	virtual TT   Get(const int idx) const=0;
	virtual TT   Get(const int i, const int j, const int k) const=0;
	virtual TT   Get(const Vec3& p) const=0;
};


template class FIELD<Vec3>;
template class FIELD<float>;
template class FIELD<double>;