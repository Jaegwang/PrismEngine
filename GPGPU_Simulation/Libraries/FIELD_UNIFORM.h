#pragma once 

#include "GRID.h"
#include "FIELD.h"

template<class TT>
class FIELD_UNIFORM : public FIELD<TT>
{
private:

	GRID grid_;
	TT*  arr_;

public:

	FIELD_UNIFORM() : arr_(0)
	{}
	~FIELD_UNIFORM()
	{
		Finalize();
	}

public:

	void Initialize(const GRID& grid_input);
	void Finalize();

	void Set(const int idx, const TT& data);
	void Set(const int i, const int j, const int k, const TT& data);

	TT   Get(const int idx);
	TT   Get(const int i, const int j, const int k);
	TT   Get(const Vec3& p);
};