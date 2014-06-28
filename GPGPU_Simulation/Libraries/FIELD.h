#pragma once 

#include "GRID.h"
#include "ARRAY.h"

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
	virtual TT   Get(const TV3& p) const=0;

	virtual void Rebuild()=0;
};


template class FIELD<TV3>;
template class FIELD<float>;
template class FIELD<double>;


template<class TT>
void Rasterize(FIELD<TT>& val_field, FIELD<TS>& weight_field, const ARRAY<TV3>& pos_array, const ARRAY<TT>& val_array)
{
	int size = pos_array.Size();
	int i,j,k;

	int i_start, i_end;
	int j_start, j_end;
	int k_start, k_end;

	GRID grid = val_field.Grid();

	for(int p=0; p<size; p++)
	{
		TV3 pos = pos_array.Get(p);
		TT   val = val_array.Get(p);

		grid.CellCenterIndex(pos,i,j,k);		

		i_start = MAX(0, i-1);
		j_start = MAX(0, j-1);
		k_start = MAX(0, k-1);

		i_end = MIN(grid.i_res_-1, i+1);
		j_end = MIN(grid.j_res_-1, j+1);
		k_end = MIN(grid.k_res_-1, k+1);

		for(int i=i_start; i<=i_end; i++) for(j=j_start; j<=j_end; j++) for(k=k_start; k<=k_end; k++)
		{
			int ix = grid.Index3Dto1D(i,j,k);

			TV3 cell_center = grid.CellCenterPosition(i,j,k);

			TS deviation = cell_center-pos;

			TS w = QuadBSplineKernel(deviation.x*grid.one_over_dx_)*
				    QuadBSplineKernel(deviation.y*grid.one_over_dx_)*
					QuadBSplineKernel(deviation.z*grid.one_over_dx_);

			weight_field.Set(ix, weight_field.Get(ix)+w);

			val_field.Set(ix, val_field.Get(ix)+val*w);
		}
	}	

	for(int p=0; p<grid.ijk_res_; p++)
	{
		TS weight = weight_field.Get(p);
		
		if(weight >= (TS)1e-06)
			val_field.Set(val_field.Get(p)/weight);
		else
			val_field.Set(TT());
	}
}