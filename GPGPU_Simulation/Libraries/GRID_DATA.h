
#include <vector>
#include <stack>
#include "MATH_DEFINITION.h"

template<class TT>
class GRID_DYNAMIC
{
private:

	TT** compressed_arr_;

	std::stack<TT*> free_blocks_;
	std::stack<TT*> used_blocks_;

	TT default_data_;

	int num_blocks_;

public:

	Vec3 min_, max_;

	FLT dx_, dy_, dz_;

	FLT gx_, gy_, gz_;

	FLT one_over_dx_, one_over_dy_, one_over_dz_;

	int i_res_, j_res_, k_res_;
	int ij_res_, ijk_res_;

	int ghost_width_;

	int block_i_res_;
	int block_size_;

public:

	GRID_DYNAMIC()
	{};
	~GRID_DYNAMIC()
	{};

public:

	void Initialize(const Vec3& min, const Vec3& max, int i, int j, int k, int g, int block_i=16);
	void Finalize();

	void Reset()
	{
		while(used_blocks_.empty() == false)
		{
			free_blocks_.push(used_blocks_.top());
			used_blocks_.pop();
		}	
	}

	Vec3 CellCenterPosition(const int i, const int j, const int k)
	{
		return min_ + Vec3(((FLT)i+(FLT)0.5)*dx_, ((FLT)j+(FLT)0.5)*dy_, ((FLT)k+(FLT)0.5)*dz_);
	}

	TT GetData(const int i, const int j, const int k)
	{
		int b = i/block_size_;
		int b_ix = k*(block_i_res_*j_res_) + j*block_i_res_ + b;

		if(compressed_arr_[b_ix] == 0) return default_data_;

		int p = i - b*block_size_;
		return (compressed_arr_[b_ix])[p];
	}

	void SetData(const int i, const int j, const int k, const TT& data)
	{
		int b = i/block_size_;
		int b_ix = k*(block_i_res_*j_res_) + j*block_i_res_ + b;

		if(compressed_arr_[b_ix] == 0)
		{
			compressed_arr_[b_ix] = free_blocks_.top();
			free_blocks_.pop();
			used_blocks_.push(compressed_arr_[b_ix]);
		}

		int p = i - b*block_size_;
		(compressed_arr_[b_ix])[p] = data;
	}

	void Render();
};