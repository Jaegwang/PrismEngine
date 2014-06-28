
#include <vector>
#include <stack>
#include "MATH_DEFINITION.h"

class GRID_DYNAMIC
{
private:

	int* compressed_arr_;

public:

	TV3 min_, max_;

	TS dx_, dy_, dz_, dw_;

	TS gx_, gy_, gz_, gw_;

	TS one_over_dx_, one_over_dy_, one_over_dz_;

	int i_res_, j_res_, k_res_;
	int ij_res_, ijk_res_;
	int b_res_;

	int ghost_width_;

	int block_size_;

	std::atomic<int> mem_ptr_;

public:

	GRID_DYNAMIC() : compressed_arr_(0), mem_ptr_(0)
	{};
	~GRID_DYNAMIC()
	{};

public:

	void Initialize(const TV3& min, const TV3& max, int i, int j, int k, int g, int block_i=16);
	void Initialize(const TV3& min, const TV3& max, const TS dw, const int g, const int block_size=16);
	void Finalize();

	TV3 CellCenterPosition(const int i, const int j, const int k)
	{
		return min_ + TV3(((TS)i+(TS)0.5)*dx_, ((TS)j+(TS)0.5)*dy_, ((TS)k+(TS)0.5)*dz_);
	}

	template<class TT>
	TT GetData(const int i, const int j, const int k, TT* arr, const TT& data)
	{
		int b = i/block_size_;
		int b_ix = k*(b_res_*j_res_) + j*b_res_ + b;

		if(compressed_arr_[b_ix] < 0) return data;

		int p = i - b*block_size_;
		return arr[compressed_arr_[b_ix]+p];
	}

	template<class TT>
	void SetData(const int i, const int j, const int k, TT* arr, const TT& data)
	{
		int b = i/block_size_;
		int b_ix = k*(b_res_*j_res_) + j*b_res_ + b;

		if(compressed_arr_[b_ix] < 0)
		{
			compressed_arr_[b_ix] = mem_ptr_;
			mem_ptr_ += block_size_;
		}

		int p = i - b*block_size_;
		arr[compressed_arr_[b_ix]+p] = data;
	}

	void Render();
};