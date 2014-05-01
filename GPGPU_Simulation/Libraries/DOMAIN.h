#pragma once


class DOMAIN_BLOCK
{ // it has no data
public:

	int start_idx_;
	int end_idx_;

	int start_ptr_;
	int end_ptr_;

	DOMAIN_BLOCK* prev_block_;
	DOMAIN_BLOCK* next_block_;

	DOMAIN_BLOCK()
		: start_idx_(-1), end_idx_(-1), start_ptr_(-1), end_ptr_(-1), prev_block_(0), next_block_(0)
	{}
	~DOMAIN_BLOCK()
	{}
};

template<class TT>
class DOMAIN_DYNAMIC
{
public:

	int i_res_, j_res_, k_res_;
	int ij_res_, jk_res_;
	int ijk_res_;

	TT* arr_;

	DOMAIN_BLOCK** jk_blocks_;

public:

	DOMAIN_BLOCK* SearchBlock(const int i, const int j, const int k)
	{
		int jk_ix = k*j_res_ + j;

		for(DOMAIN_BLOCK* b = jk_blocks_[jk_ix]; b != 0; b = b->next_block_)
		{
			if(block->start_idx_ <= i && block->end_idx_ >= i) return b;
		}

		return 0;
	}

	TT* SearchData(const int i, const int j, const int k)
	{
		DOMAIN_BLOCK* block = SearchBlock(i,j,k);
		if(block == 0) return 0;

		return arr_ + block->start_ptr_ + (i-block->start_idx_);
	}

	void Insert(const int i, const int j, const int k, const TT& data)
	{
		int jk_ix = k*j_res_ + j;

		DOMAIN_BLOCK* block = 0;
		for(block = jk_blocks_[jk_ix]; block != 0; block = block->next_block_)
		{
//			if(block->start_idx_ <= i && block->end_idx_ >= i) break;


		}

	}
};
