#pragma once

#include <atomic>
#include <iostream>
#include "GRID_UNIFORM_3D.h"
#include "FIELD_BLOCK.h"
#include "ARRAY_ATOMIC.h"

template<class TT>
class FIELD_DYNAMIC
{
public:

	GRID_UNIFORM_3D grid_;

	int i_res_, j_res_, k_res_;
	int ij_res_, jk_res_;
	int ijk_res_;

	ARRAY_ATOMIC<TT>* arr_;
	ARRAY_ATOMIC<TT>* arr_temp_;

	ARRAY_ATOMIC<FIELD_BLOCK*> stack_block_;
	
	FIELD_BLOCK** jk_blocks_;
	FIELD_BLOCK** jk_blocks_temp_;

	TT default_data_;

public:

	FIELD_DYNAMIC() : arr_(0), arr_temp_(0), jk_blocks_(0), jk_blocks_temp_(0)
	{}
	~FIELD_DYNAMIC()
	{
		Finalize();	
	}

private:

	void Insert (const int i, const int j, const int k, const TT& data);
	TT   Find   (const int i, const int j, const int k);

	FIELD_BLOCK* Rearrange(FIELD_BLOCK* block, const TT& default_data);

public:

	void Initialize(const GRID_UNIFORM_3D& grid_input, const TT& default_data_input);
	void Finalize();
	
	void Set(const int i, const int j, const int k, const TT& data);
	TT   Get(const int i, const int j, const int k);
	TT   Get(const Vec3& p);

	void RebuildField();
	void ReloadBlocks();

	TT TriLinearInterpolate(const Vec3& p);

	// for test
	void Render();
};