#pragma once

#include <atomic>
#include <iostream>
#include "GRID.h"
#include "CACHE_BLOCK.h"
#include "FIELD.h"

template<class TT>
class FIELD_ENCODED : public FIELD<TT>
{
private:

	int i_res_, j_res_, k_res_;
	int ij_res_, ijk_res_;

	int table_k_res_;
	int block_k_res_;

	std::atomic<CACHE_BLOCK<TT>*>* table_;
	
	TT default_data_;

	std::atomic<int> num_blocks_;

public:

	FIELD_ENCODED()	: table_(0), num_blocks_(0)
	{}

	~FIELD_ENCODED()
	{
		Finalize();	
	}

private:

	TT TriLinearInterpolate(const TV3& p) const;

public:

	void Initialize(const GRID& grid_input, const TT& default_data_input);
	void Finalize();
	
	void Set(const int idx, const TT data);
	void Set(const int i, const int j, const int k, const TT data);

	TT   Get(const int idx) const;
	TT   Get(const int i, const int j, const int k) const;
	TT   Get(const TV3& p) const;

	void Rebuild(){};

	// for test
	void Render();
};