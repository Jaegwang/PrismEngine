
#pragma once

class FIELD_BLOCK
{ 
public:

	int start_idx_;
	int end_idx_;

	int start_ptr_;
	int end_ptr_;

	FIELD_BLOCK* prev_block_;
	FIELD_BLOCK* next_block_;

public:

	FIELD_BLOCK()
		: start_idx_(-1), end_idx_(-1), start_ptr_(-1), end_ptr_(-1), prev_block_(0), next_block_(0)
	{}
	~FIELD_BLOCK()
	{}
};