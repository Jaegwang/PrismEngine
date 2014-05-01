#pragma once

#include <atomic>
#include <iostream>
#include "GRID_UNIFORM_3D.h"

class DOMAIN_BLOCK
{ // it has no data
public:

	GRID_UNIFORM_3D grid_;

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

	GRID_UNIFORM_3D grid_;

	int i_res_, j_res_, k_res_;
	int ij_res_, jk_res_;
	int ijk_res_;

	TT* arr_;
	TT* arr_temp_;

	std::atomic<int> arr_ptr_;
	std::atomic<int> arr_temp_ptr_;
	
	DOMAIN_BLOCK** jk_blocks_;

	TT default_data_;

public:

	DOMAIN_DYNAMIC() : arr_ptr_(-1), arr_temp_ptr_(-1)
	{}
	~DOMAIN_DYNAMIC()
	{}

	void Initialize(const GRID_UNIFORM_3D& grid_input, const TT& default_data_input)
	{
		grid_ = grid_input;

		i_res_ = grid_.i_res_;
		j_res_ = grid_.j_res_;
		k_res_ = grid_.k_res_;

		ij_res_ = grid_.i_res_ * grid_.j_res_;
		jk_res_ = grid_.j_res_ * grid_.k_res_;

		ijk_res_ = grid_.ijk_res_;

		arr_ = new TT[1000000];
		arr_temp_ = new TT[1000000];

		arr_ptr_ = -1;
		arr_temp_ptr_ = -1;

		jk_blocks_ = new DOMAIN_BLOCK*[jk_res_];

		for(int i=0; i<jk_res_; i++) jk_blocks_[i] = 0;

		default_data_ = default_data_input;
	}

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
		DOMAIN_BLOCK* tail_block = 0;

		for(DOMAIN_BLOCK* block=jk_blocks_[jk_ix]; block != 0; block = block->next_block_)
		{
			if(block) tail_block = block;

			if(block->start_idx_ <= i && block->end_idx_ >= i)
			{
				arr_[block->start_ptr_+(i-block->start_idx_)] = data;
				return;			
			}
		
			if(block->next_block_ && block->end_idx_ < i && block->next_block_->start_idx_ > i)
			{
				DOMAIN_BLOCK* new_block = new DOMAIN_BLOCK; // TODO : create

				int ptr = atomic_fetch_add(&arr_ptr_, 1)+1;

				new_block->start_idx_ = i;
				new_block->end_idx_ = i;

				new_block->start_ptr_ = ptr;
				new_block->end_ptr_ = ptr;

				arr_[ptr] = data;

				new_block->prev_block_ = block;
				new_block->next_block_ = block->next_block_;

				block->next_block_->prev_block_ = new_block;
				block->next_block_ = new_block;				

				return;
			}
		}


		{
			DOMAIN_BLOCK* new_block = new DOMAIN_BLOCK;
			
			int ptr = atomic_fetch_add(&arr_ptr_, 1)+1;

			new_block->start_idx_ = i;
			new_block->end_idx_ = i;

			new_block->start_ptr_ = ptr;
			new_block->end_ptr_ = ptr;

			arr_[ptr] = data;

			if(tail_block)
			{
				new_block->prev_block_ = tail_block;
				tail_block->next_block_ = new_block;		
			}
			else
			{
				jk_blocks_[jk_ix] = new_block;
			}
		}
	}

	void RebuildDomain()
	{
		arr_ptr_ = -1;
		arr_temp_ptr_ = -1;		

		for(int p=0; p<jk_res_; p++)
		{
			jk_blocks_[p] = Rearrange(jk_blocks_[p], default_data_);		
		}

		TT* temp = arr_;
		arr_ = arr_temp_;
		arr_temp_ = temp;

		arr_ptr_ = (int)arr_temp_ptr_;
	}

	DOMAIN_BLOCK* Rearrange(DOMAIN_BLOCK* block, const TT& default_data)
	{
		if(!block) return 0;

		DOMAIN_BLOCK* block_head = 0;
		DOMAIN_BLOCK* block_tail = 0;

		int total_count = 0;

		for(DOMAIN_BLOCK* b=block; b != 0; b = b->next_block_)
		{
			for(int p=b->start_ptr_; p<=b->end_ptr_; p++)
			{
				if(arr_[p] != default_data)
				{
					total_count++;
				}
			}
		}

		int start_ptr = atomic_fetch_add(&arr_ptr_, total_count)+1;
		int count_pts = start_ptr;

		int start_idx = 0;

		int count = 0;
		int acc_count = 0;
		
		for(DOMAIN_BLOCK* b=block; b != 0; b = b->next_block_)
		{
			for(int p=b->start_ptr_; p<=b->end_ptr_; p++)
			{
				if(arr_[p] != default_data)
				{
					count++; acc_count++;
					arr_temp_[count_pts++] = arr_[p];
				}

				if(count == 1)
				{
					start_idx = b->start_idx_ + (p - b->start_ptr_);				
				}

				if(count > 0 && arr_[p] == default_data ||
				   count > 0 && p == b->end_ptr_ && b->next_block_ && b->end_idx_+1 < b->next_block_->start_idx_ ||
				   count > 0 && p == b->end_ptr_ && b->next_block_ == 0)
				{
					DOMAIN_BLOCK* new_block = new DOMAIN_BLOCK;

					new_block->start_idx_ = start_idx;
					new_block->end_idx_ = start_idx+count-1;

					new_block->start_ptr_ = start_ptr;
					new_block->end_ptr_ = start_ptr+count-1;

					start_ptr = new_block->end_ptr_+1;

					if(block_head) 
					{
						new_block->prev_block_ = block_tail;
						block_tail->next_block_ = new_block;
						block_tail = new_block;
					}
					else
					{
						block_head = new_block;
						block_tail = new_block;
					}

					count = 0;
				}
			}
		}

		if(total_count != acc_count) { std::cout<<"error : block countin critical !!!"<<std::endl; }

		return block_head;
	}

	void Render()
	{
		glDisable(GL_LIGHTING);
		glPointSize(3.0f);
		glBegin(GL_POINTS);
		glColor3f(1.0, 0.0, 0.0);
		for(int k=0; k<k_res_; k++) for(int j=0; j<j_res_; j++)
		{
			int ix = k*j_res_+j;
			for(DOMAIN_BLOCK* b=jk_blocks_[ix]; b != 0; b = b->next_block_)
			{
				int start_i = b->start_idx_;
				int end_i = b->end_idx_;

				Vec3 v0 = grid_.CellCenterPosition(start_i, j, k);		
				Vec3 v1 = grid_.CellCenterPosition(end_i, j, k);

				glVertex3fv(&v0.x);
				glVertex3fv(&v1.x);
			}
		}	
		glEnd();

		glBegin(GL_LINES);
		for(int k=0; k<k_res_; k++) for(int j=0; j<j_res_; j++)
		{
			int ix = k*j_res_+j;
			for(DOMAIN_BLOCK* b=jk_blocks_[ix]; b != 0; b = b->next_block_)
			{
				int start_i = b->start_idx_;
				int end_i = b->end_idx_;

				Vec3 v0 = grid_.CellCenterPosition(start_i, j, k);		
				Vec3 v1 = grid_.CellCenterPosition(end_i, j, k);

				glVertex3fv(&v0.x);
				glVertex3fv(&v1.x);
			}
		}	
		glEnd();

		glEnable(GL_LIGHTING);
	}
};
