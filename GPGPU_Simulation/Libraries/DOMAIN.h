#pragma once

#include <atomic>
#include <iostream>
#include "GRID_UNIFORM_3D.h"
#include "ARRAY_ATOMIC.h"

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

	GRID_UNIFORM_3D grid_;

	int i_res_, j_res_, k_res_;
	int ij_res_, jk_res_;
	int ijk_res_;

	ARRAY_ATOMIC<TT>* arr_;
	ARRAY_ATOMIC<TT>* arr_temp_;

	ARRAY_ATOMIC<DOMAIN_BLOCK*> stack_block_;
	
	DOMAIN_BLOCK** jk_blocks_;
	DOMAIN_BLOCK** jk_blocks_temp_;

	TT default_data_;

public:

	DOMAIN_DYNAMIC()
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

		int num = 8;

		arr_ = new ARRAY_ATOMIC<TT>;
		arr_temp_ = new ARRAY_ATOMIC<TT>;

		arr_->Initialize(num);
		arr_temp_->Initialize(num);

		stack_block_.Initialize(num);
		int ix = stack_block_.Push(num);

		for(int p=0; p<num; p++)
		{
			stack_block_(ix++) = new DOMAIN_BLOCK();
		}

		jk_blocks_ = new DOMAIN_BLOCK*[jk_res_];
		jk_blocks_temp_ = new DOMAIN_BLOCK*[jk_res_];

		for(int i=0; i<jk_res_; i++) jk_blocks_[i] = 0;

		default_data_ = default_data_input;
	}

	TT Find(const int i, const int j, const int k)
	{
		int jk_ix = k*j_res_ + j;
		for(DOMAIN_BLOCK* block=jk_blocks_[jk_ix]; block != 0; block = block->next_block_)
		{
			if(block->start_idx_ <= i && block->end_idx_ >= i)
			{
				return (*arr_)(block->start_ptr_+(i-block->start_idx_));
			}
		}	

		return default_data_;
	}

	void Insert(const int i, const int j, const int k, const TT& data)
	{
		ReloadBlocks();

		int jk_ix = k*j_res_ + j;
		DOMAIN_BLOCK* tail_block = 0;

		for(DOMAIN_BLOCK* block=jk_blocks_[jk_ix]; block != 0; block = block->next_block_)
		{
			if(block) tail_block = block;

			if(block->start_idx_ <= i && block->end_idx_ >= i)
			{
				(*arr_)(block->start_ptr_+(i-block->start_idx_)) = data;
				return;			
			}
		
			if(block->next_block_ && block->end_idx_ < i && block->next_block_->start_idx_ > i)
			{
				DOMAIN_BLOCK* new_block = stack_block_(stack_block_.Pop());
				(*new_block) = DOMAIN_BLOCK();

				int ptr = arr_->Push(1);
				(*arr_)(ptr) = data;

				new_block->start_idx_ = i;
				new_block->end_idx_ = i;

				new_block->start_ptr_ = ptr;
				new_block->end_ptr_ = ptr;

				new_block->prev_block_ = block;
				new_block->next_block_ = block->next_block_;

				block->next_block_->prev_block_ = new_block;
				block->next_block_ = new_block;				

				return;
			}
		}

		// for tail block
		{
			DOMAIN_BLOCK* new_block = stack_block_(stack_block_.Pop());
			(*new_block) = DOMAIN_BLOCK();			

			int ptr = arr_->Push(1);
			(*arr_)(ptr) = data;

			new_block->start_idx_ = i;
			new_block->end_idx_ = i;

			new_block->start_ptr_ = ptr;
			new_block->end_ptr_ = ptr;

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

	void ReloadBlocks()
	{
		if(stack_block_.Size() > 0) return;

		int num = stack_block_.Total();
		int ix = stack_block_.Push(num);

		#pragma omp parallel for
		for(int i=ix; i<=ix+num-1; i++) stack_block_(i) = new DOMAIN_BLOCK();
	}

	void RebuildDomain()
	{
		arr_temp_->Reset();

		DOMAIN_BLOCK** jk_temp;
		SWAP(jk_blocks_, jk_blocks_temp_, jk_temp);

		#pragma omp parallel for
		for(int p=0; p<jk_res_; p++)
		{
			jk_blocks_[p] = Rearrange(jk_blocks_temp_[p], default_data_);		
		}

		#pragma omp parallel for
		for(int p=0; p<jk_res_; p++)
		{
			for(DOMAIN_BLOCK* block = jk_blocks_temp_[p]; block != 0; block = block->next_block_)
			{
				stack_block_(stack_block_.Push()) = block;
			}
		}

		ARRAY_ATOMIC<TT>* tt_temp;
		SWAP(arr_, arr_temp_, tt_temp);
	}

	DOMAIN_BLOCK* Rearrange(DOMAIN_BLOCK* block, const TT& default_data)
	{
		if(!block) return 0;

		ReloadBlocks();

		DOMAIN_BLOCK* block_head = 0;
		DOMAIN_BLOCK* block_tail = 0;

		int total_count = 0;

		for(DOMAIN_BLOCK* b = block; b != 0; b = b->next_block_)
		{
			for(int p=b->start_ptr_; p<=b->end_ptr_; p++)
			{
				if((*arr_)(p) != default_data) total_count++;
			}
		}

		int start_ptr = arr_temp_->Push(total_count);
		int count_ptr = start_ptr;

		int start_idx = 0;

		int count = 0;
		int acc_count = 0;
		
		for(DOMAIN_BLOCK* b = block; b != 0; b = b->next_block_)
		{
			for(int p=b->start_ptr_; p<=b->end_ptr_; p++)
			{
				if((*arr_)(p) != default_data)
				{
					count++; acc_count++;
					(*arr_temp_)(count_ptr++) = (*arr_)(p);
				}

				if(count == 1)
				{
					start_idx = b->start_idx_ + (p - b->start_ptr_);				
				}

				if(count > 0 && (*arr_)(p) == default_data ||
				   count > 0 && p == b->end_ptr_ && b->next_block_ && b->end_idx_+1 < b->next_block_->start_idx_ ||
				   count > 0 && p == b->end_ptr_ && b->next_block_ == 0)
				{
					DOMAIN_BLOCK* new_block = stack_block_(stack_block_.Pop());
					(*new_block) = DOMAIN_BLOCK();	

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
