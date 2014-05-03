
#include "FIELD_DYNAMIC.h"

template<class TT>
void FIELD_DYNAMIC<TT>::Initialize(const GRID& grid_input, const TT& default_data_input)
{
	Finalize();

	grid_ = grid_input;

	i_res_ = grid_.i_res_; j_res_ = grid_.j_res_; k_res_ = grid_.k_res_;

	jk_res_ = grid_.j_res_ * grid_.k_res_;

	ijk_res_ = grid_.ijk_res_;

	arr_ = new ARRAY_DYNAMIC<TT>;
	arr_temp_ = new ARRAY_DYNAMIC<TT>;

	arr_->Initialize(jk_res_);
	arr_temp_->Initialize(jk_res_);

	stack_block_.Initialize(jk_res_);
	int ix = stack_block_.Push(jk_res_);

	for(int p=0; p<jk_res_; p++) stack_block_(ix++) = new FIELD_BLOCK();

	jk_blocks_ = new FIELD_BLOCK*[jk_res_];
	jk_blocks_temp_ = new FIELD_BLOCK*[jk_res_];

	for(int i=0; i<jk_res_; i++) jk_blocks_[i] = jk_blocks_temp_[i] = 0;

	default_data_ = default_data_input;
}

template<class TT>
void FIELD_DYNAMIC<TT>::Finalize()
{
	int tot = stack_block_.Size();
	for(int i=0; i<tot; i++)
	{
		FIELD_BLOCK* b = stack_block_(i);
		delete b;
	}

	if(arr_)      delete arr_;
	if(arr_temp_) delete arr_temp_;

	if(jk_blocks_)      delete[] jk_blocks_;
	if(jk_blocks_temp_) delete[] jk_blocks_temp_;
}

template<class TT>
void FIELD_DYNAMIC<TT>::Insert(const int i, const int j, const int k, const TT& data)
{
	ReloadBlocks();

	int jk_ix = k*j_res_ + j;
	FIELD_BLOCK* tail_block = 0;

	for(FIELD_BLOCK* block=jk_blocks_[jk_ix]; block != 0; block = block->next_block_)
	{
		if(block) tail_block = block;

		if(block->start_idx_ <= i && block->end_idx_ >= i)
		{
			(*arr_)(block->start_ptr_+(i-block->start_idx_)) = data;
			return;			
		}

		if(data == default_data_) continue;		

		if(block->next_block_ && block->end_idx_ < i && block->next_block_->start_idx_ > i)
		{
			FIELD_BLOCK* new_block = stack_block_(stack_block_.Pop());
			(*new_block) = FIELD_BLOCK();

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
		if(data == default_data_) return;

		FIELD_BLOCK* new_block = stack_block_(stack_block_.Pop());
		(*new_block) = FIELD_BLOCK();			

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

template<class TT>
TT FIELD_DYNAMIC<TT>::Find(const int i, const int j, const int k)
{
	int jk_ix = k*j_res_ + j;
	for(FIELD_BLOCK* block=jk_blocks_[jk_ix]; block != 0; block = block->next_block_)
	{
		if(block->start_idx_ <= i && block->end_idx_ >= i)
		{
			return (*arr_)(block->start_ptr_+(i-block->start_idx_));
		}
	}
	return default_data_;
}

template<class TT>
void FIELD_DYNAMIC<TT>::ReloadBlocks()
{
	if(stack_block_.Size() > 0) return;

	int num = stack_block_.Total();
	int ix = stack_block_.Push(num);

	#pragma omp parallel for
	for(int i=ix; i<=ix+num-1; i++) stack_block_(i) = new FIELD_BLOCK();
}

template<class TT>
void FIELD_DYNAMIC<TT>::RebuildField()
{
	arr_temp_->Reset();

	FIELD_BLOCK** jk_temp;
	SWAP(jk_blocks_, jk_blocks_temp_, jk_temp);

	#pragma omp parallel for
	for(int p=0; p<jk_res_; p++)
	{
		jk_blocks_[p] = Rearrange(jk_blocks_temp_[p], default_data_);		
	}

	#pragma omp parallel for
	for(int p=0; p<jk_res_; p++)
	{
		for(FIELD_BLOCK* block = jk_blocks_temp_[p]; block != 0; block = block->next_block_)
		{
			stack_block_(stack_block_.Push()) = block;
		}
	}

	ARRAY_DYNAMIC<TT>* tt_temp;
	SWAP(arr_, arr_temp_, tt_temp);
}

template<class TT>
FIELD_BLOCK* FIELD_DYNAMIC<TT>::Rearrange(FIELD_BLOCK* block, const TT& default_data)
{
	if(!block) return 0;

	ReloadBlocks();

	FIELD_BLOCK* block_head = 0;
	FIELD_BLOCK* block_tail = 0;

	int total_count = 0;

	for(FIELD_BLOCK* b = block; b != 0; b = b->next_block_)
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
		
	for(FIELD_BLOCK* b = block; b != 0; b = b->next_block_)
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

			if( count > 0 && (*arr_)(p) == default_data ||
				count > 0 && p == b->end_ptr_ && b->next_block_ && b->end_idx_+1 < b->next_block_->start_idx_ ||
				count > 0 && p == b->end_ptr_ && b->next_block_ == 0 )
			{
				FIELD_BLOCK* new_block = stack_block_(stack_block_.Pop());
				(*new_block) = FIELD_BLOCK();	

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

template<class TT>
void FIELD_DYNAMIC<TT>::Set(const int idx, const TT& data)
{
	int i,j,k;
	grid_.Index1Dto3D(idx,i,j,k);
	Insert(i,j,k, data);
}

template<class TT>
void FIELD_DYNAMIC<TT>::Set(const int i, const int j, const int k, const TT& data)
{
	Insert(i,j,k,data);
}

template<class TT>
TT FIELD_DYNAMIC<TT>::Get(const int idx)
{
	int i,j,k;
	grid_.Index1Dto3D(idx,i,j,k);
	return Find(i,j,k);
}

template<class TT>
TT FIELD_DYNAMIC<TT>::Get(const int i, const int j, const int k)
{
	return Find(i,j,k);
}

template<class TT>
TT FIELD_DYNAMIC<TT>::Get(const Vec3& p)
{
	return TriLinearInterpolate(p);
}

template<class TT>
TT FIELD_DYNAMIC<TT>::TriLinearInterpolate(const Vec3& p)
{ //http://en.wikipedia.org/wiki/Trilinear_interpolation
	Vec3 cp = grid_.Clamp(p);

	int i,j,k;
	grid_.LeftBottomIndex(cp,i,j,k);

	Vec3 p_0 = grid_.CellCenterPosition(i,j,k);

	FLT x_d = (cp.x-p_0.x)/grid_.dx_;
	FLT y_d = (cp.y-p_0.y)/grid_.dx_;
	FLT z_d = (cp.z-p_0.z)/grid_.dx_;

	TT c_00 = Get(i,j,k)*((FLT)1-x_d)     + Get(i+1,j,k)*(x_d);
	TT c_10 = Get(i,j+1,k)*((FLT)1-x_d)   + Get(i+1,j+1,k)*(x_d);
	TT c_01 = Get(i,j,k+1)*((FLT)1-x_d)   + Get(i+1,j,k+1)*(x_d);
	TT c_11 = Get(i,j+1,k+1)*((FLT)1-x_d) + Get(i+1,j+1,k+1)*(x_d);

	TT c_0 = c_00*((FLT)1-y_d) + c_10*y_d;
	TT c_1 = c_01*((FLT)1-y_d) + c_11*y_d;

	return c_0*((FLT)1-z_d) + c_1*z_d;
}

template<class TT>
void FIELD_DYNAMIC<TT>::Render()
{
	glDisable(GL_LIGHTING);
	glPointSize(3.0f);
	glBegin(GL_POINTS);
	glColor3f(1.0, 0.0, 0.0);
	for(int k=0; k<k_res_; k++) for(int j=0; j<j_res_; j++)
	{
		int ix = k*j_res_+j;
		for(FIELD_BLOCK* b=jk_blocks_[ix]; b != 0; b = b->next_block_)
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
	glColor3f(0.0, 1.0, 0.0);
	for(int k=0; k<k_res_; k++) for(int j=0; j<j_res_; j++)
	{
		int ix = k*j_res_+j;
		for(FIELD_BLOCK* b=jk_blocks_[ix]; b != 0; b = b->next_block_)
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

template class FIELD_DYNAMIC<float>;
template class FIELD_DYNAMIC<double>;
template class FIELD_DYNAMIC<int>;
template class FIELD_DYNAMIC<Vec3>;