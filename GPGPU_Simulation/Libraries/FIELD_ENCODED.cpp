
#include "FIELD_ENCODED.h"
#include "GENERAL_MACRO.h"

template<class TT>
void FIELD_ENCODED<TT>::Initialize(const GRID& grid_input, const TT& default_data_input)
{
	Finalize();
	table_k_res_ = 16;

	int km = grid_input.k_res_/table_k_res_;
	km += (grid_input.k_res_%table_k_res_ > 0) ? 1 : 0;

	k_res_ = table_k_res_*km;
	block_k_res_ = k_res_/table_k_res_;
	
	grid_.Initialize(grid_input.min_, grid_input.max_, k_res_-grid_input.ghost_width_*2, grid_input.ghost_width_);

	assert(grid_.k_res_ == k_res_);

	i_res_ = grid_.i_res_;
	j_res_ = grid_.j_res_;

	ij_res_ = i_res_*j_res_;
	ijk_res_ = ij_res_*k_res_;

	default_data_ = default_data_input;

	int table_size = table_k_res_*j_res_*i_res_;
	table_ = new std::atomic<CACHE_BLOCK<TT>*>[table_size];
	block_array_ = new CACHE_BLOCK<TT>*[table_size];

	for(int t=0; t<table_size; t++)
	{
		table_[t] = 0;
		block_array_[t] = 0;
	}
}

template<class TT>
void FIELD_ENCODED<TT>::Finalize()
{
	FOR_EACH_PARALLEL(k, 0, table_k_res_-1)
	for(int j=0; j<j_res_; j++) for(int i=0; i<i_res_; i++)
	{
		int t_ix = k*ij_res_ + j*i_res_ + i;

		CACHE_BLOCK<TT>* tb = table_[t_ix];
		if(tb) delete tb;
	}

	if(table_) { delete table_; table_ = 0; };
	if(block_array_) { delete block_array_; block_array_ = 0; };

	num_blocks_ = 0;
	count_blocks_ = 0;
}

template<class TT>
void FIELD_ENCODED<TT>::Set(const int idx, const TT data)
{
	int i,j,k;
	grid_.Index1Dto3D(idx,i,j,k);

	Set(i,j,k,data);
}

template<class TT>
void FIELD_ENCODED<TT>::Set(const int i, const int j, const int k, const TT data)
{
	const int t_k = k / block_k_res_;
	const int b_k = k % block_k_res_;

	const int t_idx = t_k*ij_res_ + j*i_res_ + i;

	CACHE_BLOCK<TT>* block = table_[t_idx];
	CACHE_BLOCK<TT>* null_block = 0;

	if(!block && data != default_data_)
	{
		CACHE_BLOCK<TT>* new_block = new CACHE_BLOCK<TT>(block_k_res_, default_data_);

		new_block->i_ = i;
		new_block->j_ = j;
		new_block->k_ = t_k;

		if(table_[t_idx].compare_exchange_weak(null_block, new_block, std::memory_order_release, std::memory_order_relaxed) == true)
		{
			std::atomic_fetch_add(&num_blocks_, 1);		
		}
		else
		{
			delete new_block;
		}

		block = table_[t_idx];
	}

	if(block) block->array_[b_k] = data;
}

template<class TT>
TT FIELD_ENCODED<TT>::Get(const int idx) const
{
	int i,j,k;
	grid_.Index1Dto3D(idx,i,j,k);

	return Get(i,j,k);
}

template<class TT>
TT FIELD_ENCODED<TT>::Get(const int i, const int j, const int k) const
{
	const int t_k = k / block_k_res_;
	const int b_k = k % block_k_res_;

	const int t_idx = t_k*ij_res_ + j*i_res_ + i;
	CACHE_BLOCK<TT>* block = table_[t_idx];

	if(!block) return default_data_;

	return block->array_[b_k];
}

template<class TT>
TT FIELD_ENCODED<TT>::Get(const TV3& p) const
{
	return TriLinearInterpolate(p);
}

template<class TT>
void FIELD_ENCODED<TT>::Rebuild()
{
	count_blocks_ = 0;

	FOR_EACH_PARALLEL(k, 0, table_k_res_-1)
	for(int j=0; j<j_res_; j++) for(int i=0; i<i_res_; i++)
	{
		int t_ix = k*ij_res_ + j*i_res_ + i;
		bool is_release = true;

		CACHE_BLOCK<TT>* tb = table_[t_ix];
		if(!tb) continue;

		for(int t=0; t<tb->size_; t++)
		{
			if(tb->array_[t] != default_data_)
			{
				is_release = false;
				break;
			}							
		}

		if(is_release == true)
		{
			delete tb;
			table_[t_ix] = 0;

			std::atomic_fetch_sub(&num_blocks_,1);
		}
		else
		{
			int c_ix = std::atomic_fetch_add(&count_blocks_, 1);
			block_array_[c_ix] = tb;
		}
	}	

	assert(num_blocks_ == count_blocks_);

	if(num_blocks_ != count_blocks_)
	{	
		std::cout<<"Error count"<<std::endl;
	}
}

template<class TT>
TT FIELD_ENCODED<TT>::TriLinearInterpolate(const TV3& p) const
{ //http://en.wikipedia.org/wiki/Trilinear_interpolation
	TV3 cp = grid_.Clamp(p);

	int i,j,k;
	grid_.LeftBottomIndex(cp,i,j,k);

	TV3 p_0 = grid_.CellCenterPosition(i,j,k);

	TS x_d = (cp.x-p_0.x)/grid_.dx_;
	TS y_d = (cp.y-p_0.y)/grid_.dx_;
	TS z_d = (cp.z-p_0.z)/grid_.dx_;

	TT c_00 = Get(i,j,k)*((TS)1-x_d)     + Get(i+1,j,k)*(x_d);
	TT c_10 = Get(i,j+1,k)*((TS)1-x_d)   + Get(i+1,j+1,k)*(x_d);
	TT c_01 = Get(i,j,k+1)*((TS)1-x_d)   + Get(i+1,j,k+1)*(x_d);
	TT c_11 = Get(i,j+1,k+1)*((TS)1-x_d) + Get(i+1,j+1,k+1)*(x_d);

	TT c_0 = c_00*((TS)1-y_d) + c_10*y_d;
	TT c_1 = c_01*((TS)1-y_d) + c_11*y_d;

	return c_0*((TS)1-z_d) + c_1*z_d;
}

template<class TT>
void FIELD_ENCODED<TT>::Render()
{
	glDisable(GL_LIGHTING);

	glBegin(GL_LINES);
	glColor3f(0.0, 1.0, 0.0);
	for(int k=0; k<table_k_res_; k++) 
	for(int j=0; j<j_res_; j++)  
	for(int i=0; i<i_res_; i++)
	{
		int ix = k*ij_res_ + j*i_res_ + i;

		CACHE_BLOCK<TT>* block = table_[ix];

		if(block == 0) continue;

		TV3 v0 = grid_.CellCenterPosition(block->i_, block->j_, block->k_*block_k_res_);
		TV3 v1 = grid_.CellCenterPosition(block->i_, block->j_, (block->k_+1)*block_k_res_);

		glVertex3fv(&v0.x);
		glVertex3fv(&v1.x);
		
	}	
	glEnd();

	glEnable(GL_LIGHTING);
}

template class FIELD_ENCODED<float>;
template class FIELD_ENCODED<double>;
template class FIELD_ENCODED<int>;
template class FIELD_ENCODED<TV3>;
template class FIELD_ENCODED<bool>;