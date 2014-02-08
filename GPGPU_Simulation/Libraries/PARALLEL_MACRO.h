
#pragma once

#include <condition_variable>
#include <assert.h>
#include <amp.h>
#include <vector>
#include <thread>
#include <algorithm>

using namespace concurrency;
using namespace std;

static const int _AMP_TILE_SIZE = 1024;
static const int _AMP_MAX_THREADS = 30000000;

static int _CPU_MAX_THREADS = std::thread::hardware_concurrency();
static std::condition_variable _thread_cv;
static std::mutex _mutex_cv;

#define BEGIN_CPU_THREADS_1D(_i_res) { vector<std::thread> threads; int _num_threads = MIN(_CPU_MAX_THREADS, (int)_i_res); const int _quotient = (int)_i_res/_num_threads; \
												int _remainder=(int)_i_res%_num_threads, _ix_base=0, _thread_size=0; atomic<int> _sync_threads=0; \
												for(int id = 0; id < _num_threads; ++id) { _thread_size = id < _remainder ? (_quotient + 1) : _quotient; \
													auto _thread_func = [&](int tid, int ix_begin, int ix_end){
													// Thread Work.
#define END_CPU_THREADS_1D }; threads.push_back(std::thread(_thread_func, id, _ix_base, _ix_base+_thread_size-1)); _ix_base+=_thread_size; } for(auto& thread : threads) { thread.join(); }}


#define THREAD_LOOPS_1D(_i) for(int _i = ix_begin; _i <= ix_end; _i++)


#define SYNC_CPU_THREADS {std::unique_lock<std::mutex> _lk(_mutex_cv); _sync_threads++;\
						  if (_sync_threads == _num_threads) { _sync_threads = 0; _thread_cv.notify_all(); } else { _thread_cv.wait(_lk); }}


#define BEGIN_PARALLEL_FOR_EACH_1D(_i_res, _i) { int amp_res = (int)_i_res; int amp_num = (amp_res-1)/_AMP_MAX_THREADS+1; int amp_start_idx=0; \
													for (int n = 0; n < amp_num; n++) { \
														int size = MIN(_AMP_MAX_THREADS, amp_res-amp_start_idx); concurrency::extent<1> ext(size); \
														parallel_for_each(ext.tile<_AMP_TILE_SIZE>().pad(), [=](tiled_index<_AMP_TILE_SIZE> tidx) restrict(amp) { \
															int _i = amp_start_idx + tidx.global[0]; if (tidx.global[0] >= size) return;
#define END_PARALLEL_FOR_EACH_1D			   }); amp_start_idx += size; }};