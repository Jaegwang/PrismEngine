
#pragma once

#include <condition_variable>
#include <assert.h>
#include <amp.h>
#include <vector>
#include <thread>

using namespace concurrency;

static int _CPU_MAX_THREADS = std::thread::hardware_concurrency();
static std::condition_variable _thread_cv;
static std::mutex _mutex_cv;

#define BEGIN_CPU_THREADS_1D(_i_res) { vector<std::thread> threads; int num_elements = _i_res; int _num_threads = MIN(_CPU_MAX_THREADS, num_elements); \
												const int quotient = num_elements / _num_threads, remainder = num_elements%_num_threads; atomic<int> _sync_threads = 0; \
												for (int ix_base = 0, id = 0; id < _num_threads; ++id) {\
												threads.push_back(std::thread([&](){ \
												const int ix_offset = id < remainder ? (quotient + 1) : quotient; \
												const int ix_begin = ix_base; ix_base += ix_offset; \
												const int ix_end = ix_begin + ix_offset - 1; \
												const int thread_id = id;
#define END_CPU_THREADS_1D		     })); } for(auto& thread : threads) { thread.join(); }}


#define THREAD_LOOPS_1D(_i) for(int _i = ix_begin; _i <= ix_end; _i++)


#define SYNC_CPU_THREADS {std::unique_lock<std::mutex> _lk(_mutex_cv); _sync_threads++;\
						  if (_sync_threads == _num_threads) { _sync_threads = 0; _thread_cv.notify_all(); }\
						  else { _thread_cv.wait(_lk); }}


#define BEGIN_PARALLEL_FOR_EACH_1D(_i_res, _i) { concurrency::extent<1> view_ext(_i_res);  parallel_for_each(view_ext, [=](index<1> idx) restrict(amp) { int _i = idx[0];
#define END_PARALLEL_FOR_EACH_1D               });}

