
#pragma once

#define SAFE_DELETE(_arr) {if(_arr){delete _arr; _arr=0;}}
#define SAFE_DELETE_ARRAY(_arr) {if(_arr){delete[] _arr; _arr=0;}}

#define FOR_EACH_PARALLEL(_it, _i0, _i1) __pragma(omp parallel for) for(int _it=_i0; _it<=_i1; _it++)