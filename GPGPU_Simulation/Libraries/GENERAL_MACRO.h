
#pragma once

#define SAFE_DELETE(_arr) {if(_arr){delete _arr; _arr=0;}}
#define SAFE_DELETE_ARRAY(_arr) {if(_arr){delete[] _arr; _arr=0;}}