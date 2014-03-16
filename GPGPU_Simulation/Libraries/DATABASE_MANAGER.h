
#pragma once
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

class DATABASE_MANAGER
{
private:

	std::map<std::string, int> idx_map_;

	std::vector<void*> ptr_arr_;
	std::vector<int>   num_arr_;
	std::vector<int>   size_arr_;

public:

	DATABASE_MANAGER();
	~DATABASE_MANAGER();

public:

	void* Insert(const std::string key_str, const int size, const int num)
	{
		void* ptr = malloc(size*num);

		num_arr_.push_back(num);
		size_arr_.push_back(size);
		ptr_arr_.push_back(ptr);

		idx_map_.insert(std::pair<std::string, int>(key_str, (int)(ptr_arr_.size()-1)));	

		return ptr;
	}

	void Erase(const std::string key_str)
	{
		auto it = idx_map_.find(key_str);
		if(it != idx_map_.end())
		{
			int ix = it->second;

			free(ptr_arr_[ix]);

			num_arr_[ix] = 0;
			size_arr_[ix] = 0;
			ptr_arr_[ix] = 0;		
		}	
	}

	void* Find(const std::string key_str)
	{
		auto it = idx_map_.find(key_str);
		if(it != idx_map_.end())
		{
			return ptr_arr_[it->second];
		}
		return 0;
	}
};