#pragma once

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <io.h>
#include <string.h>
#include <fstream>

#include <windows.h>


class FILE_UTILITIES
{
public:
	FILE_UTILITIES()
	{}
	~FILE_UTILITIES()
	{}

public:

	int Execute_Process(const std::string &cmd)
	{
		return system(cmd.c_str());
	}

	bool Directory_Exists(const std::string& dirname)
	{DWORD attr=GetFileAttributes(dirname.c_str());return((attr!=-1)&&(attr&FILE_ATTRIBUTE_DIRECTORY));}

	bool Create_Directory(const std::string& dirname)
	{
		if(!Directory_Exists(dirname))
		{
			//std::cout<<"Creating directory using CreateDirectory...";
			CreateDirectory(dirname.c_str(),0);
			if(GetLastError()==ERROR_ALREADY_EXISTS) std::cout<<"ERROR_ALREADY_EXISTS  : "<<dirname.c_str()<<std::endl;//debug
			if(GetLastError()==ERROR_PATH_NOT_FOUND) std::cout<<"ERROR_PATH_NOT_FOUND  : "<<dirname.c_str()<<std::endl;//debug
			if(!Directory_Exists(dirname))
			{				
				std::cerr<<dirname.c_str()<<" is NOT created Successfully."<<std::endl;
				return false; 
			}
			std::cout<<dirname.c_str()<<" is created Successfully."<<std::endl;
		}
		return true;
	}
};

