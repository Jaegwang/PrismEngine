/********************************************************************
	CAPTURE_MANAGER.cpp
*********************************************************************/

#include "CAPTURE_MANAGER.h"
#include "READ_WRITE_PNG.h"
#include "FILE_UTILITIES.h"
#include <string>

CAPTURE_MANAGER::CAPTURE_MANAGER()
	: last_frame_(-1)
	, auto_run_(false)
	, auto_copy_script_(false)
	, auto_capture_(false)
	, auto_video_(false)
	, auto_delete_image_(false)
	, auto_exit_(false)
	, output_common_path_(std::string())
	, output_image_path_(std::string())
	, output_video_path_(std::string())
	, output_log_path_(std::string())
	, output_script_path_(std::string())
	, output_common_file_basename_(std::string())
	, output_image_file_basename_(std::string())
	, output_video_file_basename_(std::string())
	, output_log_file_basename_(std::string())
	, output_script_file_basename_(std::string())
{
}

CAPTURE_MANAGER::~CAPTURE_MANAGER()
{
}

void CAPTURE_MANAGER::Initialize(std::string& script_abs_path)
{
	char buffer[MAX_PATH];
	absolute_path_ = _fullpath(buffer, 0, MAX_PATH);

	output_image_path_ = absolute_path_;
	output_image_path_.append("\\images");	
	output_image_file_basename_ = "image";

	output_video_path_ = absolute_path_;
	output_video_path_.append("\\video");
	output_video_file_basename_ = "video";

	//make folders for image, video, log, script
	FILE_UTILITIES file_utilities;
	if (!output_image_path_.empty())
	if (!file_utilities.Directory_Exists(output_image_path_))
		file_utilities.Create_Directory(output_image_path_);

	if (!output_video_path_.empty())
	if (!file_utilities.Directory_Exists(output_video_path_))
		file_utilities.Create_Directory(output_video_path_);
}

void CAPTURE_MANAGER::ResetCapture()
{
//	Initialize(PROJECT_INFO::GetScriptAbsPath());
}

void CAPTURE_MANAGER::CopyScript()
{	
	/*
	std::string output_script_file_basename = output_script_file_basename_;
	output_script_file_basename += "_" + PROJECT_INFO::GetBeginDateTime();

	std::string script_target;		
	script_target.append(output_script_path_);
	script_target.append(output_script_file_basename);
	script_target.append(".sim");

	if(CopyFile(PROJECT_INFO::GetScriptAbsPath().c_str(), script_target.c_str(), false))
		std::cout << "Copying script file to \"" << script_target << "\" has been succeed." << std::endl;		
	else
		std::cout << "Copying script file to \"" << script_target << "\" has been failed." << std::endl;
	*/
}

void CAPTURE_MANAGER::AutoCopyScript()
{
	if(auto_copy_script_)
		CopyScript();
}

void CAPTURE_MANAGER::CaptureImage(int current_frame, int width, int height)
{
	//NOTE : This causes images are not saved when opengl window is minimized.
	if(last_frame_ > 0 && current_frame <= last_frame_)
	{
		char image_saving_pos[MAX_PATH] = {0,};
		sprintf(image_saving_pos, "%s\\%s_%04d.png", 
			output_image_path_.c_str(),
			output_image_file_basename_.c_str(), 
			current_frame);
		READ_WRITE_PNG::PngSaveImage(image_saving_pos, width, height);	
	}else if(last_frame_ < 0)
	{
		char image_saving_pos[MAX_PATH] = {0,};
		sprintf(image_saving_pos, "%s\\%s_%04d.png", 
			output_image_path_.c_str(),
			output_image_file_basename_.c_str(), 
			current_frame);
		READ_WRITE_PNG::PngSaveImage(image_saving_pos, width, height);	
	}
}

void CAPTURE_MANAGER::MakeVideoAtLastFrame()
{
	if(auto_video_) 
		MakeVideo();

	if(auto_exit_) 
		exit(0);
}

void CAPTURE_MANAGER::MakeVideo()
{
//	using namespace boost::posix_time;
//	static std::locale loc(std::cout.getloc(), new time_facet("%Y%m%d_%H%M%S"));

	FILE_UTILITIES file_utilities;
	std::string video_encoding_command, removing_image_target;
	std::string video_file_basename = output_video_file_basename_;

//	video_file_basename += "_" + PROJECT_INFO::GetBeginDateTime();
	
	video_encoding_command = absolute_path_;
	video_encoding_command.append("\\ffmpeg -i \"");	
	video_encoding_command.append(output_image_path_);	
	video_encoding_command.append("\\");
	video_encoding_command.append(output_image_file_basename_);
	video_encoding_command.append("_%04d.png\" -r 30 -g 1 -qscale 1 \"");	
	video_encoding_command.append(output_video_path_);		
	video_encoding_command.append("\\");
	video_encoding_command.append(video_file_basename); 
	video_encoding_command.append(".avi\"");		
	file_utilities.Execute_Process(video_encoding_command);		

	std::cout<<"-end of video encoding"<<std::endl;// debug

	//removing images in images folder
	if(auto_delete_image_ == true)
	{
		char ch[100];
		for(int i=1; i <= last_frame_;i++)
		{
			removing_image_target = output_image_path_;
			removing_image_target.append(output_image_file_basename_);
			sprintf(ch, "_%04d.png", i);
			removing_image_target.append(ch);

			remove(removing_image_target.c_str());
		}
		std::cout<<"-end of deleting images"<<std::endl;// debug
	}
}

void CAPTURE_MANAGER::SetBatchMode()
{
	auto_run_ = true;
	auto_capture_ = true;
	auto_video_ = true;
	auto_delete_image_ = true;
	auto_exit_ = true;
}