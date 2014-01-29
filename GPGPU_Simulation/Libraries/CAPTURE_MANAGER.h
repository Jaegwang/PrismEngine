/********************************************************************
	CAPTURE_MANAGER.h
*********************************************************************/

#pragma once

#include <string>

class CAPTURE_MANAGER
{
public:
    CAPTURE_MANAGER();
	~CAPTURE_MANAGER();

	void Initialize(std::string& script_abs_path);
	void ResetCapture();

	void CopyScript();
	void AutoCopyScript();

	bool IsAutoCaptureImage() { return auto_capture_; }
	void CaptureImage(int current_frame, int width, int height);

	bool IsAutoCaptureMovieAtLastFrame() { return auto_video_; }
	bool IsAutoDeleteImage() { return auto_delete_image_; }

	bool IsAutoRun() { return auto_run_;}

	void MakeVideo();
	void MakeVideoAtLastFrame();

	void SetBatchMode();

private:
	int last_frame_;
	bool auto_run_;
	bool auto_copy_script_;
	bool auto_capture_;
	bool auto_video_;
	bool auto_delete_image_;
	bool auto_exit_;

	std::string absolute_path_;
	std::string output_common_path_;
	std::string output_image_path_;
	std::string output_video_path_;
	std::string output_log_path_;
	std::string output_script_path_; 

	std::string output_common_file_basename_;
	std::string output_image_file_basename_;
	std::string output_video_file_basename_;
	std::string output_log_file_basename_;
	std::string output_script_file_basename_;
};

