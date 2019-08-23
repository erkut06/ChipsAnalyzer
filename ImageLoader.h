#pragma once
#include <string>

class vtkImageData;
class ImageLoader
{
public:
	ImageLoader();
	~ImageLoader();
	void SetFolderName(std::string);
	void Read();
	vtkImageData* GetOutput();
private:
	std::string folderName;
	vtkImageData* data = NULL;
};

