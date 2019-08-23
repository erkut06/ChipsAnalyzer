#include "ImageLoader.h"
#include <vector>
#include <vtkBMPReader.h>
#include "vtkImageData.h"
#include <algorithm>
#include "ThreadPool.h"
#include "vtkExtractVOI.h"
#include "vtkSmartPointer.h"
#ifdef _WIN32
#include <Windows.h>
#else

#include <sys/types.h>
#include <dirent.h>

#endif


std::vector<std::string> GetAllFileNamesFromDirectory(std::string folder)
{
	std::vector<std::string> names;
	std::string search_path = folder + "/*.bmp";
#ifdef _WIN32
	WIN32_FIND_DATA fd;
	HANDLE hFind = ::FindFirstFile(search_path.c_str(), &fd);
	if (hFind != INVALID_HANDLE_VALUE) {
		do {
			
			if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
				size_t found = std::string(fd.cFileName).find("spr");
				if (found == std::string::npos)
				{
					names.push_back(folder + "\\" + fd.cFileName);
				}
			}
		} while (::FindNextFile(hFind, &fd));
		::FindClose(hFind);
	}
#else
	DIR* dirFile = opendir(folder.c_str());
	if (dirFile)
	{
		struct dirent* hFile;
		errno = 0;
		while ((hFile = readdir(dirFile)) != NULL)
		{
			if (!strcmp(hFile->d_name, ".")) continue;
			if (!strcmp(hFile->d_name, "..")) continue;
			// dirFile.name is the name of the file. Do whatever string comparison
			// you want here. Something like:
			if (strstr(hFile->d_name, ".bmp") && !strstr(hFile->d_name, "spr"))
				names.push_back(folder + "\/" + hFile->d_name);
		}
		closedir(dirFile);
	}
#endif
	return names;
}

ImageLoader::ImageLoader()
{
}

ImageLoader::~ImageLoader()
{
	folderName.clear();
	if (data != NULL)
	{
		data->Delete();
	}
}

void ImageLoaderThreadMethod(const char* fileName, vtkImageData* destination, int sliceNumber, int sliceSize)
{
	vtkBMPReader* reader = vtkBMPReader::New();
	reader->SetFileName(fileName);
	reader->SetAllow8BitBMP(1);
	reader->Update();
	vtkImageData* result = reader->GetOutput();
	unsigned char* dPtr = (unsigned char*)destination->GetScalarPointer();
	dPtr += static_cast<size_t>((static_cast<size_t>(sliceSize) * static_cast<size_t>(sliceNumber)));
	unsigned char* cPtr = (unsigned char*)reader->GetOutput()->GetScalarPointer();
	memcpy(dPtr, cPtr, sliceSize);
	reader->Delete();
}

void ImageLoader::Read()
{
	std::vector<std::string> fileList;
	fileList = GetAllFileNamesFromDirectory(folderName);
	std::sort(fileList.begin(), fileList.end());
	//read the first one
	vtkBMPReader* reader = vtkBMPReader::New();
	reader->SetAllow8BitBMP(1);
	reader->SetFileName(fileList[0].c_str());
	reader->Update();

	int* dims = reader->GetOutput()->GetDimensions();
	int scalarType = reader->GetOutput()->GetScalarType();
	int numberOfScalarComponents = reader->GetOutput()->GetNumberOfScalarComponents();
	int scalarSize = reader->GetOutput()->GetScalarSize();
	size_t numberOfBytes = dims[0] * dims[1] * numberOfScalarComponents * scalarSize;
	numberOfBytes *= fileList.size();

	if (data != NULL)
	{
		data->Delete();
		data = NULL;
	}
	data = vtkImageData::New();
	data->SetDimensions(dims[0], dims[1], static_cast<int>(fileList.size()));
	data->SetScalarType(scalarType, data->GetInformation());
	data->SetNumberOfScalarComponents(1, data->GetInformation());
	data->AllocateScalars(data->GetInformation());
	unsigned char* dPtr = (unsigned char*)data->GetScalarPointer();
	int sliceSize = dims[0] * dims[1];

	ThreadPool pool(8);
	std::vector< std::future<void> > results;

	for (int i = 0; i < fileList.size(); i++)
	{
		std::string currentFileName = fileList[i];
		results.emplace_back(
			pool.enqueue([currentFileName, this, i, sliceSize]() {
			ImageLoaderThreadMethod(currentFileName.c_str(), data, i, sliceSize);
		})
		);
	}

	for (auto && result : results)
		result.get();

	reader->Delete();
	dims = data->GetDimensions();
	double yScaleUp = 30.0 / 1144.0;
	double yScaleDown = 65.0 / 1144.0;
	double xScale = 50 / 1144.0;

	int xMin = xScale * static_cast<double>(dims[0]);
	int xMax = dims[0] - xMin;

	int yMin = yScaleUp * static_cast<double>(dims[1]);
	int yMax = dims[1] - yScaleDown * static_cast<double>(dims[1]);

	vtkSmartPointer<vtkExtractVOI> extractVoi = vtkSmartPointer<vtkExtractVOI>::New();
	extractVoi->SetInputData(data);
	extractVoi->SetVOI(xMin, xMax, yMin, yMax, 0, dims[2]);
	extractVoi->Update();

	data->ShallowCopy(extractVoi->GetOutput());

}

void ImageLoader::SetFolderName(std::string folder)
{
	folderName = folder;
}

vtkImageData* ImageLoader::GetOutput()
{
	return data;
}
