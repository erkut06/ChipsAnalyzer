#include "ChipsSegmentation.h"
#include <vtkImageData.h>
#include <vtkJPEGWriter.h>
#include <vtkSmartPointer.h>
#include <vtkImageFlip.h>
#include <string>
#include <vtkImageDilateErode3D.h>
#include "ThreadPool.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkVTKImageToImageFilter.h"
#include <itkImageToVTKImageFilter.h>
#include "itkConnectedComponentImageFilter.h"
#include "itkMesh.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkMeshFileWriter.h"
#include "itkLineCell.h"
#include "itkVTKPolyDataWriter.h"
#include <vtkMarchingCubes.h>
#include "vtkPolyDataWriter.h"
#include "DataStructures.h"
#include "Vector3.h"
#include "BoundaryModel2D.h"
#include "vtkMetaImageWriter.h"
#include <vector>
#include "QFileDialog.h"
#include "QCoreApplication.h"
#include "vtkBMPWriter.h"
#include "vtkTIFFWriter.h"
//#define WRITE_CHIPS_MASK_AS_MHD

typedef   unsigned short         InternalPixelType;
typedef itk::Image< InternalPixelType, 3 >  InternalImageType;
typedef unsigned short                         OutputPixelType;
typedef itk::Image< OutputPixelType, 3 > OutputImageType;
typedef itk::CastImageFilter< InternalImageType, OutputImageType >
CastingFilterType;
typedef itk::ConnectedThresholdImageFilter< InternalImageType,
	InternalImageType > ConnectedFilterType;

ChipsSegmentation::ChipsSegmentation()
{

}

ChipsSegmentation::~ChipsSegmentation()
{
	chipsMask->Delete();
	connectedMask->Delete();
}

void ChipsSegmentation::ExtractSliceSurface(int sliceNum, double maxRadius, unsigned short forbiddenFruit)
{
	int* dims = chipsMask->GetDimensions();
	BoundaryModel2D* model = new BoundaryModel2D();
	model->SetCenter(static_cast<double>(dims[0] / 2), static_cast<double>(dims[1] / 2));
	model->AllocatePoints(5, maxRadius);
	size_t sliceOffSet = static_cast<size_t>(sliceNum) * static_cast<size_t>(dims[0]) * static_cast<size_t>(dims[1]);
	unsigned short* dPtr = &((unsigned short*)chipsMask->GetScalarPointer())[sliceOffSet];
	model->FitBoundary(dPtr, dims, forbiddenFruit);
	delete model;
}

void ChipsSegmentation::ExtractSurface()
{

	CalculateChipsMask();
	unsigned short uShortMax = std::numeric_limits<unsigned short>::max();
	int* dims = this->inputData->GetDimensions();
	double maxRadius = sqrt(dims[0] * dims[0] + dims[1] * dims[1] + dims[2] * dims[2]);
	ThreadPool pool(8);
	std::vector< std::future<void> > results;

	for (int slice = 0; slice< dims[2]; slice++)
	{
		results.emplace_back(
			pool.enqueue([this,slice,maxRadius,uShortMax]() {
			ExtractSliceSurface(slice, maxRadius, uShortMax);
		}));
		
	}
	for (auto && result : results)
		result.get();

	CalculateInitialPorosity();

#ifdef WRITE_CHIPS_MASK_AS_MHD

	vtkMetaImageWriter* writer = vtkMetaImageWriter::New();
	writer->SetInputData(chipsMask);
	writer->SetFileName("C:\\Test\\chipsMask.mhd");
	writer->SetFileDimensionality(3);
	writer->Update();
	writer->Write();
	std::cout << "chipsMask written!!!" << std::endl;
#endif

}

void ChipsSegmentation::CalculateInitialPorosity()
{

	unsigned short* chipsPtr = (unsigned short*)chipsMask->GetScalarPointer();
	unsigned short maskValue = std::numeric_limits<unsigned short>::max();
	unsigned short borderValue = maskValue / 2;
	unsigned short outValue = maskValue / 4;
	unsigned short void_value = 0;

	int* dims = chipsMask->GetDimensions();
	size_t maxSize = static_cast<size_t>(dims[0]) * static_cast<size_t>(dims[1]) * static_cast<size_t>(dims[2]);
	size_t outsideCount = 0;
	size_t voidCount = 0;
	size_t chipsCount = 0;

	for (int i = 0; i < maxSize; i++)
	{
		if (chipsPtr[i] == void_value)
		{
			voidCount++;
		}
		if (chipsPtr[i] == borderValue || chipsPtr[i] == outValue)
		{
			outsideCount++;
		}

		if (chipsPtr[i] == maskValue)
		{
			chipsCount++;
		}
	}

	size_t _32bitMax = static_cast<size_t>(1024) * static_cast < size_t>(1024) * static_cast < size_t>(1024) * static_cast < size_t>(2) - static_cast < size_t>(1);

	size_t _maxDivident = maxSize / _32bitMax + 1;

	maxSize /= _maxDivident;
	outsideCount /= _maxDivident;
	voidCount /= _maxDivident;
	chipsCount /= _maxDivident;
	outRatio = static_cast<double>(outsideCount) / static_cast<double>(maxSize);
	initialPorosity = static_cast<double>(voidCount) / static_cast<double>(maxSize - outsideCount);
	initialPorosity *= 100.0;
}

void ChipsSegmentation::FillChipsMask()
{
	unsigned short* chipsPtr = (unsigned short*)chipsMask->GetScalarPointer();
	unsigned char* iPtr = (unsigned char*)inputData->GetScalarPointer();
	unsigned short* mPtr = (unsigned short*)chipsPtr;
	int* dims = chipsMask->GetDimensions();
	int sliceSize = dims[0] * dims[1];

	unsigned short maskValue = std::numeric_limits<unsigned short>::max();
	unsigned short borderValue = maskValue / 2;
	unsigned short outsideValue = maskValue / 4;
	for (int z = 0; z < dims[2]; z++)
	{
		unsigned short* maskSlice = &mPtr[z * sliceSize];

		for (int y = 0; y < dims[1]; y++)
		{
			//left check
			for (int x = 0; x < dims[0]; x++)
			{
				size_t sliceIndex = static_cast<size_t>(y) * static_cast<size_t>(dims[0]) + static_cast<size_t>(x);

				if (maskSlice[sliceIndex] == borderValue)
				{
					break;
				}
				else
				{
					maskSlice[sliceIndex] = outsideValue;
				}
			}

			//right check

			for (int x = dims[0] - 1; x >= 0; x--)
			{
				size_t sliceIndex = static_cast<size_t>(y) * static_cast<size_t>(dims[0]) + static_cast<size_t>(x);
				if (maskSlice[sliceIndex] == borderValue)
				{
					break;
				}
				else
				{
					maskSlice[sliceIndex] = outsideValue;
				}

			}

		}

		for (int x = 0; x < dims[0]; x++)
		{
			//up check
			for (int y = 0; y < dims[1]; y++)
			{
				size_t sliceIndex = static_cast<size_t>(y) * static_cast<size_t>(dims[0]) + static_cast<size_t>(x);
				if (maskSlice[sliceIndex] == borderValue)
				{
					break;
				}
				else
				{
					maskSlice[sliceIndex] = outsideValue;
				}
			}
			//down check
			for (int y = dims[1] - 1; y >= 0; y--)
			{
				size_t sliceIndex = static_cast<size_t>(y) * static_cast<size_t>(dims[0]) + static_cast<size_t>(x);
				if (maskSlice[sliceIndex] == borderValue)
				{
					break;
				}
				else
				{
					maskSlice[sliceIndex] = outsideValue;
				}
			}

		}
	}

}

void ChipsSegmentation::CalculateChipsMask()
{

	int* dims = inputData->GetDimensions();
	if (chipsMask != NULL)
	{
		chipsMask->Delete();
	}
	
	chipsMask = vtkImageData::New();

	chipsMask->SetScalarType(5, chipsMask->GetInformation());
	chipsMask->SetNumberOfScalarComponents(1, chipsMask->GetInformation());
	chipsMask->SetDimensions(dims);
	chipsMask->AllocateScalars(chipsMask->GetInformation());
	unsigned char* iPtr = (unsigned char*)inputData->GetScalarPointer();
	unsigned short* cPtr = (unsigned short*)chipsMask->GetScalarPointer();
	unsigned short uShortMax = std::numeric_limits<unsigned short>::max();
	unsigned short uShortMin = std::numeric_limits<unsigned short>::min();
	size_t volumeSize = static_cast<size_t>(dims[0]) * static_cast<size_t>(dims[1]) * static_cast<size_t>(dims[2]);
	memset(cPtr, static_cast < unsigned short>(0), volumeSize * sizeof(unsigned short));

	unsigned short tMin = static_cast<unsigned short>(chipsMin);
	unsigned short tMax = static_cast<unsigned short>(chipsMax);

	for (size_t i = 0; i <volumeSize; i++)
	{
		if (iPtr[i] >= tMin && iPtr[i] <= tMax)
		{
			cPtr[i] = uShortMax;
		}
		else
		{
			cPtr[i] = uShortMin;
		}
	}
}

void WriteToJpeg(const char* path, vtkImageData* image)
{
	vtkJPEGWriter* writer = vtkJPEGWriter::New();
	//vtkSmartPointer<vtkImageFlip> flipYFilter =
	//	vtkSmartPointer<vtkImageFlip>::New();
	//flipYFilter->SetFilteredAxis(1); // flip y axis
	//flipYFilter->SetInputData(image);
	//flipYFilter->Update();
	writer->SetInputData(image);
	writer->SetFileName(path);
	writer->Update();
	writer->Write();
	writer->Delete();
}

void GetPointFromIndex(size_t startIndex, int* dims, Point3DInt& p)
{
	size_t sliceSize = static_cast<size_t>(dims[0]) * static_cast<size_t>(dims[1]);
	p.z = startIndex / sliceSize;
	p.y = (startIndex - p.z * sliceSize) / static_cast<size_t>(dims[0]);
	p.x = startIndex - p.z * sliceSize - p.y * static_cast<size_t>(dims[0]);
}

void GrowItk(OutputImageType::ConstPointer myitkImage,vtkImageData* data, Point3DInt& point, unsigned short label, int neighbourHood)
{

	unsigned short ushortMax = std::numeric_limits<unsigned short>::max();
	const     unsigned int    Dimension = 3;
	CastingFilterType::Pointer caster = CastingFilterType::New();
	ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
	//myitkImage->Print(std::cout);
	connectedThreshold->SetInput(myitkImage);
	connectedThreshold->SetLower(ushortMax);
	connectedThreshold->SetUpper(ushortMax);
	connectedThreshold->SetReplaceValue(label);

	if (neighbourHood == 1)
	{
		connectedThreshold->SetConnectivity(ConnectedFilterType::ConnectivityEnumType::FullConnectivity);
	}
	else
	{
		connectedThreshold->SetConnectivity(ConnectedFilterType::ConnectivityEnumType::FaceConnectivity);
	}

	InternalImageType::IndexType  index;
	index[0] = point.x;
	index[1] = point.y;
	index[2] = point.z;
	connectedThreshold->SetSeed(index);
	connectedThreshold->SetNumberOfThreads(8);
	connectedThreshold->Update();
	typedef itk::ImageToVTKImageFilter<OutputImageType>       ConnectorType;
	ConnectorType::Pointer connector = ConnectorType::New();
	connector->SetInput(connectedThreshold->GetOutput());
	connector->Update();
	vtkImageData* result = connector->GetOutput();
	unsigned short* rPtr = (unsigned short*)result->GetScalarPointer();
	int* dims = data->GetDimensions();
	size_t volumeSize = static_cast<size_t>(dims[0]) * static_cast<size_t>(dims[1]) * static_cast<size_t>(dims[2]);
	unsigned short* dPtr = (unsigned short*)data->GetScalarPointer();
	for (int i = 0; i < volumeSize ; i++)
	{
		if (rPtr[i] == label)
		{
			dPtr[i] = rPtr[i];
		}
	}
}

vtkImageData* ChipsSegmentation::GetMask()
{
	return this->chipsMask;
}

void ChipsSegmentation::LabelAndAnalyze()
{
	int* dims = chipsMask->GetDimensions();
	unsigned short* mPtr = (unsigned short*)chipsMask->GetScalarPointer();
	size_t volumeSize = static_cast<size_t>(dims[0]) * static_cast<size_t>(dims[1]) * static_cast<size_t>(dims[2]);
	connectedMask = vtkImageData::New();
	connectedMask->SetScalarType(5, connectedMask->GetInformation());
	connectedMask->SetNumberOfScalarComponents(1, chipsMask->GetInformation());
	connectedMask->SetDimensions(dims);
	connectedMask->AllocateScalars(connectedMask->GetInformation());
	unsigned short* connectedMap = (unsigned short*)connectedMask->GetScalarPointer();
	for (int i = 0; i < volumeSize; i++)
	{
		connectedMap[i] = 0;
	}

	for (int i = 0; i < volumeSize; i++)
	{
		if (mPtr[i] < 1)
		{
			connectedMap[i] = std::numeric_limits<unsigned short>::max();
		}
	}

	unsigned short componentLabel = 0;


	unsigned short ushortMax = std::numeric_limits<unsigned short>::max();

	using FilterType = itk::VTKImageToImageFilter< OutputImageType >;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput(connectedMask);
	filter->SetNumberOfThreads(8);
	filter->Update();
	OutputImageType::ConstPointer myitkImage = filter->GetOutput();
	while (true)
	{
		//check if any node left
		bool carryOn = false;
		size_t startIndex = 0;

		for (size_t i = 0; i < volumeSize; i++)
		{
			if (connectedMap[i] == ushortMax)
			{
				carryOn = true;
				startIndex = i;
				break;
			}
		}
		if (carryOn == false)
		{
			break;
		}
		componentLabel++;
		Point3DInt p;
		GetPointFromIndex(startIndex, dims, p);
		GrowItk(myitkImage,connectedMask, p, static_cast<unsigned short>(componentLabel), neighbourHood);
		//GrowOnMask(ushortMax, connectedMap, startIndex, dims, componentLabel);
	}

	std::cout << componentLabel << " void component found"<<std::endl;

	double outsideVoidCount = 0;
	for (int i = 0; i < volumeSize; i++)
	{
		if (mPtr[i] == 0)
		{
			outsideVoidCount += 1.0;
		}
	}

	double chipsVolume = (double)volumeSize - outsideVoidCount;

	std::vector<double> voidComponentsVolumeList;
	double biggestVoidVolume = std::numeric_limits<double>::min();
	int biggestVoidIndex = 1;
	for (unsigned short i = 1; i <= componentLabel; i++)
	{
		double labelSize = 0;
		for (int ii = 0; ii < volumeSize; ii++)
		{
			if (connectedMap[ii] == i)
			{
				labelSize++;
			}
		}
		voidComponentsVolumeList.push_back(labelSize);
		if (biggestVoidVolume < labelSize)
		{
			biggestVoidVolume = labelSize;
			biggestVoidIndex = i;
		}
	}

	double totalVoidVolume = 0;

	for (int i = 0; i < voidComponentsVolumeList.size(); i++)
	{
		totalVoidVolume += voidComponentsVolumeList[i];
	}

	connectivityIndex = (biggestVoidVolume / totalVoidVolume);
	//std::cout << "connectivity index is " << connectivityIndex << std::endl;

	unsigned short mult = (ushortMax-1) / componentLabel;

	for (int i = 0; i < volumeSize; i++)
	{
		if (connectedMap[i] == biggestVoidIndex)
		{
			connectedMap[i] = ushortMax;
		}
		else
		{
			connectedMap[i] *= mult;
		}
	}

	/*vtkMetaImageWriter* writer = vtkMetaImageWriter::New();
	writer->SetInputData(connectedMask);
	writer->SetFileName("C:\\Test\\connectedMap.mhd");
	writer->SetFileDimensionality(3);
	writer->Update();
	writer->Write();
	std::cout << "connectedMap written!!!" << std::endl;*/
}

int GetNumberOfDigits2(int num)
{
	std::string numStr = std::to_string(num);
	return static_cast<int>(numStr.size());
}

void ChipsSegmentation::WriteChipsBoundary(std::string currentFolderName)
{
	QString appDir = QCoreApplication::applicationDirPath();
	std::string folderName = appDir.toUtf8().constData();
	folderName += "/" + currentFolderName;
	QDir dir(folderName.c_str());
	if (!dir.exists())
		dir.mkpath(".");
	folderName += "/Boundary/";
	QDir dir2(folderName.c_str());
	if (!dir2.exists())
		dir2.mkpath(".");

	/*foreach(QString dirFile, dir.entryList())
	{
		dir.remove(dirFile);
	}*/

	unsigned short* chipsPtr = (unsigned short*)chipsMask->GetScalarPointer();
	unsigned short maskValue = std::numeric_limits<unsigned short>::max();
	unsigned short borderValue = maskValue / 2;
	unsigned short outValue = maskValue / 4;
	unsigned short void_value = 0;

	int* dims = chipsMask->GetDimensions();
	size_t maxSize = static_cast<size_t>(dims[0]) * static_cast<size_t>(dims[1]) * static_cast<size_t>(dims[2]);
	size_t outsideCount = 0;
	size_t voidCount = 0;
	size_t chipsCount = 0;

	vtkImageData* currentSlice = vtkImageData::New();
	currentSlice->SetScalarType(VTK_UNSIGNED_CHAR, currentSlice->GetInformation());
	currentSlice->SetDimensions(dims[0], dims[1], 1);
	currentSlice->SetNumberOfScalarComponents(1, currentSlice->GetInformation());
	currentSlice->AllocateScalars(currentSlice->GetInformation());
	vtkBMPWriter* writer = vtkBMPWriter::New();

	size_t volumeIndex = 0;
	for (size_t z = 0; z < dims[2]; z++)
	{
		size_t imageIndex = 0;
		unsigned char* cPtr = (unsigned char*)currentSlice->GetScalarPointer();
		memset(cPtr, 0, dims[0] * dims[1]);
		for (size_t y = 0; y < dims[1] ; y++)
		{
			for (size_t x = 0 ; x < dims[0] ; x++)
			{
				if (chipsPtr[volumeIndex] == borderValue || chipsPtr[volumeIndex] == outValue)
				{
					cPtr[imageIndex] = 128;
				}

				if (chipsPtr[volumeIndex] == maskValue)
				{
					cPtr[imageIndex] = 255;
				}
				imageIndex++;
				volumeIndex++;
			}
		}
		writer->SetInputData(currentSlice);
		std::string filename = folderName;
		int zeroCount = 5 - GetNumberOfDigits2(z);
		for (int k = 0; k < zeroCount; k++)
		{
			filename += "0";
		}
		filename += std::to_string(z);
		filename += ".bmp";
		writer->SetFileName(filename.c_str());
		writer->Write();

	}
	currentSlice->Delete();
	writer->Delete();
}

void ChipsSegmentation::WriteConnectedComponents(std::string currentFolderName)
{
	QString appDir = QCoreApplication::applicationDirPath();
	std::string folderName = appDir.toUtf8().constData();
	folderName += "/" + currentFolderName;
	QDir dir(folderName.c_str());
	if (!dir.exists())
		dir.mkpath(".");
	folderName += "/ConnectedComponents/";
	QDir dir2(folderName.c_str());
	if (!dir2.exists())
		dir2.mkpath(".");

	/*foreach(QString dirFile, dir.entryList())
	{
		dir.remove(dirFile);
	}*/

	unsigned short* maskPtr = (unsigned short*)connectedMask->GetScalarPointer();
	unsigned short maskValue = std::numeric_limits<unsigned short>::max();
	unsigned short borderValue = maskValue / 2;
	unsigned short outValue = maskValue / 4;
	unsigned short void_value = 0;

	int* dims = connectedMask->GetDimensions();
	size_t maxSize = static_cast<size_t>(dims[0]) * static_cast<size_t>(dims[1]) * static_cast<size_t>(dims[2]);
	size_t outsideCount = 0;
	size_t voidCount = 0;
	size_t chipsCount = 0;

	vtkImageData* currentSlice = vtkImageData::New();
	currentSlice->SetScalarType(VTK_UNSIGNED_SHORT, currentSlice->GetInformation());
	currentSlice->SetDimensions(dims[0], dims[1], 1);
	currentSlice->SetNumberOfScalarComponents(1, currentSlice->GetInformation());
	currentSlice->AllocateScalars(currentSlice->GetInformation());
	vtkTIFFWriter* writer = vtkTIFFWriter::New();

	size_t volumeIndex = 0;
	for (size_t z = 0; z < dims[2]; z++)
	{
		size_t imageIndex = 0;
		unsigned short* cPtr = (unsigned short*)currentSlice->GetScalarPointer();
		memcpy(cPtr, &maskPtr[volumeIndex], dims[0] * dims[1] * 2);
		
		writer->SetInputData(currentSlice);
		std::string filename = folderName;
		int zeroCount = 5 - GetNumberOfDigits2(z);
		for (int k = 0; k < zeroCount; k++)
		{
			filename += "0";
		}
		filename += std::to_string(z);
		filename += ".bmp";
		writer->SetFileName(filename.c_str());
		writer->Write();

		volumeIndex += dims[0] * dims[1];

	}
	currentSlice->Delete();
	writer->Delete();
}

void ChipsSegmentation::SetInputData(vtkImageData* i)
{
	this->inputData = i;
}

void ChipsSegmentation::SetAirIntensities(unsigned char min, unsigned char max)
{
	this->airMin = min;
	this->airMax = max;
}

void ChipsSegmentation::SetChipsIntensities(unsigned char cMn, unsigned char cMx)
{
	this->chipsMin = cMn;
	this->chipsMax = cMx;
}