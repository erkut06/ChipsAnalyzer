#include "ChipsAnalyzer.h"
#include "ui_ChipsAnalyzer.h"
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkCubeSource.h>
#include "vtkImageData.h"
#include "vtkOpenGLGPUVolumeRayCastMapper.h"
#include "vtkVolume.h"
#include "vtkColorTransferFunction.h"
#include "vtkPiecewiseFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCamera.h"
#include "vtkSmartVolumeMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkFileOutputWindow.h"
#include "vtkSphereSource.h"
#include "vtkMetaImageReader.h"
#include "vtkExtractVOI.h"
#include "vtkImageHistogramStatistics.h"
#include "ThreadPool.h"
#include <future>
#include "vtkIdTypeArray.h"
#include "vtkImageResample.h"
#include "ImageLoader.h"
#include "QFileDialog.h"
#include "QCheckBox.h"
#include "ChipsSegmentation.h"
#include "vtkBMPWriter.h"
#include <chrono>

void ChipsAnalyzer::PreProcessData(vtkImageData* data)
{

	int* dims = data->GetDimensions();
	double dataSize = static_cast<double>(dims[0]) * static_cast<double>(dims[1]) * static_cast<double>(dims[2]);
	size_t dataSizeInt = static_cast<size_t>(dims[0]) * static_cast<size_t>(dims[1]) * static_cast<size_t>(dims[2]);
	//gamma correction
	double gamma = ui->gammaCorrection->value();

	if (gamma != 1.0)
	{
		double* gammaCorrected = new double[dataSizeInt];
		double exponent = 1.0 / gamma;
		unsigned char* d = (unsigned char*)data->GetScalarPointer();
		double min = std::numeric_limits<double>::max();
		double max = std::numeric_limits<double>::min();

		for (size_t currentIndex = 0; currentIndex < dataSizeInt; currentIndex++)
		{
			double currentValue = static_cast<double>(d[currentIndex]);
			gammaCorrected[currentIndex] = std::pow(currentValue, exponent);

			if (min > gammaCorrected[currentIndex])
			{
				min = gammaCorrected[currentIndex];
			}

			if (max < gammaCorrected[currentIndex])
			{
				max = gammaCorrected[currentIndex];
			}

		}

		double multiplier = 255 / max;

		for (size_t currentIndex = 0; currentIndex < dataSizeInt; currentIndex++)
		{
			gammaCorrected[currentIndex] *= multiplier;
			d[currentIndex] = static_cast<unsigned char>(gammaCorrected[currentIndex]);
		}
		delete[] gammaCorrected;
	}

	std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
	std::cout << "Histogram analysis started" << std::endl;
	stats = vtkImageHistogramStatistics::New();
	stats->SetInputData(data);
	stats->SetNumberOfThreads(8);
	stats->SetAutoRangePercentiles(0, 100);
	stats->Update();

	std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
	long long duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

	std::cout << "Histogram analysis finished in " << duration << " ms" << std::endl;
	vtkIdTypeArray* his = stats->GetHistogram();
	vtkIdType size = his->GetSize();
	vtkIdType* hisPtr = his->GetPointer(0);
	
	std::vector<double> changeRatios;
	double* range = stats->GetAutoRange();
	std::cout << "Data range is " << range[0] << " " << range[1] << std::endl;
	double mean = stats->GetMean();
	std::cout << "Data mean is " << mean << std::endl;
	double median = stats->GetMedian();
	std::cout << "Data median is " << median << std::endl;

	int maxIndex = 0;
	vtkIdType maxHis = std::numeric_limits<long long>::min();
	/*std::cout << "Histogram : \n";

	for (int i = 0; i < size; i++)
	{
		std::cout << hisPtr[i] << " %" << (static_cast<double>(hisPtr[i]) / dataSize) * 100 << std::endl;
		if (maxHis < hisPtr[i])
		{
			maxHis = hisPtr[i];
			maxIndex = i;
		}
	}

	std::cout << "Max histogram intensity is " << maxIndex << std::endl;*/

	for (int i = maxIndex + 1; i < size - 1; i++)
	{
		double current = static_cast<double>(hisPtr[i]);
		double next = static_cast<double>(hisPtr[i + 1]);
		double changeRatio = 0;
		if (current == 0)
		{
			current = std::numeric_limits<double>::min();
		}
		changeRatio = (next - current) / current;
		changeRatios.push_back((changeRatio));
	}

	airIntensities[0] = 0;
	airIntensities[1] = maxIndex;
	chipsIntensities[0] = maxIndex + 1;

	double firstChange = changeRatios[maxIndex];
	double noiseChangeRate = firstChange / 10.0;

	for (int i = maxIndex; i < changeRatios.size(); i++)
	{
		if (std::abs(changeRatios[i]) < std::abs(noiseChangeRate))
		{
			airIntensities[1] = i;
			chipsIntensities[0] = i + 1;
			break;
		}
	}

	for (int i = 0; i < changeRatios.size(); i++)
	{
		if (std::abs(changeRatios[i]) > 0.5)
		{
			break;
		}
		chipsIntensities[1] = i + 1;
	}
	chipsIntensities[1] = range[1];
	std::cout << "Air & noise intensities calculated between " << airIntensities[0] << " " << airIntensities[1] << std::endl;
	std::cout << "Chips intensities calculated between " << chipsIntensities[0] << " " << chipsIntensities[1] << std::endl;
	size_t volumeSize = static_cast<size_t>(dims[0]) * static_cast<size_t>(dims[1]) * static_cast<size_t>(dims[2]);
	unsigned char* ptr = (unsigned char*)data->GetScalarPointer();
	for (size_t i = 0; i < volumeSize; i++)
	{
		if (ptr[i] > chipsIntensities[1])
		{
			ptr[i] = 0;
		}
		if (ptr[i] < chipsIntensities[0])
		{
			ptr[i] = 0;
		}
	}

	//calculate initial porosity
	int min = chipsIntensities[0];
	int max = chipsIntensities[1];

	double totalpercentage = 0.0;

	for (int i = min; i <= max; i++)
	{
		double currentpercentage = static_cast<double>(hisPtr[i]) / dataSize;
		currentpercentage *= 100;
		totalpercentage += currentpercentage;
	}
	double porosity = 100.0 - totalpercentage;
	QString instantstring = "%" + QString::number(porosity, 'g', 8);
	this->ui->initialPorosity->setText(instantstring);
}

ChipsAnalyzer::ChipsAnalyzer()
{
	this->ui = new Ui_ChipsAnalyzer;
	this->ui->setupUi(this);
	connect(ui->minSpinBox, SIGNAL(valueChanged(int)), this, SLOT(on_min_spinbox_changed(int)));
	connect(ui->maxSpinBox, SIGNAL(valueChanged(int)), this, SLOT(on_max_spinbox_changed(int)));
	connect(ui->AirOpacitySpinBox, SIGNAL(valueChanged(double)), this, SLOT(on_airopacity_spinbox_changed(double)));
	connect(ui->ChipsOpacitySpinBox, SIGNAL(valueChanged(double)), this, SLOT(on_chipsopacity_spinbox_changed(double)));
	connect(ui->renderQualitySpinBox, SIGNAL(valueChanged(int)), this, SLOT(on_renderquality_spinbox_changed(int)));
	connect(ui->FolderSelector, SIGNAL(pressed()), this, SLOT(on_folder_selection()));
	// Set up action signals and slots
	connect(this->ui->actionExit, SIGNAL(triggered()), this, SLOT(slotExit()));

}

void ChipsAnalyzer::CalculateTransferFunction()
{

	int minIntensity = ui->minSpinBox->value();
	int maxIntensity = ui->maxSpinBox->value();

	double airOpacity = ui->AirOpacitySpinBox->value();
	double chipsOpacity = ui->ChipsOpacitySpinBox->value();

	colorFunction->RemoveAllPoints();
	opacityFunction->RemoveAllPoints();

	colorFunction->AddRGBPoint(airIntensities[0], 40.0 / 255.0, 137.0 / 255.0, 189.0 / 255.0);
	colorFunction->AddRGBPoint(airIntensities[1], 40.0 / 255.0, 137.0 / 255.0, 189.0 / 255.0);
	colorFunction->AddRGBPoint(chipsIntensities[0], 0.9, 0.8, 0.5);
	colorFunction->AddRGBPoint(chipsIntensities[1], 0.9, 0.8, 0.5);
	opacityFunction->AddPoint(0, 0.0);

	if (minIntensity <= airIntensities[1])
	{
		opacityFunction->AddPoint(minIntensity - 0.0001, 0.0);
		opacityFunction->AddPoint(minIntensity, airOpacity);

		if (maxIntensity <= airIntensities[1])
		{
			opacityFunction->AddPoint(maxIntensity, airOpacity);

		}
		else
		{
			opacityFunction->AddPoint(airIntensities[1], airOpacity);
			opacityFunction->AddPoint(airIntensities[1] + 0.0001, chipsOpacity);
			opacityFunction->AddPoint(maxIntensity, chipsOpacity);
		}

	}
	else
	{
		opacityFunction->AddPoint(minIntensity - 0.0001, 0.0);
		opacityFunction->AddPoint(minIntensity, chipsOpacity);
		opacityFunction->AddPoint(maxIntensity, chipsOpacity);

	}
	if (maxIntensity < 255)
	{
		opacityFunction->AddPoint(maxIntensity + 0.0001, 0.0);
		opacityFunction->AddPoint(255, 0.0);
	}

}

void ChipsAnalyzer::slotExit()
{
	stats->Delete();
	reader->Delete();
	renderer->Delete();
	volume->Delete();
	mapper->Delete();
	colorFunction->Delete();
	opacityFunction->Delete();
	ChipsData->Delete();
	SegmentationMask->Delete();
	qApp->exit();
}

void ChipsAnalyzer::CalculateCurrentPorosity()
{
	vtkIdTypeArray* his = stats->GetHistogram();
	vtkIdType size = his->GetSize();
	vtkIdType* hisPtr = his->GetPointer(0);
	int* dims = ChipsData->GetDimensions();
	double dataSize = dims[0] * dims[1] * dims[2];

	int min = ui->minSpinBox->value();
	int max = ui->maxSpinBox->value();

	double totalPercentage = 0.0;

	for (int i = min; i <= max; i++)
	{
		double currentPercentage = static_cast<double>(hisPtr[i]) / dataSize;
		currentPercentage *= 100;
		totalPercentage += currentPercentage;
	}
	double porosity = 100.0 - totalPercentage;

	QString instantString = "%" + QString::number(porosity, 'g', 8);

	this->ui->instantPorosity->setText(instantString);
}

void ChipsAnalyzer::on_min_spinbox_changed(int value)
{

	int maxValue = ui->maxSpinBox->value();
	if (value > maxValue)
	{
		ui->minSpinBox->setValue(maxValue);
	}
	CalculateTransferFunction();
	CalculateCurrentPorosity();
	this->ui->qvtkWidget->GetRenderWindow()->Render();

}

void ChipsAnalyzer::on_max_spinbox_changed(int value)
{
	int minValue = ui->minSpinBox->value();
	if (value < minValue)
	{
		ui->maxSpinBox->setValue(minValue);
	}
	CalculateTransferFunction();
	CalculateCurrentPorosity();

	this->ui->qvtkWidget->GetRenderWindow()->Render();
}

void ChipsAnalyzer::on_airopacity_spinbox_changed(double value)
{
	CalculateTransferFunction();
	this->ui->qvtkWidget->GetRenderWindow()->Render();
}

void ChipsAnalyzer::on_chipsopacity_spinbox_changed(double value)
{
	CalculateTransferFunction();
	this->ui->qvtkWidget->GetRenderWindow()->Render();
}

void ChipsAnalyzer::on_renderquality_spinbox_changed(int value)
{

	switch (value)
	{
	case 1:
		sampleDistance = 0.1;
		break;
	case 2:
		sampleDistance = 0.5;
		break;
	case 3:
		sampleDistance = 1.0;
		break;
	case 4:
		sampleDistance = 3.0;
		break;
	case 5:
		sampleDistance = 5.0;
		break;
	default:
		sampleDistance = 5.0;
		break;
	}
	mapper->SetSampleDistance(sampleDistance);
	this->ui->qvtkWidget->GetRenderWindow()->Render();

}

int GetNumberOfDigits(int num)
{
	std::string numStr = std::to_string(num);
	return static_cast<int>(numStr.size());
}


void ChipsAnalyzer::WriteChipsAsBinaryImage()
{
	int* dims = ChipsData->GetDimensions();
	unsigned char* ChipsMask = new unsigned char[dims[0] * dims[1] * dims[2]];
	unsigned char* cm = ChipsMask;
	unsigned char* cptr = (unsigned char*)ChipsData->GetScalarPointer();
	size_t totalSize = static_cast<size_t>(dims[0]) * static_cast<size_t>(dims[1]) * static_cast<size_t>(dims[2]);
	memset(ChipsMask, 0, totalSize);
	for (size_t i = 0; i < totalSize; i++)
	{
		if (*cptr >= chipsIntensities[0] && *cptr <= chipsIntensities[1])
		{
			*ChipsMask = 255;
		}
		ChipsMask++;
		cptr++;
	}

	vtkImageData* bitmap = vtkImageData::New();
	bitmap->SetScalarType(VTK_UNSIGNED_CHAR, bitmap->GetInformation());
	bitmap->SetDimensions(dims[0], dims[1], 1);
	bitmap->AllocateScalars(bitmap->GetInformation());

	unsigned char* bptr = (unsigned char*)bitmap->GetScalarPointer();
	ChipsMask = cm;
	int sliceSize = dims[0] * dims[1];

	vtkBMPWriter* writer = vtkBMPWriter::New();

	QString appDir = QCoreApplication::applicationDirPath();
	std::string folderName = appDir.toUtf8().constData();
	folderName += "/" + currentFolderName;
	QDir dir(folderName.c_str());
	if (!dir.exists())
		dir.mkpath(".");
	folderName += "/Bitmap/";
	QDir dir2(folderName.c_str());
	if (!dir2.exists())
		dir2.mkpath(".");

	/*foreach(QString dirFile, dir.entryList())
	{
		dir.remove(dirFile);
	}*/

	for (int z = 0; z < dims[2]; z++)
	{
		memcpy(bptr, ChipsMask, sliceSize);
		writer->SetInputData(bitmap);
		std::string filename = folderName;
		int zeroCount = 5 - GetNumberOfDigits(z);
		for (int k = 0; k < zeroCount; k++)
		{
			filename += "0";
		}
		filename += std::to_string(z);
		filename += ".bmp";
		writer->SetFileName(filename.c_str());
		writer->Write();
		ChipsMask += sliceSize;
	}
	delete[] cm;
	writer->Delete();
	bitmap->Delete();

}

void ChipsAnalyzer::on_folder_selection()
{

	QFileDialog dialog;
	dialog.setFileMode(QFileDialog::DirectoryOnly);
	dialog.setOption(QFileDialog::ShowDirsOnly, false);
	QString directory = dialog.getExistingDirectory();
	std::string folderName = directory.toUtf8().constData();
	if (folderName == "")
	{
		return;
	}
	size_t fIndex = folderName.find_last_of("/\\");
	currentFolderName = folderName.substr(fIndex + 1);
	std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

	std::cout << "Reading chips data started" << std::endl;
	ImageLoader* loader = new ImageLoader();
	loader->SetFolderName(folderName);
	loader->Read();
	std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();

	long long duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

	std::cout << "Reading chips data finished in " << duration << " ms" << std::endl;

	start = std::chrono::high_resolution_clock::now();
	std::cout << "Volume resampling started " << std::endl;
	double reductionFactor = ui->reductionFactor->value();
	vtkImageResample *resample = vtkImageResample::New();
	if (reductionFactor < 0.98)
	{
		resample->SetInputData(loader->GetOutput());
		resample->SetAxisMagnificationFactor(0, reductionFactor);
		resample->SetAxisMagnificationFactor(1, reductionFactor);
		resample->SetAxisMagnificationFactor(2, reductionFactor);
		resample->SetNumberOfThreads(8);
		resample->Update();
	}
	end = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << "Volume resampling finished in " << duration << " ms" << std::endl;
	if (ChipsData != NULL)
	{
		ChipsData->Delete();
		ChipsData = NULL;
	}
	ChipsData = vtkImageData::New();
	if (reductionFactor < 0.99)
	{
		ChipsData->ShallowCopy(resample->GetOutput());
	}
	else
		ChipsData->ShallowCopy(loader->GetOutput());

	delete loader;
	resample->Delete();
	PreProcessData(ChipsData);

	if (ui->writeChipsMask->checkState() == Qt::Checked)
	{
		std::cout << "Writing Chips Mask Binary Images started" << std::endl;
		WriteChipsAsBinaryImage();
		std::cout << "Writing Chips Mask Binary Images finished" << std::endl;
	}

	start = std::chrono::high_resolution_clock::now();
	std::cout << "Surface extraction started" << std::endl;

	int neighbourHood = ui->neighbourHood->currentIndex();

	Segmentation = new ChipsSegmentation();
	Segmentation->SetAirIntensities(airIntensities[0], airIntensities[1]);
	Segmentation->SetChipsIntensities(chipsIntensities[0], chipsIntensities[1]);
	Segmentation->SetInputData(ChipsData);
	Segmentation->SetNeighbourHood(neighbourHood);
	Segmentation->ExtractSurface();

	end = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << "Surface extraction finished in " << duration << " ms" << std::endl;

	start = std::chrono::high_resolution_clock::now();
	std::cout << "Connectivity Index calculation started" << std::endl;

	Segmentation->LabelAndAnalyze();

	end = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << "Connectivity Index calculation finished in " << duration << " ms" << std::endl;

	/*if (SegmentationMask == NULL)
	{
		SegmentationMask = vtkImageData::New();
	}

	this->SegmentationMask->DeepCopy(Segmentation->GetMask());*/

	QString instantString = "%" + QString::number(Segmentation->GetInitialPorosity(), 'g', 8);
	this->ui->segmentedPorosity->setText(instantString);

	instantString = QString::number(Segmentation->GetConnectivityIndex(), 'g', 8);
	this->ui->connectivityIndex->setText(instantString);

	if (ui->writeChipsBoundaryMask->checkState() == Qt::Checked)
	{
		start = std::chrono::high_resolution_clock::now();
		std::cout << "WriteChipsBoundary started" << std::endl;
		Segmentation->WriteChipsBoundary(currentFolderName);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "WriteChipsBoundary finished in " << duration << " ms" << std::endl;
	}

	if (ui->WriteConnectedComponents->checkState() == Qt::Checked)
	{
		start = std::chrono::high_resolution_clock::now();
		std::cout << "WriteConnectedComponents started" << std::endl;
		Segmentation->WriteConnectedComponents(currentFolderName);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "WriteConnectedComponents finished in " << duration << " ms" << std::endl;
	}

	//delete Segmentation;

	if (renderer != NULL)
	{
		renderer->Delete();
		renderer = NULL;
	}

	if (volume != NULL)
	{
		volume->Delete();
		volume = NULL;
	}

	if (mapper != NULL)
	{
		mapper->Delete();
		mapper = NULL;
	}

	if (colorFunction != NULL)
	{
		colorFunction->Delete();
		colorFunction = NULL;
	}

	if (opacityFunction != NULL)
	{
		opacityFunction->Delete();
		opacityFunction = NULL;
	}

	if (property != NULL)
	{
		property->Delete();
		property = NULL;
	}

	renderer = vtkRenderer::New();

	volume = vtkVolume::New();
	mapper = vtkOpenGLGPUVolumeRayCastMapper::New();
	mapper->SetMaxMemoryFraction(0.95);
	mapper->SetInputData(ChipsData);

	colorFunction = vtkColorTransferFunction::New();
	opacityFunction = vtkPiecewiseFunction::New();

	property = vtkVolumeProperty::New();
	property->SetIndependentComponents(true);
	property->SetColor(colorFunction);
	property->SetScalarOpacity(opacityFunction);
	property->SetInterpolationTypeToLinear();
	property->ShadeOn();
	property->SetDiffuse(1.0);
	property->SetAmbient(0.15);
	//property->SetSpecular(0.005);
	//property->SetSpecularPower(10);
	volume->SetProperty(property);
	volume->SetMapper(mapper);

	CalculateTransferFunction();

	/*colorFunction->AddRGBPoint(0, 0.0, 0.0, 0.0);
	colorFunction->AddRGBPoint(2, 0.0, 0.0, 0.0);
	colorFunction->AddRGBPoint(3, 0.0, 0.2, 0.3);
	colorFunction->AddRGBPoint(30, 0.0, 0.2, 0.3);
	colorFunction->AddRGBPoint(80, 0.9, 0.8, 0.5);
	colorFunction->AddRGBPoint(255, 0.9, 0.8, 0.5);*/
	/*colorFunction->AddRGBPoint(255, 0.9, 0.02, 0.05);*/

	/*opacityFunction->AddPoint(0, 0.00);
	opacityFunction->AddPoint(30, 0.0);
	opacityFunction->AddPoint(31, 0.1);
	opacityFunction->AddPoint(79, 0.1);
	opacityFunction->AddPoint(80, 0.10);
	opacityFunction->AddPoint(255, 0.1);*/

	mapper->SetBlendModeToComposite();
	mapper->AutoAdjustSampleDistancesOff();
	mapper->SetSampleDistance(sampleDistance);

	//mapper->AutoAdjustSampleDistancesOn();
	/*mapper->SetSampleDistance(0.1);
	mapper->SetImageSampleDistance(1);*/
	renderer->AddVolume(volume);

	renderer->ResetCamera();

	this->ui->qvtkWidget->GetRenderWindow()->AddRenderer(renderer);
	ui->minSpinBox->setValue(chipsIntensities[0]);
	ui->maxSpinBox->setValue(chipsIntensities[1]);
	this->ui->qvtkWidget->GetRenderWindow()->Render();
}