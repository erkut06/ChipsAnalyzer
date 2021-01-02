#pragma once

#include <vtkSmartPointer.h>
#include <QMainWindow>

#include <memory>

class Ui_ChipsAnalyzer;
class vtkRenderer;
class vtkColorTransferFunction;
class vtkPiecewiseFunction;
class vtkImageData;
class vtkMetaImageReader;
class vtkVolume;
class vtkOpenGLGPUVolumeRayCastMapper;
class vtkVolumeProperty;
class vtkImageHistogramStatistics;
class ChipsSegmentation;

class ChipsAnalyzer : public QMainWindow
{
	Q_OBJECT
public:

	ChipsAnalyzer();
	~ChipsAnalyzer() {};
	void PreProcessData(vtkImageData*);
public slots:

	virtual void slotExit();

private slots:
	void on_min_spinbox_changed(int);
	void on_max_spinbox_changed(int);
	void on_airopacity_spinbox_changed(double);
	void on_chipsopacity_spinbox_changed(double);
	void on_renderquality_spinbox_changed(int);
	void on_folder_selection();
private:

	vtkMetaImageReader* reader = nullptr;
	vtkVolume* volume = nullptr;
	vtkRenderer* renderer = nullptr;
	vtkOpenGLGPUVolumeRayCastMapper* mapper = nullptr;
	vtkColorTransferFunction* colorFunction = nullptr;
	vtkPiecewiseFunction* opacityFunction = nullptr;
	vtkVolumeProperty* property = nullptr;
	vtkImageData* chipsData = nullptr;
	vtkImageHistogramStatistics* stats = nullptr;
	vtkImageData* segmentationMask = nullptr;
	std::shared_ptr<ChipsSegmentation> segmentation = nullptr;
	double airIntensities[2];
	double chipsIntensities[2];
	std::string currentFolderName = "";
	//1 0.1
	//2 0.5
	//3 1.0
	//4 3.0
	//5 5.0
	double sampleDistance = 5;
	Ui_ChipsAnalyzer *ui;
	void CalculateTransferFunction();
	void CalculateCurrentPorosity();
	void WriteChipsAsBinaryImage();
};
