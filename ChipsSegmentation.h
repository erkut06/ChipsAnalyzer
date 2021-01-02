#pragma once
#include <cstddef>
#include <string>
class vtkImageData;

class ChipsSegmentation
{
public:

	ChipsSegmentation();
	~ChipsSegmentation();
	void SetInputData(vtkImageData* input);
	void SetAirIntensities(unsigned char airMin, unsigned char airMax);
	void SetChipsIntensities(unsigned char cMn, unsigned char cMx);
	void SetNeighbourHood(int n) { neighbourHood = n; }
	void ExtractSurface();
	double GetInitialPorosity() { return initialPorosity; }
	double GetConnectivityIndex() { return connectivityIndex; }
	double GetOutRatio() { return outRatio; }
	void WriteChipsBoundary(std::string fN);
	void WriteConnectedComponents(std::string fN);
	void LabelAndAnalyze();
	vtkImageData* GetMask();

private:

	void ExtractSliceSurface(int sliceNum, double maxRadius, unsigned short forbiddenFruit);
	void CalculateChipsMask();
	void FillChipsMask();
	void CalculateInitialPorosity();
	vtkImageData* inputData = nullptr;
	vtkImageData* chipsMask = nullptr;
	vtkImageData* connectedMask = nullptr;
	unsigned char airMin = 0;
	unsigned char airMax = 1;
	unsigned char chipsMin = 0;
	unsigned char chipsMax = 1;
	double initialPorosity = 100;
	double connectivityIndex = 0;
	double outRatio = 0;

	//0 6-n
	//1 26-n
	int neighbourHood = 0;
};

