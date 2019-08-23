#pragma once
#include "BoundaryPoint2D.h"
class BoundaryModel2D
{
public:

	BoundaryModel2D();
	~BoundaryModel2D();
	void SetCenter(double cx, double cy);
	void AllocatePoints(int depth, double);
	void FitBoundary(unsigned short* dPtr, int* dims, unsigned short& forbiddenFruit);
	void MarchAll(unsigned short* dPtr, int* dims, unsigned short&);

private:

	void CheckBoundary(unsigned short* dPtr, int* dims, unsigned short& forbiddenFruit, unsigned char* pointMap);
	void ClearPointsIfAllocated();
	BoundaryPoint2D** points = NULL;
	int pointSize = 0;
	Point2DDouble center;

};

