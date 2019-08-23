#pragma once
#include "DataStructures.h"
#include "Point2DDouble.h"
class BoundaryPoint2D:public Point2DDouble
{
public:

	BoundaryPoint2D();
	~BoundaryPoint2D();

	void UpdatePoint(double& cx, double& cy);
	void UpdateCenter(double cx, double cy);
	void SetRadius(double r);
	void SetAngle(double a);
	void DecrementRadius();
	void IncrementRadius();
	bool fixed = false;
	double radius = 1.0;
private:
	
	double teta = 0.0;
	double sinTeta = 0.0;
	double cosTeta = 0.0;
	Point2DDouble center;
	bool angleChanged = true;
	bool radiusChanged = true;
};

