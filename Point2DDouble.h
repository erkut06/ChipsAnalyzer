#pragma once

class Point2DDouble
{
public:
	Point2DDouble();
	~Point2DDouble();
	double x = 0;
	double y = 0;
	int GetImageIndex(int* dims);
private:

};

