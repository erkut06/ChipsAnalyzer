#include "Point2DDouble.h"

Point2DDouble::Point2DDouble()
{

}

Point2DDouble::~Point2DDouble()
{

}

int Point2DDouble::GetImageIndex(int* dims)
{
	int xInt = static_cast<int>(x);
	int yInt = static_cast<int>(y);
	if (xInt < 0 || xInt >= dims[0] || yInt < 0 || yInt >= dims[1])
	{
		return -1;
	}
	return dims[0] * yInt + xInt;
}