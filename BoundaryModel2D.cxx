#include "BoundaryModel2D.h"
#include "Vector3.h"
#include "vtkImageData.h"
#include <algorithm>
#include <math.h>

#define PI 3.14159265

BoundaryModel2D::BoundaryModel2D()
{

}

BoundaryModel2D::~BoundaryModel2D()
{
	points.clear();
}

void BoundaryModel2D::SetCenter(double cx, double cy)
{
	center.x = cx;
	center.y = cy;
}

void BoundaryModel2D::ClearPointsIfAllocated()
{
	points.clear();
}

bool IsMarchingDone(unsigned char* PointMap, int size)
{
	for (auto i = 0; i < size ; i++)
	{
		if (PointMap[i] == 0)
		{
			return false;
		}
	}
	return true;
}

void BoundaryModel2D::MarchAll(unsigned short* dPtr, int* dims, unsigned short& forbiddentFruit)
{
	auto maxIteration = sqrt(dims[0] * dims[0] + dims[1] * dims[1]) * 10;
	auto pointMap = new unsigned char[pointSize];
	memset(pointMap, 0, pointSize);
	unsigned char* tempPointMap = new unsigned char[pointSize];
	memset(tempPointMap, 0, pointSize);

	int iterationCount = 0;

	while (!IsMarchingDone(pointMap, pointSize))
	{
		for (int i = 0; i < pointSize; i++)
		{
			if (pointMap[i] == 0)
			{
				auto& current = points[i];
				current->DecrementRadius();
			}
		}
		memcpy(tempPointMap, pointMap, pointSize);
		CheckBoundary(dPtr, dims, forbiddentFruit, pointMap);

		for (int i = 0; i < pointSize; i++)
		{
			if (tempPointMap[i] != pointMap[i] )
			{
				auto& current = points[i];
				current->IncrementRadius();
				current->fixed = true;
			}
		}
		iterationCount++;
		if (iterationCount > maxIteration)
		{
			break;
		}
	}

	int incrementCount = 5;
	for (auto& current: points)
	{
		current->fixed = false;
		for (int inc = 0; inc < incrementCount; inc++)
		{
			current->IncrementRadius();
		}
	}
}

void BoundaryModel2D::FitBoundary(unsigned short* dPtr, int* dims, unsigned short& forbiddenFruit)
{
	MarchAll(dPtr, dims, forbiddenFruit);

	Point2DDouble lerp;
	int currentIndex = 0;
	for (int i = 0; i < pointSize; i++)
	{
		auto& currentPoint = points[i];
		auto& nextPoint = points[(i + 1) % pointSize];

		Vector3 way(nextPoint->x - currentPoint->x, nextPoint->y - currentPoint->y, 0);
		double length = static_cast<double>(way.Length());

		unsigned short changeValue = forbiddenFruit / 2;

		way.Normalize();
		for (double l = 0; l <= length; l++)
		{
			lerp.x = way.X * l + currentPoint->x;
			lerp.y = way.Y * l + currentPoint->y;
			currentIndex = lerp.GetImageIndex(dims);
			if (currentIndex != -1)
			{
				dPtr[currentIndex] = changeValue;
			}
		}
	}

	unsigned short borderValue = forbiddenFruit / 2;
	unsigned short outsideValue = forbiddenFruit / 4;
	unsigned short* maskSlice = dPtr;

	for (int y = 0; y < dims[1]; y++)
	{
		//left check
		for (int x = 0; x < dims[0]; x++)
		{
			int sliceIndex = y * dims[0] + x;

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
			int sliceIndex = y * dims[0] + x;
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
			int sliceIndex = y * dims[0] + x;
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
			int sliceIndex = y * dims[0] + x;
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

void BoundaryModel2D::CheckBoundary(unsigned short* dPtr, int* dims, unsigned short& forbiddenFruit, unsigned char* pointMap)
{

	unsigned short currentIntensity = 0;
	int currentIndex = 0;
	Point2DDouble lerp;

	for (int i = 0; i < pointSize; i++)
	{
		auto& currentPoint = points[i];
		auto& nextPoint = points[(i + 1) % pointSize];

		std::shared_ptr<BoundaryPoint2D> previousPoint = NULL;
		if (i != 0)
		{
			previousPoint = points[i - 1];
		}
		else
		{
			previousPoint = points[pointSize - 1];
		}
		Vector3 way(nextPoint->x - currentPoint->x, nextPoint->y - currentPoint->y, 0);
		int length = static_cast<int>(way.Length());
		if (length <= 2)
		{
			currentIndex = currentPoint->GetImageIndex(dims);
			if (currentIndex != -1)
			{
				currentIntensity = dPtr[currentIndex];
				if (currentIntensity == forbiddenFruit)
				{
					pointMap[i] = 1;
				}
			}
		}
		else
		{
			way.Normalize();
			for (double l = 1; l < length; l++)
			{
				lerp.x = way.X * l + currentPoint->x;
				lerp.y = way.Y * l + currentPoint->y;
				currentIndex = lerp.GetImageIndex(dims);
				if (currentIndex != -1)
				{
					currentIntensity = dPtr[currentIndex];
					if (currentIntensity == forbiddenFruit)
					{
						pointMap[i] = 1;
					}
				}
			}
		}

		Vector3 line1(currentPoint->x - previousPoint->x, currentPoint->y - previousPoint->y, 1);
		Vector3 line2(currentPoint->x - nextPoint->x, currentPoint->y - nextPoint->y, 1);
		line1.Normalize();
		line2.Normalize();

		double radAngle = acos(line1.Dot(line2));

		if (radAngle < 0)
		{
			radAngle += (2 * PI);
		}

		if (radAngle < (PI / 1.1))
		{
			pointMap[i] = 1;
		}

	}
}

void BoundaryModel2D::AllocatePoints(int depth, double radius)
{
	if (depth < 3)
	{
		depth = 3;
	}
	if (depth > 16)
	{
		depth = 16;
	}
	ClearPointsIfAllocated();

#ifdef __APPLE__
	pointSize = pow(2, depth);
#else
	pointSize = std::pow(2, depth);
#endif

	double pi = 22.0 / 7.0;
	double increment = (2 * pi) / static_cast<double>(pointSize);
	int pointIndex = 0;
	for (double teta = 0.0; teta < 2 * pi; teta += increment)
	{
		if (pointIndex >= pointSize)
		{
			break;
		}
		auto newBoundaryPoint = std::shared_ptr<BoundaryPoint2D>(new BoundaryPoint2D());
		points.push_back(newBoundaryPoint);
		newBoundaryPoint->SetAngle(teta);
		newBoundaryPoint->SetRadius(radius);
		newBoundaryPoint->UpdatePoint(center.x, center.y);
		pointIndex++;
	}
}
