#include "BoundaryPoint2D.h"
#include <cmath>
BoundaryPoint2D::BoundaryPoint2D()
{
}

BoundaryPoint2D::~BoundaryPoint2D()
{
}

void BoundaryPoint2D::SetRadius(double r)
{
	if (radius != r)
	{
		radius = r;
		radiusChanged = true;
	}
}

void BoundaryPoint2D::SetAngle(double a)
{
	if (teta != a)
	{
		teta = a;
		angleChanged = true;
	}
}

void BoundaryPoint2D::IncrementRadius()
{
	if (!fixed)
	{
		radius += 1.0;
		radiusChanged = true;
		UpdatePoint(center.x, center.y);
	}
}

void BoundaryPoint2D::DecrementRadius()
{
	if (!fixed)
	{
		radius -= 1.0;
		radiusChanged = true;
		UpdatePoint(center.x, center.y);
	}
}

void BoundaryPoint2D::UpdatePoint(double& cx, double& cy)
{
	if (angleChanged)
	{
		cosTeta = std::cos(teta);
		sinTeta = std::sin(teta);
		x = cx + radius * cosTeta;
		y = cy + radius * sinTeta;
		angleChanged = false;
		radiusChanged = false;
	}
	else if (radiusChanged)
	{
		x = cx + radius * cosTeta;
		y = cy + radius * sinTeta;
		radiusChanged = false;
		
	}
	else if(center.x != cx || center.y != cy)
	{
		x += cx - center.x;
		y += cy - center.y;
	}
	center.x = cx;
	center.y = cy;
}

void BoundaryPoint2D::UpdateCenter(double cx, double cy)
{
	double currentX = x;
	double currentY = y;
}