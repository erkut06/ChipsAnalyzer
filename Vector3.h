#pragma once
#include <math.h>

class Vector3
{
public:
    Vector3();
    Vector3(double x, double y, double z);
    Vector3 Cross(Vector3& );
    double Dot(Vector3 );
    void Normalize();
    double Length();
    Vector3 operator +(Vector3& );
	double X;
	double Y;
	double Z;
};

