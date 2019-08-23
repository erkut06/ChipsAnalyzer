#include "Vector3.h"

Vector3::Vector3()
{
	X = 0.0;
	Y = 0.0;
	Z = 0.0;
}

Vector3::Vector3(double x, double y, double z)
{
	X = x;
	Y = y;
	Z = z;
}

Vector3 Vector3::Cross(Vector3& crossWith)
{
  Vector3 result;
  result.X = Y * crossWith.Z - Z * crossWith.Y;
  result.Y = Z * crossWith.X - X * crossWith.Z;
  result.Z = X * crossWith.Y - Y * crossWith.X;
  result.Normalize();
  return result;
}

double Vector3::Dot(Vector3 dotV)
{
  return (X * dotV.X) + (Y * dotV.Y) + (Z * dotV.Z);
}

void Vector3::Normalize()
{
  double length = sqrt(X * X + Y * Y + Z * Z);

  if (length != 0.0)
  {
    X /= length;
    Y /= length;
    Z /= length;
  }
}
double Vector3::Length()
{
	return sqrt(X * X + Y * Y + Z * Z);
}

Vector3 Vector3::operator +(Vector3& vectorToAdd)
{
  return Vector3(X + vectorToAdd.X, Y + vectorToAdd.Y, Z + vectorToAdd.Z);
}


