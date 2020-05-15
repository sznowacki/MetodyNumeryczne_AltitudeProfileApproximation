#include "Point.h"



Point::Point()
{
}

Point::Point(double x, double y)
{
	this->x = x;
	this->y = y;
}

double Point::getXValue()
{
	return this->x;
}

double Point::getYValue()
{
	return this->y;
}


Point::~Point()
{
}
