#pragma once
class Point
{
private:
	double x = 0;
	double y = 0;
public:
	Point();
	Point(double x, double y);
	double getXValue();
	double getYValue();
	~Point();
};

