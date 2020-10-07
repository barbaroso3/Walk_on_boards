#include <iostream>
#include <array>
#include <vector>
#include <math.h>
#include <fstream>

double RAND_01()
{
	double U;
	do U = rand() / (double)RAND_MAX;
	while ((U == 0) || (U == 1));

	return U;
}

class Point
{
private:
	double m_x;
	double m_y;
	double m_z;

public:
	double dist_to(Point b) const {
		
		return sqrt(pow(m_x - b.get_x(), 2) + pow(m_y - b.get_y(), 2) + pow(m_z - b.get_z(), 2));
	}

	Point(): m_x(0), m_y(0), m_z(0) {}
	
	Point(double x, double y, double z): m_x(x), m_y(y), m_z(z) {}
	
	double get_x () const { return m_x; }
	
	double get_y () const { return m_y; }
	
	double get_z () const { return m_z; }
	 
	void set_x (double x) { m_x = x; }
	
	void set_y (double y) { m_y = y; }
	
	void set_z (double z) { m_z = z; }
	
	friend std::ostream& operator<< (std::ostream &out, const Point &point)
	{
		out << "(" << point.m_x << ", " << point.m_y << ", " << point.m_z << ")";
		
		return out;
	}
};

struct DomainPrism
{
	// h - ширина призмы
	double h = 0;
	// l - длина призмы
	double l = 0;
};

struct WoB_answer
{
	double expect1 = 0;
	double expect2 = 0;
	double variance = 0;
	double stat_error = 0;
};


int main() {

	
	
	return 0;
}
