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

Point operator+ (const Point &a, const Point &b)
{
	return Point(a.get_x() + b.get_x(), a.get_y() + b.get_y(), a.get_z() + b.get_z());
}

Point operator- (const Point &a, const Point &b)
{
	return Point(a.get_x() - b.get_x(), a.get_y() - b.get_y(), a.get_z() - b.get_z());
}

class Vector3
{
private:
	double m_x;
	double m_y;
	double m_z;
	
public:
	Vector3(): m_x(0), m_y(0), m_z(0) {}
	
	Vector3(double x, double y, double z): m_x(x), m_y(y), m_z(z) {}
	
	Vector3(Point x): m_x(x.get_x()), m_y(x.get_y()), m_z(x.get_z()) {}
	
	double get_x () const { return m_x; }
	
	double get_y () const { return m_y; }
	
	double get_z () const { return m_z; }
	 
	void set_x (double x) { m_x = x; }
	
	void set_y (double y) { m_y = y; }
	
	void set_z (double z) { m_z = z; }
	
	friend std::ostream& operator<< (std::ostream &out, const Vector3 &vector)
	{
		out << "(" << vector.m_x << ", " << vector.m_y << ", " << vector.m_z << ")";
		
		return out;
	}
};

double Scholar(const Vector3 &a, const Vector3 &b)
{
	return a.get_x() * b.get_x() + a.get_y() * b.get_y() + a.get_z()*b.get_z();
}

Vector3 operator* (const Vector3 &a, double x)
{
	return Point(a.get_x()*x, a.get_y()*x, a.get_z()*x);
}

Point operator+(const Vector3 &vec, const Point &p)
{
	return Point(vec.get_x() + p.get_x() ,vec.get_y() + p.get_y(), vec.get_z() + p.get_z());
}

class DomainPrism
{
public:
	double const h;
	double const l;
	
	struct Prism_plane
	{
		Point p0;
		Vector3 normal;
		std::string name;
	};
	
	std::array<Prism_plane, 6> plane;
		
	DomainPrism(double h1, double l1): h(h1), l(l1)
	{
		plane[0].name = "z_0";
		plane[0].p0 = Point(h/2, h/2 , 0);
		plane[0].normal = Vector3(0, 0, 1);
		
		plane[1].name = "z_l";
		plane[1].p0 = Point(h/2, h/2 , l);
		plane[1].normal = Vector3(0, 0, 1);
		
		plane[2].name = "x_0";
		plane[2].p0 = Point(0, h/2 , l/2);
		plane[2].normal = Vector3(1, 0, 0);
		
		plane[3].name = "x_h";
		plane[3].p0 = Point(h, h/2 , l/2);
		plane[3].normal = Vector3(1, 0, 0);
		
		plane[4].name = "y_0";
		plane[4].p0 = Point(h/2, 0 , l/2);
		plane[4].normal = Vector3(0, 1, 0);
		
		plane[5].name = "y_h";
		plane[5].p0 = Point(h/2, h , l/2);
		plane[5].normal = Vector3(0, 1, 0);
	}
};

struct WOB_answer
{
	double expect1 = 0;
	double expect2 = 0;
	double variance = 0;
	double stat_error = 0;
};

class WOB_parameters
{
private:
	int N;
	Point x0;
	int num_walks;
	
public:
	WOB_parameters(): N(0) {}
	
	void set_N(int N1) { N = N1; }
	int get_N() const {return N; }
	
	void set_x0 (Point x) {x0 = x; }
	Point get_x0 () const { return x0; }
	
	void set_num_walks (int num_walks1) {num_walks = num_walks1;}
	int get_num_walks () const {return num_walks; }
	
};

Vector3 Simulate_vector3_isotropic ()
{
	double cosT = 1 - 2 * RAND_01();
	double fi = 2 * M_PI * RAND_01();
	
	return Vector3(sin(fi)*sin(acos(cosT)), cos(fi)*sin(acos(cosT)), cosT);
}

Point move_point_to_border (DomainPrism prism, Point near_border_point)
{
	Point border_point = near_border_point;
	
	double eps = pow(10, -15);
	
	if (border_point.get_x() >= -eps && border_point.get_x() <= eps) border_point.set_x(0);
	else if (border_point.get_x() >= prism.h - eps && border_point.get_x() <= prism.h + eps) border_point.set_x(prism.h);
	
	if (border_point.get_y() >= -eps && border_point.get_y() <= eps) border_point.set_y(0);
	else if (border_point.get_y() >= prism.h - eps && border_point.get_y() <= prism.h + eps) border_point.set_y(prism.h);
	
	if (border_point.get_z() >= -eps && border_point.get_z() <= eps) border_point.set_z(0);
	else if (border_point.get_z() >= prism.l - eps && border_point.get_z() <= prism.l + eps) border_point.set_z(prism.l);
	
	
	
	return border_point;
}

Point Simulate_first_border_point(Point const &start_point, DomainPrism const &prism)
{
	Vector3 ray_direc(Simulate_vector3_isotropic());
	
	double t = sqrt(pow(prism.h, 2) + pow(prism.h, 2) + pow(prism.l, 2));
	for(auto const & plane:prism.plane)
	{
		double t_version = Scholar (Vector3(plane.p0 - start_point), plane.normal) / Scholar(ray_direc, plane.normal);
		if(t_version > 0 && t_version < t) t = t_version;
	}
	Point t_dist(ray_direc.get_x()*t, ray_direc.get_y()*t, ray_direc.get_z()*t);
	
	Point border_point = move_point_to_border(prism,  start_point + t_dist);
	
	return border_point;
	
}

bool belong_to_prsim(Point const &point, DomainPrism const &prism)
{
	double eps = pow(10, -15);
	
	if (point.get_x() >= -eps && point.get_x() <= prism.h + eps)
		if (point.get_y() >= -eps && point.get_y() <= prism.h + eps)
			if (point.get_z() >= -eps && point.get_z() <= prism.l + eps)
				return true;
	
	return false;
}

Point Simulate_next_border_point(Point const &start_point, DomainPrism const &prism)
{
	Vector3 ray_direc(Simulate_vector3_isotropic());
	
	Point next_border_point;
	
	for(auto const & plane:prism.plane)
	{
		double t = Scholar (Vector3(plane.p0 - start_point), plane.normal) / Scholar(ray_direc, plane.normal);
		Point t_dist(ray_direc.get_x()*t, ray_direc.get_y()*t, ray_direc.get_z()*t);
		
		Point possible_point = start_point + t_dist;
		if(belong_to_prsim(possible_point, prism) && t != 0)
			next_border_point = possible_point;
	}
	
	next_border_point = move_point_to_border(prism, next_border_point);
	
	return next_border_point;
}

double Test_function(Point p)
{
	return pow(p.get_x(), 2) + pow(p.get_y(), 2) - 2 * pow(p.get_z(), 2) + 5;
}

WOB_answer WOB_algorithm(WOB_parameters parameters, DomainPrism prism)
{
	WOB_answer answer;
	
	Point source_point(parameters.get_x0());
	
	for(int it_N = 0; it_N < parameters.get_N(); ++it_N)
	{
		if (it_N % static_cast<int>(parameters.get_N() / 100.0) == 0)
			std::cout << (double)it_N/parameters.get_N() * 100 << "%" << std::endl;
		
		double sum = 0;
		Point current_point(Simulate_first_border_point(source_point, prism));
		
		if(current_point.get_z() == prism.l)
		{
			sum += 2;
//			sum += 2 * Test_function(current_point);
		}

		for(int it_walks = 1; it_walks <= parameters.get_num_walks(); ++it_walks)
		{
			current_point = Simulate_next_border_point(current_point, prism);

			
//			double weight = 2 * pow(-1, it_walks) * Test_function(current_point);
//			if(it_walks == parameters.get_num_walks())
//				weight *= 0.5;
//
//			sum += weight;
			
			if(current_point.get_z() == prism.l)
			{
				double weight = 2 * pow(-1, it_walks);
				if(it_walks == parameters.get_num_walks())
					weight *= 0.5;

				sum += weight;
			}
		}
		
		answer.expect1 += sum;
		answer.expect2 += pow(sum, 2);
	}
	
	answer.expect1 /= parameters.get_N();
	answer.expect2 /= parameters.get_N();
	answer.variance = answer.expect2 - pow(answer.expect1, 2);
	answer.stat_error = 3 * sqrt( answer.variance / parameters.get_N());
	
	return answer;
}

int main() {
	
	srand(static_cast<unsigned int>(time(nullptr)));
	
	double h = 1;
	double L = 3;
	DomainPrism prism (h, L);
	
	WOB_parameters parameters;
	parameters.set_N(pow(10, 8));
	parameters.set_x0(Point(h/2, h/2, h/2));
	parameters.set_num_walks(10);
	
	std::ofstream answ_File;
	answ_File.open("WOB_l3_N8.txt", std::ios_base::app);
	answ_File << "num_walks answ stat_eror time" << std::endl;
	
	for(int i = 4; i <= 8; ++i)
	{
		parameters.set_num_walks(i * 4);
		
		clock_t start = clock();
		WOB_answer answer = WOB_algorithm(parameters, prism);
		double time = (double)(clock() - start)/(CLOCKS_PER_SEC);
		
		std::cout <<
		parameters.get_num_walks() << " " <<
		answer.expect1 << " " <<
		answer.stat_error << " " <<
		time << " " <<
		std::endl;
		
		answ_File <<
		parameters.get_num_walks() << " " <<
		answer.expect1 << " " <<
		answer.stat_error << " " <<
		time << " " <<
		std::endl;
	}
	
	return 0;
}
