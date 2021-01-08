/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Shreya Bhaumik
 * *************************
*/
struct Vec3 // 3d vector
{
	double x;
	double y;
	double z;

	Vec3() {x = double(0.0); y = double(0.0); z = double(0.0);} // default constructor
	Vec3(double _x, double _y, double _z) {x=_x; y=_y; z=_z;} // parameterized constructor

	double dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }
	double len2() const { return (x*x + y*y + z*z); }
	double len() const { return std::sqrt(len2()); }
	Vec3& Normalize()
	{
	  double d = len2();
	  if (d > 0) {
	    double dinv = (double)1.0 / std::sqrt(d);
	    x *= dinv;
	    y *= dinv;
	    z *= dinv;
	  }
	  return *this;
	}
	Vec3& cross(const Vec3& v1, const Vec3& v2)
	{
	  x = v1.y * v2.z - v1.z * v2.y;
	  y = v1.z * v2.x - v1.x * v2.z;
	  z = v1.x * v2.y - v1.y * v2.x;
	  return *this;
	}
	// -- Unary operators --
	Vec3 operator-() const { return Vec3(-x, -y, -z); }
	// -- Binary operators --
	Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
	Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
	Vec3 operator*(double scalar) const { return Vec3(x * scalar, y * scalar, z * scalar); }
	Vec3& operator+=(const Vec3& v) 
	{
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}
	Vec3& operator-=(const Vec3& v) 
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}
};

struct Color // RGB
{
	double r;
	double g;
	double b;

	Color() {r = 0.0; g = 0.0; b = 0.0;} // default constructor
	Color(double _r, double _g, double _b) {r=_r; g=_g; b=_b;} // parameterized constructor

	// functions
	Color& colClamp()
	{
		if(r > 1.0) r = 1.0;
		else if(r < 0.0) r = 0.0;
		if(g > 1.0) g = 1.0;
		else if(g < 0.0) g = 0.0;
		if(b > 1.0) b = 1.0;
		else if(b < 0.0) b = 0.0;
		return *this;
	}
	Color& operator+=(Color const &c) { r+=c.r; g+=c.g; b+=c.b; return *this; }
	Color operator*(double scalar) const { return Color(r * scalar, g * scalar, b * scalar); }
	Color& operator*=(double scalar) { r*=(double)scalar; g*=(double)scalar; b*=(double)scalar; return *this; }
	Color& operator/=(double scalar) { r/=(double)scalar; g/=(double)scalar; b/=(double)scalar; return *this; }
	// -- Binary operators --
	Color operator+(Color const &c) const { return Color(r+c.r, g+c.g, b+c.b).colClamp(); }
};