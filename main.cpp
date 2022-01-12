#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"




class Vector{
	public: 
	explicit Vector(double x=0 , double y=0, double z=0){
		data[0] = x;
		data[1] = y; 
		data[2] = z;
		
	}
	double norm2() const {
		return data[0]*data[0] + data[1]*data[1] + data[2]*data[2];
	}
	double norm() const {
		return sqrt(norm2());
	}
	void normalize(){
		double n = norm();

		data[0] /= n;
		data[1] /= n;
		data[2] /=n;
	}



	double operator[](int i) const {return data[i];}; //lecture seule
	double& operator[](int i){return data[i];}; //modification
	double data[3];
};

Vector operator+(const Vector& a, const Vector& b){
	return Vector(a[0]+b[0], a[1]+b[1], a[2]+b[2]);
}

Vector operator-(const Vector&a, const Vector& b){
	return Vector(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}

Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class Ray {
public:
	Ray(const Vector& origin, const Vector&  direction) : C(origin),u(direction){

};

Vector C,u;

};
 
class Sphere {
public:
	Sphere( const Vector& O, const double R) : O(O),R(R){

};

	bool intersect(const Ray& r , Vector& P, Vector& N) {
		//solve a*t^2 +b*t +c = 0
		double a =1;
		double b= 2*dot(r.u , r.C - O );
		double c= (r.C - O).norm2()- R*R;

		double delta = b*b - 4*a*c;
		if (delta < 0){
			return false ;
		}
		double t1 = (-b - sqrt(delta))/(2*a);
		double t2 = (-b + sqrt(delta))/(2*a);

		if (t2 <0){
			return false;
		}

		double t;
		if (t1>0){
			t = t1;
		}
		else{
			t=t2;
		}

		P= r.C + t*r.u;
		N=(P - O);
		N.normalize();

		return true;

	}
		double R=R;
		Vector O=O;

};
 

int main() {
	int W = 512;
	int H = 512;
	Vector C(0,0,55);

	Sphere s(Vector(0,0,0),10);
	double fov=60*3.14/180;
	double tanfov2 = tan(fov/2);

	Vector lumiere_position(-10,20,40);
	double intensite_lumiere = 3000000;

	Vector rho(1,0,0); //albedo rouge (couleur de la sphère)

	std::vector<unsigned char> image2(W * H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			Vector u(j-W/2+0.5, H-i-H/2+0.5,-W/(2*tanfov2));
			u.normalize();
			Ray r(C,u);
			Vector P,N; 
			bool inter = s.intersect(r,P,N);

			double intensite_pixel=0.;

			Vector l = (lumiere_position - P);
			l.normalize();


			if (inter){
				intensite_pixel = intensite_lumiere * std::max(0.,dot(l,N)) / (4 * M_PI * (lumiere_position - P).norm2());
				
				image2[(i * W + j) * 3 + 0] = std::min(255.,std::max(0.,intensite_pixel))*rho[0];
				image2[(i * W + j) * 3 + 1] = std::min(255.,std::max(0.,intensite_pixel))*rho[1];
				image2[(i * W + j) * 3 + 2] = std::min(255.,std::max(0.,intensite_pixel))*rho[2];
			}
		

			//image[(i * W + j) * 3 + 0] = 255;
			//image[(i * W + j) * 3 + 1] = 0;
			//image[(i * W + j) * 3 + 2] = 0;
		}
	}
	stbi_write_png("image2.png", W, H, 3, &image2[0], 0);

	return 0;
}