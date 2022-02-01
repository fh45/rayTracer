#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <list>

#include <iostream>




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
	Sphere( const Vector& O, const double R, const Vector& rho, const bool &reflect) : O(O),R(R),rho(rho),reflect(reflect){

};

	

	bool intersect(const Ray& r , Vector& P, Vector& N, double& t) {
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
		Vector rho=rho;
		bool reflect=reflect;

};

class Scene{
	public:
	Scene(){};
	void addSphere(const Sphere& s) {spheres.push_back(s);}

	bool intersect_scene(const Ray& r , Vector& P, Vector& N, int &rang_sphere, double &t_min){

		bool inter_scene = false;
		t_min=1E99;

		for (int i=0;i<spheres.size();i++){
			Vector P_sphere, N_sphere;
			double t;
			bool inter_sphere = spheres[i].intersect(r,P_sphere,N_sphere,t);
			if(inter_sphere){
				inter_scene=true;
				if(t<t_min){
					t_min=t;
					P = P_sphere;
					N= N_sphere;
					rang_sphere=i;
				}
			}
		}

		return inter_scene;

	};


	Vector getColor(Ray &r, int rebond, Vector &lumiere_position, double &intensite_lumiere){
	
	if(rebond<=0){
		return Vector(0,0,0);
	}

	else{
		Vector P,N;
		int rang_sphere;
		double t_min;

		bool inter = intersect_scene(r,P,N,rang_sphere,t_min);

		double intensite_pixel=0.;
		Vector l = (lumiere_position - P);
		l.normalize(); 

		if (inter){

			if(spheres[rang_sphere].reflect){

				//std::cout<<"reflet";
				Vector vector_reflect = r.u - 2*dot(r.u,N)*N;
				Ray rayon_reflect = Ray(P+N*0.001,vector_reflect);
				return getColor(rayon_reflect,rebond-1, lumiere_position, intensite_lumiere);
			}

			else{

				//std::cout<<"inter";

				Vector P_lum,N_lum;
				double t_min_lum;
				int rang_sphere_lum;


				if(intersect_scene(Ray(P+0.001*N,l),P_lum,N_lum,rang_sphere_lum,t_min_lum))
				{
					
					if (t_min_lum < sqrt((lumiere_position - P).norm2())){
						intensite_pixel=0.;
						
					}

					else{
						intensite_pixel = intensite_lumiere * std::max(0.,dot(l,N)) / (4 * M_PI * (lumiere_position - P).norm2());
						
					}
				}
				else{

					intensite_pixel = intensite_lumiere * std::max(0.,dot(l,N)) / (4 * M_PI * (lumiere_position - P).norm2());
				
				}
				
				//intensite_pixel = intensite_lumiere * std::max(0.,dot(l,N)) / (4 * M_PI * (lumiere_position - P).norm2());
				intensite_pixel = pow(intensite_pixel,1/2.2);
				
				double a = std::min(255.,std::max(0.,intensite_pixel))*spheres[rang_sphere].rho[0];
				double b = std::min(255.,std::max(0.,intensite_pixel))*spheres[rang_sphere].rho[1];
				double c = std::min(255.,std::max(0.,intensite_pixel))*spheres[rang_sphere].rho[2];

				return Vector(a,b,c);


			}
		
		}

		else{
			return Vector(0,0,0);
			
		}

	}
};

	std::vector<Sphere> spheres;
	

};





 

int main() {
	int W = 512;
	int H = 512;
	Vector C(0,0,55);

	Sphere s1(Vector(0,10,0),5,Vector(0,1,1),true); //albedo vert (couleur de la sphère)
	//Sphere s2(Vector(20,20,20),5,Vector(1,0,0),false); //albedo rouge (couleur de la sphère)

	
	Sphere s3(Vector(0,0,-1000),1000-80,Vector(1,0,1),false); //mur du fond 
	Sphere s5(Vector(0,-1000,0),1000-10,Vector(1,1,0),false); //sol  
	Sphere s6(Vector(1000,0,0),1000-40,Vector(0,0,1),false); //mur de droite
	Sphere s7(Vector(-1000,0,0),1000-50,Vector(0,1,1),false); //mur de gauche
	Sphere s8(Vector(0,1000,0),1000-50,Vector(1,0,0),false); //plafond

	Scene scene;
	scene.addSphere(s1);
	//scene.addSphere(s2);
	scene.addSphere(s3);
	//scene.addSphere(s4);
	scene.addSphere(s5);
	scene.addSphere(s6);
	scene.addSphere(s7);
	scene.addSphere(s8);

	//std::cout<<scene.spheres[1].reflect;

	double fov=60*3.14/180;
	double tanfov2 = tan(fov/2);

	Vector lumiere_position(-20,40,80);
	double intensite_lumiere = 8000000000;

	int nb_rebond_max = 10;//scene.spheres.size();

	std::vector<unsigned char> image2(W * H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			Vector u(j-W/2+0.5, H-i-H/2+0.5,-W/(2*tanfov2));
			u.normalize();
			Ray r(C,u);

			/*
			Vector P,N;
			int rang_sphere;
			double t_min;

			bool inter = scene.intersect_scene(r,P,N,rang_sphere,t_min);

			double intensite_pixel=0.;

			Vector l = (lumiere_position - P);
			l.normalize();

			
			if (inter){ 

				Vector P_lum,N_lum;
				double t_min_lum;
				int rang_sphere_lum;

				if(scene.intersect_scene(Ray(P+0.001*N,l),P_lum,N_lum,rang_sphere_lum,t_min_lum))
				{
					
					if (t_min_lum < (lumiere_position - P).norm2()){
						intensite_pixel=0.;
						
					}

					else{
						intensite_pixel = intensite_lumiere * std::max(0.,dot(l,N)) / (4 * M_PI * (lumiere_position - P).norm2());
					}
				}
				else{
					intensite_pixel = intensite_lumiere * std::max(0.,dot(l,N)) / (4 * M_PI * (lumiere_position - P).norm2());
				}
				
			
				intensite_pixel = pow(intensite_pixel,1/2.2);


				
				image2[(i * W + j) * 3 + 0] = std::min(255.,std::max(0.,intensite_pixel))*scene.spheres[rang_sphere].rho[0];
				image2[(i * W + j) * 3 + 1] = std::min(255.,std::max(0.,intensite_pixel))*scene.spheres[rang_sphere].rho[1];
				image2[(i * W + j) * 3 + 2] = std::min(255.,std::max(0.,intensite_pixel))*scene.spheres[rang_sphere].rho[2];
			}*/
		
			Vector pixel = scene.getColor(r,nb_rebond_max,lumiere_position,intensite_lumiere);
			image2[(i * W + j) * 3 + 0] = pixel[0];
			image2[(i * W + j) * 3 + 1] = pixel[1];
			image2[(i * W + j) * 3 + 2] = pixel[2];
		}
	}
	stbi_write_png("image2.png", W, H, 3, &image2[0], 0);

	return 0;
}