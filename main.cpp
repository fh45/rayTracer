/*VERSION DU 07/02

---Fonctionnalités---
x Affichage d'une sphère
x Affichage d'une scène de sphères avec des couleurs différentes
x Eclairage de la sphère par une source
x Ombre portée des sphères
X Surface diffuse
X Surface Speculaire 


-----Bugs-----
x Surfaces transparentes

*/


#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <list>
#include <iostream>

#include <random>

std::default_random_engine engine; 
std::uniform_real_distribution<double> uniform(0,1);

//DEFINITION DES CLASSES


//VECTOR (X,Y,Z)
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

Vector operator*(const Vector& a, const Vector& b){
	return Vector(a[0]*b[0],a[1]*b[1],a[2]*b[2]);
}




//RAYON (ORIGINE,DIRECTION)
class Ray {
public:
	Ray(const Vector& origin, const Vector&  direction) : C(origin),u(direction){

};

Vector C,u;

};
 

//SPHERE (CENTRE, RAYON, ALBEDO, MIROIR, TRANSPARENCE) 
class Sphere {
public:
	Sphere( const Vector& O, const double R, const Vector& rho, const bool &reflect, const bool &transparent) : O(O),R(R),rho(rho),reflect(reflect),transparent(transparent){

};

	
	//Routine d'intersection qui teste si un rayon touche la sphère
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
		bool transparent; 

};

// SCENE(SPEHRES)
class Scene{
	public:
	Scene(){};

	// AJOUT D'UNE SPHERE DANS LA SCENE
	void addSphere(const Sphere& s) {spheres.push_back(s);}


	// ROUTINE D'INTERSECTION ENTRE UN RAYON ET TOUTES LES SPHERES DE LA SCENE
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

	}


	// RECUPERE LA COULEUR RENVOYEE PAR UN RAYON
	Vector getColor(Ray &r, int rebond){

	Vector intensite_pixel(0,0,0); //On initialise le pixel en noir
	
	if(rebond<=0){ //Nombre max de rebond dépassé
		return intensite_pixel;
	}

	else{ //Nombre de rebond max non dépassé
		
		Vector P,N;
		int rang_sphere;
		double t_min;
		bool inter = intersect_scene(r,P,N,rang_sphere,t_min);

		if (inter){ //le rayon intersecte avec un sphere de la scène

			Vector l = (lumiere_position - P);
			l.normalize();

			// LA SPHERE EST MIROIR
			if(spheres[rang_sphere].reflect){

				Vector vector_reflect = r.u - 2*dot(r.u,N)*N;
				Ray rayon_reflect = Ray(P+N*0.001,vector_reflect);
				return getColor(rayon_reflect,rebond-1);
			}

			else{

				// .A SPHERE EST TRANSPARENTE
				if(spheres[rang_sphere].transparent){

					double n1 = 1;
					double n2 = 1.3;
					Vector N_transp = N;
					Vector dir = r.u;
					if(dot(r.u,N)>=0){ // Rayon sortant de la sphère, inversion de n1,n2. 
						n1 = 1.3;
						n2 = 1;
						N_transp = Vector(0,0,0)-N;
					}

					double racine = 1-(n1/n2)*(n1/n2)*(1-(dot(N_transp,dir)*(dot(N_transp,dir))));
					if(racine > 0){
						Vector vector_ref = (n1/n2)*(dir - dot(dir,N_transp)*N_transp) - N_transp*sqrt(racine);
						Ray rayon_ref = Ray(P ,vector_ref);
						return getColor(rayon_ref,rebond-1);
						
					}
					
					return Vector(0,0,0);
				}


				else{ //LA SPHERE EST QUELCONQUE

				// ECLAIRAGE DIRECT

				//On envoie un rayon vers la source de lumière pour voir si on a une ombre (objet de la scène rencontré)
				Vector P_lum,N_lum;
				double t_min_lum;
				int rang_sphere_lum;
				if(intersect_scene(Ray(P+0.001*N,l),P_lum,N_lum,rang_sphere_lum,t_min_lum))
				{
					
					if (t_min_lum < sqrt((lumiere_position - P).norm2())){ //objet devant la lumière donc ombre
						intensite_pixel = Vector(0,0,0);
						
					}

					else{ //objet derrière = pas d'ombre
						intensite_pixel = spheres[rang_sphere].rho * lumiere_intensite *dot(l,N) / (4 * M_PI * (lumiere_position - P).norm2());
						
					}
				}

				else{ //pas d'objet rencontré = pas d'ombre

					//intensite_pixel = intensite_lumiere * std::max(0.,dot(l,N)) / (4 * M_PI * (lumiere_position - P).norm2());
					intensite_pixel = spheres[rang_sphere].rho * lumiere_intensite * dot(l,N) / (4 * M_PI * (lumiere_position - P).norm2());
				
				}

				
				// ECLAIRAGE INDIRECT
				
				double r1 = uniform(engine);
				double r2 = uniform(engine);

				Vector direction_random_local(cos(2*M_PI*r1)*sqrt(1-r2),sin(2*M_PI*r1)*sqrt(1-r2),sqrt(r2))  ;
				Vector vector_random(uniform(engine)-0.5,uniform(engine)-0.5,uniform(engine)-0.5);
				Vector tangente1 = cross(N,vector_random); tangente1.normalize();
				Vector tangente2= cross(tangente1,N);

				Vector direction_random = direction_random_local[2]*N + direction_random_local[0]*tangente1 + direction_random_local[1]*tangente2;
				Ray rayon_random(P + 0.001*N, direction_random)	;			

				return intensite_pixel+getColor(rayon_random,rebond-1)*spheres[rang_sphere].rho;

			}
			}
		
		}

		else{ //Le rayon n'intersecte avec aucune sphère
			//std::cout<<"test";
			return intensite_pixel;
			
		}

	}
};

	std::vector<Sphere> spheres;
	Vector lumiere_position;
	double lumiere_intensite;

};





 

// MAIN 
int main() {

	// RESOLUTION DE L'IMAGE
	int W = 1024;
	int H = 1024;

	//Monte Carlo
	const int Nb_rayons=16;

	//NOM DE L'IMAGE
	const char* nom_image = "eclairage_indirect_8_rayons.png";

	// POSITION DE LA CAMERA
	Vector C(0,0,55);

	// CREATION DE LA SCENE

	Sphere s1(Vector(0,10,0),5,Vector(0,1,1),false,false); //albedo vert (couleur de la sphère)
	//Sphere s2(Vector(20,20,20),5,Vector(1,0,0),false); //albedo rouge (couleur de la sphère)
	Sphere s3(Vector(0,0,-1000),1000-80,Vector(1,0,1),false,false); //mur du fond 
	Sphere s5(Vector(0,-1000,0),1000-10,Vector(1,1,0),false,false); //sol  
	Sphere s6(Vector(1000,0,0),1000-40,Vector(0,0,1),false,false); //mur de droite
	Sphere s7(Vector(-1000,0,0),1000-50,Vector(0,1,1),false,false); //mur de gauche
	Sphere s8(Vector(0,1000,0),1000-50,Vector(1,0,0),false,false); //plafond

	Scene scene;

	scene.addSphere(s1);
	//scene.addSphere(s2);
	scene.addSphere(s3);
	//scene.addSphere(s4);
	scene.addSphere(s5);
	scene.addSphere(s6);
	scene.addSphere(s7);
	scene.addSphere(s8);

	// NOMBRE DE REBONDS MAX POUR LES SPHERES MIROIRS
	int nb_rebond_max = 10; //scene.spheres.size();


	// ANGLE DE CHAMPS
	double fov=60*3.14/180;
	double tanfov2 = tan(fov/2);

	
	// LUMIERE POSITION ET INTENSITE
	scene.lumiere_position = Vector(-20,40,80);
	scene.lumiere_intensite = 8000000000;


	// CALCUL DE L'IMAGE
	std::vector<unsigned char> image(W * H * 3, 0);

#pragma cmp parallel for
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			//CREATION DU RAYON POUR CHAQUE PIXEL
			Vector u(j-W/2+0.5, H-i-H/2+0.5,-W/(2*tanfov2));
			u.normalize();
			Ray r(C,u);

			//RECUPERATION DE LA COULEUR DU RAYON
			
			Vector pixel(0.,0.,0.);
			for (int k=0; k<Nb_rayons;k++)
				pixel = pixel + scene.getColor(r,nb_rebond_max)/Nb_rayons;

			//Vector intensite_pixel = scene.getColor(r,nb_rebond_max);

			//STOCKAGE DANS L'IMAGE
			image[(i * W + j) * 3 + 0] = std::min(255.,std::max(0., std::pow(pixel[0],1/2.2)));
			image[(i * W + j) * 3 + 1] = std::min(255.,std::max(0.,std::pow(pixel[1],1/2.2)));
			image[(i * W + j) * 3 + 2] = std::min(255.,std::max(0.,std::pow(pixel[2],1/2.2)));


			
			
			
		}
	}
	stbi_write_png(nom_image, W, H, 3, &image[0], 0);

	return 0;
}