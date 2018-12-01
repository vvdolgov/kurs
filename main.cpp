#include <iostream>
#include <nlohmann/json.hpp>
#include <fstream>
#include <math.h>

using json = nlohmann::json;
using namespace std;
	double A;
	double eps;
	double p;
	double q;
	double a0;
	double r0;
	int n;
	
class Atom{
	
	public:	
		double x;
		double y;
		double z;
	
		double distance(Atom atom){
			double dx = atom.x - x;
			double dy = atom.y - y;
			double dz = atom.z - z;
			if(dx > a0 * 1.5){
				dx -= 3 * a0;
			}
			else if(dx < - a0 * 1.5){
				dx += 3 * a0;
			}
			if(dy > a0 * 1.5){
				dy -= 3 * a0;
			}
			else if(dy < - a0 * 1.5){
				dy += 3 * a0;
			}
		    if(dz > a0 * 1.5){
				dz -= 3 * a0;
			}
			else if(dz < - a0 * 1.5){
				dz += 3 * a0;
			}
			return sqrt((dx)*(dx)+(dy)*(dy)+(dz)*(dz));
		}
		
		Atom(double _x, double _y, double _z){
			x = _x;
			y = _y;
			z = _z;
		}
		
		void print(){
			cout<<"Atom positions: {"<<x<<","<<y<<","<<z<<"}"<<endl;
		}
};

vector<double> rijs (vector<Atom> atoms){
	vector<double> rij;
	for(int i=0; i < atoms.size(); i++){
		rij.push_back(atoms[0].distance(atoms[i]));
	}
	return rij;
}
vector<Atom> atomsInit (){
	vector<Atom> atoms;
	for(int i=0; i < n; i++){
		for(int j=0; j < n; j++){
			for(int k=0; k < n; k++){
					Atom atom(i*a0, j*a0, k*a0);
					Atom atom1(i*a0, (j+0.5)*a0, (k+0.5)*a0);
					Atom atom2((i+0.5)*a0, j*a0, (k+0.5)*a0);
					Atom atom3((i+0.5)*a0, (j+0.5)*a0, k*a0);
					atoms.push_back(atom);
					atoms.push_back(atom1);
					atoms.push_back(atom2);
					atoms.push_back(atom3);
			}
		}
	}
	return atoms;
}

double ERi(vector<double> rij){
	double eri;
	for(int i=1; i < rij.size(); i++){
		eri+=(A*exp(-p*((rij[i]/r0)-1)));
	}
	return eri;
}

double EBi(vector<double> rij){
	double ebi = 0;
	for(int i=1; i < rij.size(); i++){
		ebi+=(eps*eps*exp(-2*q*(rij[i]/r0-1)));
	}
	return -sqrt(ebi);
}

int main(int argc, char* argv[]){
	json j;
	ifstream file ("resources/param.json");
	file>>j;
	A = j["A"];
	eps = j["eps"];
	p = j["p"];
	q = j["q"];
	a0 = j["a0"];
	n = j["n"];
	r0 = a0/sqrt(2);
	cout<<"A = "<<A<<endl;
	cout<<"eps = "<<eps<<endl;
	cout<<"p = "<<p<<endl;
	cout<<"q = "<<q<<endl;
	cout<<"a0 = "<<a0<<endl;
	cout<<"r0 = "<<r0<<endl;
	int n = 3;
	vector<Atom> atoms = atomsInit();
	vector<double> rij = rijs(atoms);
	double eri = ERi(rij);
	double ebi = EBi(rij);
	double Ecoh = eri+ebi;
	cout<<"Ecoh = "<<Ecoh<<endl;
}