#include <iostream>
#include <nlohmann/json.hpp>
#include <fstream>
#include <math.h>

using json = nlohmann::json;
using namespace std;
class Atom{
	
	public:	
		double x;
		double y;
		double z;
	
		double distance(Atom atom, double a0){
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

vector<double> rijs (double a0, int n){
	Atom atom(0, 0, 0);
	vector<double> rij;
	for(int i=0; i < n; i++){
		for(int j=0; j < n; j++){
			for(int k=0; k < n; k++){
					Atom atom2(i*a0, j*a0, k*a0);
					Atom atom3(i*a0, (j+0.5)*a0, (k+0.5)*a0);
					Atom atom4((i+0.5)*a0, j*a0, (k+0.5)*a0);
					Atom atom5((i+0.5)*a0, (j+0.5)*a0, k*a0);
					rij.push_back(atom.distance(atom2, a0));
					rij.push_back(atom.distance(atom3, a0));
					rij.push_back(atom.distance(atom4, a0));
					rij.push_back(atom.distance(atom5,a0));
			}
		}
	}
	return rij;
}

double ERi(double r0, double p, double A, vector<double> rij, int n){
	double eri;
	for(int i=1; i < rij.size(); i++){
		eri+=(A*exp(-p*((rij[i]/r0)-1)));
	}
	return eri;
}

double EBi(double r0, double q, double eps, vector<double> rij, int n){
	double ebi = 0;
	for(int i=1; i < rij.size(); i++){
		ebi+=(eps*eps*exp(-2*q*(rij[i]/r0-1)));
	}
	return -sqrt(ebi);
}

// double Ecoh(vector<double> eri, vector<double> ebi, int n){
	// double ecoh = 0;
	// for(int i=1; i < ebi.size(); i++){
		// ecoh+=eri[i];
		// ecoh+=ebi[i];
	// }
	// return ecoh;
// }
int main(int argc, char* argv[]){
	json j;
	ifstream file ("resources/param.json");
	file>>j;
	double A = j["A"];
	double eps = j["eps"];
	double p = j["p"];
	double q = j["q"];
	double a0 = j["a0"];
	double r0 = a0/sqrt(2);
	
	cout<<"A = "<<A<<endl;
	cout<<"eps = "<<eps<<endl;
	cout<<"p = "<<p<<endl;
	cout<<"q = "<<q<<endl;
	cout<<"a0 = "<<a0<<endl;
	cout<<"r0 = "<<r0<<endl;
	int n = 3;
	vector<double> rij = rijs(a0, n);
	double eri = ERi(r0, p, A, rij, n);
	double ebi = EBi(r0, q, eps, rij, n);
	double Ecoh = eri+ebi;
	cout<<"Ecoh = "<<Ecoh<<endl;
}