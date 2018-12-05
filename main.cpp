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
	double alpha;
	
class Atom{
	public:	
		double x;
		double y;
		double z;
	
		double distance(Atom atom, vector<double> defMatrix){
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
			double tempX = dx;
			double tempY = dy;
			double tempZ = dz;
			
			dx = tempX*defMatrix[0] + tempY*defMatrix[1] + tempZ*defMatrix[2];
		    dy = tempX*defMatrix[3] + tempY*defMatrix[4] + tempZ*defMatrix[5];
			dz = tempX*defMatrix[6] + tempY*defMatrix[7] + tempZ*defMatrix[8];
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

vector<double> rijs (vector<Atom> atoms, vector<double> defMatrix){	
	vector<double> rij;
	for(int i=0; i < atoms.size(); i++){
		rij.push_back(atoms[0].distance(atoms[i], defMatrix));
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

double E_coh(vector<double> rij){
	double ebi = 0;
	double eri=0;
	for(int i=1; i < rij.size(); i++){
		ebi+=(eps*eps*exp(-2*q*(rij[i]/r0-1)));
		eri+=(A*exp(-p*((rij[i]/r0)-1)));
	}
	return eri + (-sqrt(ebi));
}

double b_calculate (vector<Atom> atoms, double Ecoh){	
	vector<double> stretchMatrix = {1+alpha, 0, 0,
									0, 1+alpha, 0,
									0, 0, 1+alpha};
	//vector<Atom>stretchAtoms = defAtoms(atoms, stretchMatrix);
	vector<double> str_rij = rijs(atoms, stretchMatrix);
	double Estr = E_coh(str_rij);
	vector<double> compMatrix = {1-alpha, 0, 0,
									0, 1-alpha, 0,
									0, 0, 1-alpha};
	//vector<Atom>compAtoms = defAtoms(atoms, compMatrix);
	vector<double> comp_rij = rijs(atoms, compMatrix);
	double Ecomp = E_coh(comp_rij);
	double B = 1.602*4*(Estr - 2*Ecoh + Ecomp)/(alpha*alpha*9*a0 * a0 * a0);
	return B;
}

vector<double> c11_c12_calculate(vector<Atom> atoms, double Ecoh){
		
	vector<double> c11stretchMatrix = {1+alpha, 0, 0,
								   0, 1+alpha, 0,
								   0, 0, 1};	
	vector<double> c11compMatrix = {1-alpha, 0, 0,
								   0, 1-alpha, 0,
								   0, 0, 1};
	double stretchEC11 = E_coh(rijs(atoms, c11stretchMatrix));
	double compEC11 = E_coh(rijs(atoms, c11compMatrix));
	
	vector<double> c12stretchMatrix = {1+alpha, 0, 0,
						   0, 1-alpha, 0,
						   0, 0, 1};
	vector<double> c12compMatrix = {1-alpha, 0, 0,
						   0, 1+alpha, 0,
						   0, 0, 1};
	double stretchEC12 = E_coh(rijs(atoms, c12stretchMatrix));
	double compEC12 = E_coh(rijs(atoms, c12compMatrix));
	double dEC11 = 1.602*(stretchEC11 - 2*Ecoh + compEC11)/(alpha*alpha*a0 * a0 * a0);
	double dEC12 = 1.602*(stretchEC12 - 2*Ecoh + compEC12)/(alpha*alpha*a0 * a0 * a0);
	
	double _eps = 0.001;
	double a[2][2] = {1,1,1,-1};
	vector<double> b = {dEC11,dEC12};
	vector<double> diag(2);
	vector<double> x(2);
	vector<double> x1 = {0,0};
	
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			if(i==j){
				diag[i]=a[i][j];
				a[i][j]=-b[i];
			}
		}
	}
	
	for(int i = 0; i < 2; i++){
		x[i] = b[i]/diag[i];
	}
	
	do{
		for(int i = 0; i < 2; i++){
			for(int j = 0; j < 2; j++){
				if(i==j){
					x1[i] += (a[i][j]/-diag[i]);
				}
				else{
					x1[i] += (a[i][j]*x[j]/-diag[i]);
				}
			}
		}
		
		double max = 0;
		for(int i = 0; i < 2; i++){
			if(max<fabs(fabs(x1[i])-fabs(x[i]))){
				max=fabs(fabs(x1[i])-fabs(x[i]));
			}
		}
		if(max<_eps){break;}
		x1.swap(x);
	}while(true);
	return x1;
}

double c44_calculate(vector<Atom> atoms, double Ecoh){
	vector<double> c44stretchMatrix = {1, alpha, 0,
									   alpha, 1, 0,
									   0, 0, 1/(1-alpha*alpha)};
									   
	vector<double> c44compMatrix = {1, -alpha, 0,
									   -alpha, 1, 0,
									   0, 0, 1/(1-alpha*alpha)};
									   
	double stretchEC44 = E_coh(rijs(atoms, c44stretchMatrix));
	double compEC44 = E_coh(rijs(atoms, c44compMatrix));
	double dEC44 = 1.602*(stretchEC44 - 2*Ecoh + compEC44)/(alpha*alpha*a0 * a0 * a0);
	return dEC44;
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
	alpha = j["alpha"];
	r0 = a0/sqrt(2);
	cout<<"A = "<<A<<endl;
	cout<<"eps = "<<eps<<endl;
	cout<<"p = "<<p<<endl;
	cout<<"q = "<<q<<endl;
	cout<<"a0 = "<<a0<<endl;
	cout<<"r0 = "<<r0<<endl;
	cout<<"n = "<<n<<endl;
	cout<<"alpha = "<<alpha<<endl;
	
	
	vector<Atom> atoms = atomsInit();
	vector<double> init = {1, 0, 0,
						   0, 1, 0,
						   0, 0, 1};
	vector<double> rij = rijs(atoms, init);
	double Ecoh = E_coh(rij);
	cout<<"Ecoh = "<<Ecoh<<endl;
	
	double B = b_calculate(atoms, Ecoh);
	cout<<"B = "<<B<<endl;
	
	vector<double> C11C12result = c11_c12_calculate(atoms, Ecoh);
	double C11 = C11C12result[0];
	double C12 = C11C12result[1];
	cout<<"C11 = "<<C11<<endl;
	cout<<"C12 = "<<C12<<endl;
	double C44 = c44_calculate(atoms, Ecoh);
	cout<<"C44 = "<<C44<<endl;	
}