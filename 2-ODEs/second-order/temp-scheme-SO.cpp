#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

double round(double t, double dt) {
		return floor(t/dt + 0.5);
		}


double func(double y, double omega) {
	return (-omega*omega)*y;
}

void esq(int n, double dt, double v0, double q0, double omega, double m, vector<double>& u, vector<double>& v, vector<double>& en) {
	u[0]=q0;
	v[0]=v0;
	en[0]=0.5*m*v0*v0+0.5*m*omega*omega*q0;
	for (int it=0; it<n; it++) {
		v[it+1] = v[it] + dt*func(u[it], omega);
		u[it+1] = u[it] + dt*v[it];
		en[it+1] = 0.5*m*v[it+1]*v[it+1]+0.5*m*omega*omega*u[it+1]*u[it+1];
	}
}

void write_vector(int n, double dt, string filename, vector<double>& vin) {
	ofstream archivoout(filename.c_str());
	for (int i=0; i<=n; i++) {
		double t=i*dt;
		archivoout << t << " " << vin[i] << endl;
	}
	archivoout.close();
}

void read_input(double& T, double& DT, double& V0, double& Q0, double& OMEGA, double& M, string filename) {
	string tmpstr;
	ifstream filein(filename.c_str());

	if (!filein.good()) {
		cout << "problem opening file " << filename << endl;
		exit(1);
	}

	filein >> T;
	filein >> DT;
	filein >> V0;
	filein >> Q0;
	filein >> OMEGA;
	filein >> M;

	//cierro el archivo
	filein.close();
}

int main() {
	int N;
	double DT, V0, Q0, OMEGA, M, T;

	read_input(T, DT, V0, Q0, OMEGA, M, "cond.inp");

	N=round(T, DT);
	vector<double> u(N+2), v(N+2), en(N+2);

	esq(N, DT, V0, Q0, OMEGA, M, u, v, en);
	write_vector(N, DT, "vel.dat", v);
	write_vector(N, DT, "pos.dat", u);
	write_vector(N, DT, "energy.dat", en);

	return 0;
}
