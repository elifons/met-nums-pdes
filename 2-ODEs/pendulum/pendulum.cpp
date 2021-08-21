#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

double func(double y, double omega) {
	return (-omega*omega)*sin(y);
}

void esq(int n, double dt, double v0, double q0, double omega, vector<double>& u, vector<double>& v) {
	u[0]=q0;
	v[0]=v0;
	for (int it=0; it<n; it++) {
		v[it+1] = v[it] + dt*func(u[it], omega);
		u[it+1] = u[it] + dt*v[it];
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

void read_input(int& N, double& DT, double& V0, double& Q0, double& OMEGA, string filename) {
	string tmpstr;
	ifstream filein(filename.c_str());

	if (!filein.good()) {
		cout << "problem opening the file " << filename << endl;
		exit(1);
	}

	filein >> N;
	filein >> DT;
	filein >> V0;
	filein >> Q0;
	filein >> OMEGA;

	filein.close();
}

int main() {
	int N;
	double DT, V0, Q0, OMEGA;

	read_input(N, DT, V0, Q0, OMEGA, "cond.inp");

	vector<double> u(N+2), v(N+2);

	esq(N, DT, V0, Q0, OMEGA, u, v);
	write_vector(N, DT, "vel.dat", v);
	write_vector(N, DT, "pos.dat", u);

	return 0;
}
