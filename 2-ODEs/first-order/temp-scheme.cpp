#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

double func(double x, double gamma) {
	return -gamma*x;
}

double euler(double dt, double q_n, double gamma) {
	return q_n + dt*func(q_n, gamma);
}

double Heun(double dt, double q_n, double gamma){
	return q_n+(dt*0.5)*(func(q_n+func(q_n, gamma)*dt,gamma)+func(q_n, gamma));
}

double Leapfrog(double dt, double q_n, double q_n1, double gamma){
	return q_n1+2*dt*func(q_n, gamma);
}

// double AdamsBashforth(double dt, double q_n, double q_n1, double gamma){
// 	return q_n+(dt*0.5)*3*(func(q_n, gamma) + func (q_n1, gamma));
// }

void esq(int n, double dt, double q0, double gamma, int scheme, vector<double>& u) {

	double dtp = dt/n;
	vector<double> w(n+2);
	w[0]=u[0]=q0;
	for (int it=0; it<n; it++){
		w[it+1] = euler(dtp, w[it], gamma);
	}
	u[1]=w[n];

	for (int it=1; it<n; it++){
		if (scheme==1) { //Euler
			u[it+1] = euler(dt, u[it], gamma);
		}
		if (scheme==2) { //Heun
			u[it+1] = Heun(dt, u[it], gamma);
		}
		if (scheme==3) { // Leapfrog
			u[it+1] = Leapfrog(dt, u[it], u[it-1], gamma);
		}
// 		if (scheme==4) { //Adam-Bashforth
// 			u[it+1] = AdamBashforth(dt, u[it], u[it-1], gamma);
// 		}
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

void read_input(int& N, double& DT, double& Q0, double& GAMMA, int& SCHEME, string filename) {
	string tmpstr;

	ifstream filein(filename.c_str());


	if (!filein.good()) {
		cout << "cannot open file " << filename << endl;
		exit(1);
	}

	filein >> N;
	filein >> DT;
	filein >> Q0;
	filein >> GAMMA;
	filein >> SCHEME;

	filein.close();
}
int main() {
	int N, SCHEME;
	double DT, Q0, GAMMA;

	read_input(N, DT, Q0, GAMMA, SCHEME, "cond.inp");
	vector<double> u(N+2), w(N+2);

	esq(N, DT, Q0, GAMMA, 1, u);
	write_vector(N, DT, "Euler.dat", u);

	esq(N, DT, Q0, GAMMA, 2, u);
	write_vector(N, DT, "Heun.dat", u);

	esq(N, DT, Q0, GAMMA, 3, u);
	write_vector(N, DT, "Leapfrog.dat", u);

return 0;

}
