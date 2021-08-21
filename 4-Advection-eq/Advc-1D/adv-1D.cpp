#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

using namespace std;


void cond_gauss(int l, double dx, double a, double x0, double sigma, vector<double>& vin) {
	int n = round (l/dx);
	for (int i=0; i<n; i++) {
		vin[i]=a*exp(-(i*dx-x0)*(i*dx-x0)/sigma);
		vin[n]=vin[0];
	}
}

double der1(int i, int n, double dx, vector<double>& vin) {
	if (i == 0) return (vin[1]-vin[n])/(2.0*dx);
	if (i == n) return (vin[0]-vin[n-1])/(2.0*dx);
	return (vin[i+1]-vin[i-1])/(2.0*dx);
}

double euler(int n, int j, double dt, double dx, double vel, vector<double>& q_n){
	return q_n[j] + -vel*dt*der1(j, n, dx, q_n);
}

vector<double> euler1(int n, double dt, double dx, double l, double vel, vector<double>& q_n) {
	vector<double> e(n+1);
	for(int k=0; k<=n; k++)
		e[k] = q_n[k] + -vel*dt*der1(k, n, dx, q_n);
	return e;
}

double Heun(int n, double l, int j, double dt, double dx, double vel, vector<double>& q_n){
	vector<double> e = euler1(n, dt, dx, l, vel, q_n);
	return q_n[j]+(-dt*0.5*vel)*(der1(j, n, dx, e) + der1(j, n, dx, q_n));
}

double Lax(int n, int j, double dt, double dx, double vel, vector<double>& q_n) {
	return (q_n[j+1]+q_n[j-1])/2.0+(-vel*dt)*der1(j, n, dx, q_n);
}

double Leapfrog(double dt, vector<double> q_n, vector<double> q_np, double vel, int j, int n, double dx){
	return q_np[j] + 2.0*(-vel*dt)*der1(j, n, dx, q_n);
}

void write_vector(ofstream& f, double dx, double l, vector<double>& vin) {
	int n = round (l/dx);
	for (int j=0; j<=n; j++) {
		double x=j*dx;
		f << x << " " << vin[j] << endl;
	}
	f << endl << endl;
}

void read_input(double& DX, double& L, double& DT, double& T, double& VEL, double& A, double& X0, double& SIGMA, int& NFOTOS, string filename) {
	ifstream filein(filename.c_str());
	if (filein.good() == false) {
		cout << "problem opening input file" << endl;
		exit(1);
	}
	filein >> DX;
	filein >> L;
	filein >> DT;
	filein >> T;
	filein >> VEL;
	filein >> A;
	filein >> X0;
	filein >> SIGMA;
	filein >> NFOTOS;

	filein.close();
}

int main() {
	int N, NT, NFOTOS, scheme;
	double T, L, DT, DX, VEL, A, X0, SIGMA;
	read_input(DX, L, DT, T, VEL, A, X0, SIGMA, NFOTOS, "input.inp");
	N = round (L/DX);
	NT = round (T/DT);
	vector<double> uin(N+1), u(N+1), unew(N+1), v(N+1), w(N+1);
	cond_gauss(L, DX, A, X0, SIGMA, uin);
	v=uin;
	u=uin;
	ofstream archivoout;

	cout << "which scheme" << endl;
	cout << "1 for Euler, 2 for Heun, 3 for lax y 4 for Leapfrog" << endl;
	cin >> scheme;

	if (scheme==1) 	archivoout.open ("FTCSeuler.dat");
	if (scheme==2) archivoout.open ("FTCSheun.dat");
	if (scheme==3) archivoout.open ("FTCSlax.dat");
	if (scheme==4) archivoout.open ("FTCSleapfrog.dat");

	// Esta es la condición inicial (el primer paso) para Leapfrog
	if (scheme == 4) {
		vector<double> e(N+1);
		double dtp = DT/NT;
		for (int itp=0; itp<=NT; itp++) {
			for(int k=0; k<=N; k++) {
				e[k] = Heun(N, L, k, dtp, DX, VEL, v);
				//v[k] - VEL*dtp*der1(k, N, DX, v); Esto sirve para hacer la condición inicial con euler
			}
			v=e;
		}
		w=e;
	}

	for (int it=0; it<NT; it++) {
			for (int j=0; j<=N; j++) {
				if (scheme==1) unew[j]=euler(N, j, DT, DX, VEL, u);
				if (scheme==2) unew[j]=Heun(N, L, j, DT, DX, VEL, u);
				if (scheme==3) unew[j]=Lax(N, j, DT, DX, VEL, u);
				if (scheme==4) unew[j]=Leapfrog(DT, w, v, VEL, j, N, DX);
			}
			if (it%(NT/NFOTOS)==0) {
				cout << "t = " << it*DT << endl;
				write_vector(archivoout, DX, L, unew);
			}
			u=unew;
			v=w;
			w=unew;

	}

	archivoout.close();

	return 0;
}
