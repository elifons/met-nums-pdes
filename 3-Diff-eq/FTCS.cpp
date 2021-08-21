#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

void cond(int l, double dx, vector<double>& vin) {
	int n = round (l/dx);
	for (int i=0; i<n; i++) {
		vin[i]=0;
	}
	vin[n]=100;
}

double der2(int i, double dx, vector<double>& vin) {
	return (vin[i+1]+vin[i-1]-2.0*vin[i])/(dx*dx);
}

double euler(int j, double dt, double dx, double alpha, vector<double>& q_n){
	return q_n[j] + alpha*dt*der2(j, dx, q_n);
}

vector<double> euler1(double dt, double dx, double l, double alpha, vector<double>& q_n) {
	int n = round (l/dx);
	vector<double> e(n+1);
	e[0]=0;
	e[n]=100;
	for(int k=1; k<n; k++)
		e[k] = q_n[k] + alpha*dt*der2(k, dx, q_n);
	return e;
}

double Heun(double l, int j, double dt, double dx, double alpha, vector<double>& q_n){
	vector<double> e = euler1(dt, dx, l, alpha, q_n);
// 	return q_n[j]+(dt*0.5*alpha)*(func(q_n+func(q_n, gamma)*dt,gamma)+func(q_n, gamma));
	return q_n[j]+(dt*0.5*alpha)*(der2(j, dx, e) + der2(j, dx, q_n));
}
void write_vector(ofstream& f, double dx, double l, vector<double>& vin) {
	int n = round (l/dx);
	for (int j=0; j<=n; j++) {
		double x=j*dx;
		f << x << " " << vin[j] << endl;
	}
	f << endl;
}

void read_input(double& DX, double& L, double& DT, double& T, double& ALPHA, int& NFOTOS, string filename) {
	ifstream filein(filename.c_str());
	if (filein.good() == false) {
		cout << "problem opening file" << endl;
		exit(1);
	}
	filein >> DX;
	filein >> L;
	filein >> DT;
	filein >> T;
	filein >> ALPHA;
	filein >> NFOTOS;

	filein.close();
}

int main() {
	int N, NT, NFOTOS, scheme;
	double T, L, DT, DX, ALPHA;
	read_input(DX, L, DT, T, ALPHA, NFOTOS, "input.inp");
	N = round (L/DX);
	NT = round (T/DT);
	vector<double> uin(N+1), u(N+1), unew(N+1);
// 	double DT=T/NT;
// 	double DX=L/N;
	cond(L, DX, uin);
	u=uin;
	ofstream archivoout;

	cout << "which scheme?" << endl;
	cout << "1 for Euler y 2 for Heun" << endl;
	cin >> scheme;

	if (scheme==1) 	archivoout.open ("FTCSeuler.dat");
	if (scheme==2) archivoout.open ("FTCSheun.dat");

	for (int it=0; it<NT; it++) {
		for (int j=1; j<=N-1; j++) {
// 			u[0]=0;
// 			u[N]=100;
			unew[N]=100;
			if (scheme==1) unew[j]=euler(j, DT, DX, ALPHA, u);
			if (scheme==2) unew[j]=Heun(L, j, DT, DX, ALPHA, u);
		}
		if (it%(NT/NFOTOS)==0) write_vector(archivoout, DX, L, unew);
		u=unew;
	}

	archivoout.close();

	return 0;
}
