#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

typedef vector<vector<double> > matrix;

matrix init_matrix(int N) {
	matrix v;
	v.resize(N);
	for (int j=0; j<N; j++)
		v[j].resize(N);
	return v;
}

void cond_gauss(int l, double dx, double dy, double a, double x0, double y0, double sigma, matrix& vin) {
	int n = round(l/dx);
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
		vin[i][j]=a*exp(-pow(i*dx-x0, 2)/sigma)*exp(-pow(j*dy-y0, 2)/sigma);
		}
	}
}

double der1x(int i, int j, int n, double dx, matrix& vin) {
	if (i == 0) return (vin[1][j]-vin[n-1][j])/(2.0*dx);
	if (i == n-1) return (vin[0][j]-vin[n-2][j])/(2.0*dx);
	return (vin[i+1][j]-vin[i-1][j])/(2.0*dx);
}
double der1y(int i, int j, int n, double dy, matrix& vin) {
	if (j == 0) return (vin[i][1]-vin[i][n-1])/(2.0*dy);
	if (j == n-1) return (vin[i][0]-vin[i][n-2])/(2.0*dy);
	return (vin[i][j+1]-vin[i][j-1])/(2.0*dy);
}

double euler(int n, int j, int i, double dt, double dx, double dy, double velx, double vely, matrix& q_n){ //arreglar parámetros de entrada
	return q_n[i][j] -velx*dt*der1x(i, j, n, dx, q_n) - vely*dt*der1y(i, j, n, dy, q_n);
}

void write_vector(ofstream& f, double dx, double dy, double l, matrix& m) {
	int n = round (l/dx);
	for(int i=0; i<n; i++) {
		double x = i*dx;
		for (int j=0; j<n; j++) {
			double y=j*dy;
			f << x << " " << y << " " << m[j][i] << endl;
		}
	}
	f << endl << endl;
}

void read_input(double& DX, double& LX, double& DY, double& LY, double& DT, double& T, double& VELX, double& VELY, double& A, double& X0, double& Y0, double& SIGMA, int& NFOTOS, string filename) {
	ifstream filein(filename.c_str());
	if (filein.good() == false) {
		cout << "problem opening input file" << endl;
		exit(1);
	}
	filein >> DX;
	filein >> LX;
	filein >> DY;
	filein >> LY;
	filein >> DT;
	filein >> T;
	filein >> VELX;
	filein >> VELY;
	filein >> A;
	filein >> X0;
	filein >> Y0;
	filein >> SIGMA;
	filein >> NFOTOS;

	filein.close();
}

int main() {
	int N, NT, NFOTOS;
	double T, LX, LY, DT, DX, DY, VELX, VELY, A, X0, Y0, SIGMA;
	read_input(DX, LX, DY, LY, DT, T, VELX, VELY, A, X0, Y0, SIGMA, NFOTOS, "input.inp");
	N = round (LX/DX);
	NT = round (T/DT);
	matrix uin, u, unew;
	uin = init_matrix(N);
	u = init_matrix(N);
	unew = init_matrix(N);

	cond_gauss(LX, DX, DY, A, X0, Y0, SIGMA, uin);
	u=uin;
	ofstream archivoout;

	archivoout.open ("FTCSeuler.dat");

	for (int it=0; it<NT; it++) {
			for (int i=0; i<N; i++) {
				for (int j=0; j<N; j++) {
					unew[i][j]=euler(N, j, i, DT, DX, DY, VELX, VELY, u);
				}
			}
			if (it%(NT/NFOTOS)==0) {
				cout << "t = " << it*DT << endl;
				write_vector(archivoout, DX, DY, LX, unew);
			}
			u=unew;
	}

	archivoout.close();

	return 0;
}
