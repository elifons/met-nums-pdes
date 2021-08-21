#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

using namespace std;


//------------------------------------------------------------------------------------------
double func(double x) {
	return sin(x);
}
//------------------------------------------------------------------------------------------

void init_vector(int n, double dx, vector<double>& vin) {
	for (int i=0; i<=n; i++) {
		double x=i*dx;
		vin[i]=func(x);
	}
}
//------------------------------------------------------------------------------------------

void der1(int n, double dx, int scheme, vector<double>& vin, vector<double>& vout) {
	vout[0]=(vin[1]-vin[0])/dx;
	vout[n]=(vin[n]-vin[n-1])/dx;
	for (int i=1; i<=n-1; i++) {
		double x=i*dx;
		//forward
		if (scheme==3) {
			vout[i]=(vin[i+1]-vin[i])/dx;
		}
		//backward
		if (scheme==2) {
			vout[i]=(vin[i]-vin[i-1])/dx;
		}
		//centered
		if (scheme==1) {
			vout[i]=(vin[i+1]-vin[i-1])/(2.0*dx);
		}
	}
}
//------------------------------------------------------------------------------------------

void der2(int n, double dx, int scheme, vector<double>& vin, vector<double>& vout) {
	//Borders
	vout[0]=(vin[2]-2.0*vin[1]+vin[0])/(dx*dx);
	vout[1]=(vin[3]-2.0*vin[2]+vin[1])/(dx*dx);
	vout[n-1]=(vin[n-1]-2.0*vin[n-2]+vin[n-3])/(dx*dx);
	vout[n]=(vin[n]-2.0*vin[n-1]+vin[n-2])/(dx*dx);
	for (int i=2; i<=n-2; i++) {
		double x=i*dx;
		//forward
		if (scheme==3) {
			vout[i]=(vin[i+2]-2.0*vin[i+1]+vin[i])/(dx*dx);
		}
		//backward
		if (scheme==2) {
			vout[i]=(vin[i]-2.0*vin[i-1]+vin[i-2])/(dx*dx);
		}
		//centered
		if (scheme==1) {
			vout[i]=(vin[i+1]+vin[i-1]-2.0*vin[i])/(dx*dx);
		}
	}
}

//------------------------------------------------------------------------------------------
void print_vector(int n, double dx, vector<double>& vin) {
	cout << "( ";
	for (int i=0; i<n; i++)
		cout << vin[i] << ", ";
	cout << vin[n] << ")";
	cout << endl;
}

//------------------------------------------------------------------------------------------
void write_vector(int n, double dx, string filename, vector<double>& vin) {
	ofstream archivoout(filename.c_str());
	for (int i=0; i<=n; i++) {
		double x=i*dx;
		archivoout << x << " " << vin[i] << endl;
	}
	archivoout.close();
}
//------------------------------------------------------------------------------------------
void read_input(int& N, double& DX, int& SCHEME, int& ORD, string filename) {
	string tmpstr;
	ifstream filein(filename.c_str());

	if (!filein.good()) {
		cout << "can't open file " << filename << endl;
		exit(1);
	}

	filein >> N;
	filein >> DX;
	filein >> SCHEME;
	filein >> ORD;

	filein.close();
}
//------------------------------------------------------------------------------------------
int main() {
	int N, SCHEME, ORD;
	double L, DX;

	read_input(N, DX, SCHEME, ORD, "derive.inp");
	vector<double> grid(N+1);
	vector<double> dgrid(N+1);
	vector<double> d2grid(N+1);
	init_vector(N,DX,grid);
	der1(N,DX,SCHEME,grid,dgrid);
	der2(N,DX,SCHEME,grid,d2grid);
	write_vector(N,DX,"der.dat",dgrid);
	write_vector(N,DX,"func.dat",grid);
	write_vector(N,DX,"der2.dat",d2grid);

	return 0;
}
