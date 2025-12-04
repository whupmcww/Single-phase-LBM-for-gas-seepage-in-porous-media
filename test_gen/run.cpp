/// \reference Wang, M., Wang, J., Pan, N., & Chen, S. (2007). 
/// Mesoscopic predictions of the effective thermal conductivity for microscale random porous media. 
/// Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 75(3), 036702.
/// https://doi.org/10.1103/PhysRevE.75.036702
/// 

#include <time.h>
#include <iostream>
#include <cstring>
#include <string.h>

#include "Porous2D.h"
#include "Porous3D.h"

using namespace std;

void Generate2D();
void Generate3D();
void ReadFromFile();

int main()
{
	clock_t startTime, endTime;
	startTime = clock();

	Generate2D();
	// Generate3D();

	endTime = clock();
	cout << "Totle Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}

void Generate2D()
{
	const int M = 50;
	const int N = 50;
	const double phi = 0.371;			/// target porosity
	const double p_cd = 0.05;		/// core distribution probability
	const double z_h = 0.1;     /// Probability in orthogonal direction
	const double z_v = 0.01;     /// Probability in orthogonal direction
	const double f = 0.025;   /// Probability in oblique direction
	
	int* s = new int[M * N];
	memset(s, 0, sizeof(int)*M*N);		/// initialize to 0

	Porous2D porous(M, N, phi, p_cd, z_h, z_v, f);
	porous.Generation(s);

	//char filename1[100] = "./output/100_100_0.5_0.01_0.1_0.01_0.025";
    //char filetype1[] = ".plt";
    //strncat(filename1, filetype1, sizeof(filename1) - strlen(filename1) - 1);

	char filename2[100] = "./output/50_50_0.3_0.05_0.1_0.01_0.025";
    char filetype2[] = ".txt";
    strncat(filename2, filetype2, sizeof(filename2) - strlen(filename2) - 1);

	//porous.output2tecplot(M, N, s, filename1); 
	porous.output2datatxt(M, N, s, filename2);

	delete[] s;
}

void ReadFromFile()
{
	const int M = 200;
	const int N = 200;
	const double phi = 0.5;			/// target porosity
	const double p_cd = 0.005;		/// core distribution probability
	const double z_h = 0.05;     /// Probability in orthogonal direction
	const double z_v = 0.05;     /// Probability in orthogonal direction
	const double f = 0.0125;   /// Probability in oblique direction
	
	int* s = new int[M * N];
	memset(s, 0, sizeof(int)*M*N);		/// initialize to 0

	Porous2D porous(M, N, phi, p_cd, z_h, z_v, f);

	// If ReadTecplot2D function is used, parameters in Constructor are useless.
	porous.ReadTecplot2D(M, N, s, "Por0.2.plt");

	porous.output2tecplot(M, N, s, "PorMed.plt");
	delete[] s;
}

void Generate3D()
{
	const int Nx = 50;
	const int Ny = 50;
	const int Nz = 50;
	const double phi = 0.5;			/// target porosity
	const double p_cd = 0.005;		/// core distribution probability
	const double p_surface = 0.2;     /// probability to grow on surface
	const double p_edge = 0.0167;	  /// probability to grow on edge, default value (p_surface / 12.0)
	const double p_point = 0.0021;	/// probability to grow on point, default value (p_surface / 96.0)
	
	int* s = new int[Nx * Ny * Nz];
	memset(s, 0, sizeof(int)*Nx*Ny*Nz);		/// initialize to 0

	Porous3D porous(Nx, Ny, Nz, phi, p_cd, p_surface, p_edge, p_point);
	porous.Generation(s);

	porous.output2tecplot("PorMed_3D.plt"); 
	porous.output2datatxt("PorMed_3D_test.txt");

	delete[] s;
}