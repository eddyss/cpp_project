#include <iostream>
#include <sys/time.h>
#include "var_alea.hpp"
#include "storage.h"

using namespace std;

int main() {
	init_alea();
	struct timeval tbegin ,tend;
	double texec=0;
	gettimeofday (&tbegin,NULL);

//	uniform U(0,1);
	
	gaussian G(0,1);
//	std::cout << G() << std::endl;
	
//	t_student T(3);
//	std::cout << T() << std::endl;
//	std::cout << "densité" " " << T.pdf(1) << std::endl;

/*	chi_deux X(1);
	std::cout << X() << std::endl;
	std::cout << "densité" " " << X.pdf(1) << std::endl;
	std::cout << "répartition" " " << X.cdf(1.5) << std::endl;*/
	
/*	inverse_gaussian Y(0.5,1);
	std::cout << Y() << std::endl; 
	std::cout << "densité" " " << Y.pdf(1) << std::endl;
	std::cout << "répartition" " " << Y.cdf(1) << std::endl;

	normal_inverse_gaussian Z(0.8,0.1,0,2);
//	std::cout << "densité" " " << Z.pdf(1) << std::endl;
//	std::cout << "répartition" " " << Z.SimpleMonteCarlo(1, 500) << std::endl;*/
//	std::cout << Z() << std::endl;    
//	std::cout << Z.mean() << std::endl;
//	std::cout << Z.var() << std::endl;

//	std::cout << inverseTableau (0.001,A,B) << std::endl;
	
/*	int const n=200;
	double A[n+1];
	double B[n+1];
	remplirTableau (T, A, B, n, 500000);
	std::cout << inverseTableau (0.0001,A,B,n) << std::endl;
	std::cout << inverseTableau (0.1,A,B,n) << std::endl;
	std::cout << inverseTableau (0.2,A,B,n) << std::endl;
	std::cout << inverseTableau (0.5,A,B,n) << std::endl;
	std::cout << inverseTableau (0.8,A,B,n) << std::endl;
	std::cout << inverseTableau (0.9,A,B,n) << std::endl;
	std::cout << inverseTableau (0.9999,A,B,n) << std::endl;
	remplirTableauAnthetique (T, A, B, n, 500000);
	std::cout << inverseTableau (0.0001,A,B,n) << std::endl;
	std::cout << inverseTableau (0.1,A,B,n) << std::endl;
	std::cout << inverseTableau (0.2,A,B,n) << std::endl;
	std::cout << inverseTableau (0.5,A,B,n) << std::endl;
	std::cout << inverseTableau (0.8,A,B,n) << std::endl;
	std::cout << inverseTableau (0.9,A,B,n) << std::endl;
	std::cout << inverseTableau (0.9999,A,B,n) << std::endl;*/

	int const n=100;
	char file[] = "test.txt";
	char file_ecriture[] = "test_ecriture.txt";
	double X[n+1];
	double Y[n+1];
	remplirTableauAnthetique (G, X, Y, n, 50000);
//	for(int i=0;i<n;i++){
//	cout << X[i] << " / " << Y[i] << endl;
//	}

	store_array(file_ecriture,n,&X[0],&Y[0]);

	double A[n+1];
	double B[n+1];

	get_array(file_ecriture,n,&A[0],&B[0]);

	for (int i=0; i<=n;i++){std::cout << A[i] << " " << B[i] << std::endl;}

	gettimeofday(&tend, NULL);
	texec= ((double)(1000*(tend.tv_sec-tbegin.tv_sec)+((tend.tv_usec-tbegin.tv_usec)/1000)))/1000.;
	std::cout << texec << std::endl;
	return 0;
};
