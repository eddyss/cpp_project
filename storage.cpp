#include <iostream>
#include <fstream>
#include <cstdlib>
#include "storage.h"

using namespace std;

const string delim = ",";

/*int main(){
  char file[] = "data/test.txt";
  char file_ecriture[] = "data/test_ecriture.txt";
  double X[3];
  double Y[3];
  get_array(file,3,&X[0],&Y[0]);
  
  for(int i=0;i<3;i++){
    cout << X[i] << " / " << Y[i] << endl;
  }

  store_array(file_ecriture,3,&X[0],&Y[0]);

  return 0;
}*/

int store_array(char * file, int N, double *X, double *Y){
  ofstream fichier(file, ios::out | ios::trunc);
 
  if(fichier){
    for(int i=0;i<N;i++){
      fichier << X[i] << "," << Y[i] << endl;
    }
 
    fichier.close();
  }
  else
    return 0;
 
  return 1;
}

int get_array(char * file, int N, double *X, double *Y){
  ifstream fichier(file);  //ouverture
  if(fichier)
    {
      string contenu;
      string x, y;
      int i=0;
      
      while(i<N and getline(fichier, y)){ //lecture de la ligne
	x = y.substr(0, y.find(delim));
	y = y.substr(y.find(delim)+1,sizeof(x)/sizeof(char)-x.find(delim)-1);

	X[i] = atof(x.c_str());
	Y[i] = atof(y.c_str());
 
	i++;
      }
      fichier.close();
    }
  else{
    return 0;
  }
 
  return 1;
}
