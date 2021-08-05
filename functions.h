#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;



double dot_prod( float* a, float* b){
  double p = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  return p;
}

double dot_prod( double* a, double* b){
  double p = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  return p;
}

double dot_prod( float* a, double* b){
  double p = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  return p;
}

double len( float* n){
  double l = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  return l;
}

double len( double* n){
  double l = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  return l;
}

double area(double* vb){
  double a;
  double r01[3];
  r01[0] = vb[0] - vb[3]; r01[1] = vb[1] - vb[4]; r01[2] = vb[2] - vb[5];
  double r21[3];
  r21[0] = vb[6] - vb[3]; r21[1] = vb[7] - vb[4]; r21[2] = vb[8] - vb[5];
  a = 0.5*sqrt( pow(len(r01) * len(r21),2) -pow( dot_prod(r01,r21),2));
  return a; 
}

int save_normal( double* v, double* destination){
 // n = (P1 - P0) X (P2 - P0)
  
  double n[3];
  n[0] = (v[4] - v[1] )*( v[8] - v[2] ) - ( v[7] - v[1] )*( v[5] - v[2] );
  n[1] = (v[5] - v[2] )*( v[6] - v[0] ) - ( v[3] - v[0] )*( v[8] - v[2] );
  n[2] = (v[3] - v[0] )*( v[7] - v[1] ) - ( v[6] - v[0] )*( v[4] - v[1] );

  double x = len(n);
  double y = 1.0/x;
  double z;
  
  z = y*n[0];
  destination[0] = z;
  z = y*n[1];
  destination[1] = z;
  z = y*n[2];
  destination[2] = z;
  
  return 0;
}

void get_centre( double* centroid, double* v, int n, double* destination ){
  if( n == 0){
    destination[0] = centroid[0];destination[1] = centroid[1];destination[2] = centroid[2];
  }
  else{
    destination[0] = 0.5*( centroid[0] + v[ 3*(n-1)]);
    destination[1] = 0.5*( centroid[1] + v[ 3*(n-1) + 1]);
    destination[2] = 0.5*( centroid[2] + v[ 3*(n-1) + 2]);
  }
}

vector<double> get_gqweights(string filename = "gauss8x8_w.txt"){
  vector<double> w;
  ifstream win_file;
  win_file.open(filename);
  if (!win_file.is_open()) {
        std::cout << "Cannot read quadrature weights input file \n";
  }
  double input_double;
  while (win_file >> input_double) {
        w.push_back(input_double);
    }
  win_file.close();
  return w;
}

vector<vector<double>> get_gqpoints(string filename = "gauss8x8_x.txt"){
  vector<vector<double>> points;
  double input_double;
  ifstream xin_file;
  xin_file.open(filename);
  if (!xin_file.is_open()) {
        std::cout << "Cannot read quadrature points input file \n";
  }
  int i=0;//Counts whether the next entry is an x or y cooridnate
  int j=0;//Counts the row we're on
  vector<double> empty_row;
  while (xin_file >> input_double) {
    if( i ==0){
      points.push_back(empty_row);
      points[j].push_back(input_double);
      i=1;
    }
    else{
      points[j].push_back(input_double);
      i=0;
      j++;
    }
  }
  xin_file.close();
  return points;
}