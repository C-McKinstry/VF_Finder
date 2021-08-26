#include <math.h>
#include <iostream>
#include <fstream>
#include <MshIO/include/mshio/mshio.h>
#define _USE_MATH_DEFINES
using namespace std;

struct VFdata{
  int target;
  double vfactor;
  VFdata(int t, double f){
    target=t; vfactor = f;
  }
};


double dot_prod( double* a, double* b){
  double p = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  return p;
}

double dot_prod( vector<double> a, vector<double> b){
  double p = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  return p;
}

double dot_prod( double* a, vector<double> b){
  double p = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  return p;
}

double len( double* n){
  double l = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  return l;
}

double len( vector<double> n){
  double l = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  return l;
}

double area(vector<double> vb){
  double a;
  double r01[3];
  r01[0] = vb[0] - vb[3]; r01[1] = vb[1] - vb[4]; r01[2] = vb[2] - vb[5];
  double r21[3];
  r21[0] = vb[6] - vb[3]; r21[1] = vb[7] - vb[4]; r21[2] = vb[8] - vb[5];
  a = 0.5*sqrt( pow(len(r01) * len(r21),2) -pow( dot_prod(r01,r21),2));
  return a; 
}

vector<double> calcuate_normal( vector<double> v){
 // n = (P1 - P0) X (P2 - P0)
  
  double n[3];
  n[0] = (v[4] - v[1] )*( v[8] - v[2] ) - ( v[7] - v[1] )*( v[5] - v[2] );
  n[1] = (v[5] - v[2] )*( v[6] - v[0] ) - ( v[3] - v[0] )*( v[8] - v[2] );
  n[2] = (v[3] - v[0] )*( v[7] - v[1] ) - ( v[6] - v[0] )*( v[4] - v[1] );

  double x = len(n);
  double y = 1.0/x;
  double z;
  
  vector<double> result;
  z = y*n[0];
  result.push_back(z);
  z = y*n[1];
  result.push_back(z);
  z = y*n[2];
  result.push_back(z);
  
  return result;
}

vector<double> get_gqweights(string pointnum){
  string filename = "gauss4x4_w.txt";
  if(pointnum == "64"){
    filename = "gauss8x8_w.txt";
  }else if(pointnum == "1"){
    filename = "centroid_w.txt";
  }else if(pointnum == "6"){
    filename = "strang4_w.txt";
  }
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

vector<vector<double>> get_gqpoints(string pointnum){
  string filename = "gauss4x4_x.txt";
  if(pointnum == "64"){
    filename = "gauss8x8_x.txt";printf("64 point sampling\n");
  }else if(pointnum == "1"){
    filename = "centroid_x.txt";printf("1 point sampling\n");
  }else if(pointnum == "6"){
    filename = "strang4_x.txt";printf("6 point sampling\n");
  }else{
    printf("16 point sampling\n");
    }
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

void printdatatofile(int num_tri, vector<vector<double>> centroids, vector<vector<VFdata>> view_factors){
  ofstream OutputFile("viewfactorsdata.txt");
  OutputFile << "#Triangles" << endl;
  OutputFile << num_tri << endl;
  OutputFile << "Triangle Data" << endl;
  for(int i=0;i<num_tri;i++){

    OutputFile << i << " ";
    OutputFile << centroids[i][0] << " " << centroids[i][1] << " " << centroids[i][2];
    int num_seen = view_factors[i].size();
    //printf("Triangle %d sees %d other triangles \n", i, num_seen);
    for(int alpha = 0; alpha < num_seen;alpha++){
      OutputFile << " " << view_factors[i][alpha].target << " " << view_factors[i][alpha].vfactor;
    }
    OutputFile << endl;
  }
  OutputFile.close();
}

vector<double> getnodedata(mshio::MshSpec spec){
  vector<double> nodedata;
  for(int i=0;i<spec.nodes.num_entity_blocks;i++){
      for(int j=0;j<spec.nodes.entity_blocks[i].data.size();j++){
        nodedata.push_back(spec.nodes.entity_blocks[i].data[j]);
      }
    }
  return nodedata;
}

double sigma(double* r, vector<double> normal_a, vector<double> normal_b){
  // Equal to Cos(g1)*Cos(g2) / r*r
  double f =  dot_prod(r, normal_a) * dot_prod(r, normal_b);
  f /= pow( dot_prod(r, r), 2 );
  return f;
}

double calculate_viewfactor(vector<double> va, vector<double> vb, int gqp_num, vector<vector<double>> gqpoints,
 vector<double> gqweights, vector<double> normal_a, vector<double> normal_b, double area_b ){
  double sum_j = 0.0;
  for(int j=0;j<gqp_num;j++){
    //m,n, coefficients describing points chosen via Gaussian Quadrature
    double m = gqpoints[j][0]; double n = gqpoints[j][1];
            
    double wj = gqweights[j];
            
    // Set  r1 = r1(m,n)
    double r1[3];
    r1[0] = (1.0-m-n)*va[0] + m*va[3] + n*va[6];
    r1[1] = (1.0-m-n)*va[1] + m*va[4] + n*va[7];
    r1[2] = (1.0-m-n)*va[2] + m*va[5] + n*va[8];
            
    double sum_i = 0.0;
    for(int i=0;i<gqp_num;i++){
      //s,t, coefficients describing points chosen via Gaussian Quadrature
      double s= gqpoints[i][0];double t=gqpoints[i][1];

      //Set r2 = (1-s-t)*r20 + s*r21 + t*r22
      double r2[3];
      r2[0] = (1.0-s-t)*vb[0] + s*vb[3] + t*vb[6];
      r2[1] = (1.0-s-t)*vb[1] + s*vb[4] + t*vb[7];
      r2[2] = (1.0-s-t)*vb[2] + s*vb[5] + t*vb[8]; 
              
      //Set r_vec = r2 - r1
      double r_vec[3];
      r_vec[0] = r2[0] - r1[0];
      r_vec[1] = r2[1] - r1[1];
      r_vec[2] = r2[2] - r1[2];

      double wi = gqweights[i];
      sum_i += wi*sigma(r_vec, normal_a, normal_b);
      //printf("For m=%f,n=%f,s=%f,t=%f,wj=%f,wi=%f,sum_i=%f\n",m,n,s,t,wj,wi,sum_i);            
    }
    sum_j += wj*sum_i;
  }
  return  - area_b*sum_j / M_PI;
}

int sumonetonminusone(int n){
  int count = 0;
  for(int i=1;i<n;i++){count += i;}
  return count;
}