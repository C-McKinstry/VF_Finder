// Copyright 2009-2021 Intel Corporation
// SPDX-License-Identifier: Apache-2.0

#include <embree3/rtcore.h>
#include <limits>
#include <iostream>
#include <MshIO/include/mshio/mshio.h>
#include <typeinfo>
#define _USE_MATH_DEFINES
#include <math.h>
#include <functions.h>
using namespace std;

double factor(double* r, double* normal_a, double* normal_b){
  // Equal to Cos(g1)*Cos(g2) / r*r
  double f =  - dot_prod(r, normal_a) * dot_prod(r, normal_b);
  f /= pow( dot_prod(r, r), 2 );
  return f;
}

struct VFdata{
  int target;
  double vfactor;
  VFdata(int t, double f){
    target=t; vfactor = f;
  }
};

int main(int argc, char *argv[])
{
  if (argc < 2){std::cout <<"You haven't entered a file\n";return 0;}
  if (argc > 2){std::cout <<"One file at a time\n";return 0;}

  //Import .msh file to be read
  mshio::MshSpec spec = mshio::load_msh(argv[1]);

  if (spec.nodes.num_entity_blocks != 1){std::cout<< "Expected only one node block\n";return 0;}
  
  //Find the entity block corresponding to triangles
  int i = -1; int type;
  while(type != 2){ 
    i++;
    type = spec.elements.entity_blocks[i].element_type; 
    }


  auto& triangles = spec.elements.entity_blocks[i];
  auto& nodes = spec.nodes.entity_blocks[0];
  int num_tri = triangles.num_elements_in_block;
  int num_nodes = spec.nodes.num_nodes;

  RTCDevice device = rtcNewDevice(NULL);
  RTCScene scene   = rtcNewScene(device);
  RTCGeometry geom;

  float* vb;
  unsigned* ib;

  double centroid_array[num_tri][3];
  double normal_array[num_tri][3];
  double area_array[num_tri];
  double v[num_tri][9];
  int n1,n2,n3;
  for(int i = 0;i<num_tri;i++){

    //Create new holders for geometry data
    geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
    vb = (float*) rtcSetNewGeometryBuffer(geom,
      RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3*sizeof(float), 3);
    ib = (unsigned*) rtcSetNewGeometryBuffer(geom,
      RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3*sizeof(unsigned), 1);
    ib[0] = 0; ib[1] = 1; ib[2] = 2;

    //Node # of the triangle's vertices
    n1 = triangles.data[4*i +1];
    n2 = triangles.data[4*i +2];
    n3 = triangles.data[4*i +3];
    
    //Coordinates of the triangle
    vb[0] = nodes.data[3*(n1-1)];vb[1] = nodes.data[3*(n1-1) +1];vb[2] = nodes.data[3*(n1-1)+2];
    vb[3] = nodes.data[3*(n2-1)];vb[4] = nodes.data[3*(n2-1) +1];vb[5] = nodes.data[3*(n2-1)+2];
    vb[6] = nodes.data[3*(n3-1)];vb[7] = nodes.data[3*(n3-1) +1];vb[8] = nodes.data[3*(n3-1)+2];

    for(int j=0;j<9;j++){ //Float -> Double
      v[i][j] = vb[j];
    }

    //Calculates centre of triangle for later use
    centroid_array[i][0] = ( v[i][0] + v[i][3] + v[i][6] ) / 3 ;
    centroid_array[i][1] = ( v[i][1] + v[i][4] + v[i][7] ) / 3 ;
    centroid_array[i][2] = ( v[i][2] + v[i][5] + v[i][8] ) / 3 ;
    save_normal( v[i], normal_array[i]);
    area_array[i] = area(v[i]);

    rtcCommitGeometry(geom);
    rtcAttachGeometry(scene, geom);
    rtcReleaseGeometry(geom);
    
  }
  rtcCommitScene(scene);

  vector<vector<VFdata>> view_factors;
  vector<VFdata> empty_row;
  RTCRayHit rayhit;
  double r_vec[3];
  double r1[3]; double r2[3];
  double s; double t; double m; double n; double wi; double wj;
  double f; double x;
  double alpha_centre[3];
  double beta_centre[3];
  double f_inner;

  vector<vector<double>> gqpoints = get_gqpoints();
  vector<double> gqweights = get_gqweights();
  int gqp_num = gqweights.size();
    
  for(int k=0;k<num_tri;k++){
    view_factors.push_back(empty_row);
  }
  
  for(int a =0; a<num_tri;a++){
    for(int b=a; b<num_tri;b++){
      if(a !=b && dot_prod(normal_array[a], normal_array[b]) < 0){
        rayhit.ray.org_x = centroid_array[a][0];
        rayhit.ray.org_y = centroid_array[a][1];
        rayhit.ray.org_z = centroid_array[a][2];
        r_vec[0] = centroid_array[b][0] - centroid_array[a][0]; rayhit.ray.dir_x=r_vec[0];
        r_vec[1] = centroid_array[b][1] - centroid_array[a][1]; rayhit.ray.dir_y=r_vec[1];
        r_vec[2] = centroid_array[b][2] - centroid_array[a][2]; rayhit.ray.dir_z=r_vec[2];

        rayhit.ray.tnear  = 0.f;
        rayhit.ray.tfar   = std::numeric_limits<float>::infinity();
        rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
        RTCIntersectContext context;
        rtcInitIntersectContext(&context);
        rtcIntersect1(scene, &context, &rayhit);

        if(rayhit.hit.geomID == b){
          printf(" Triangle %d  sees triangle %d ", a, b);
          f = 0.0;
          for(int j=0;j<gqp_num;j++){
            x = 0.0;
            m = gqpoints[j][0];n=gqpoints[j][1];
            // Set  r1 = r1(m,n)
            r1[0] = (1-m-n)*v[a][0] + m*v[a][3] + n*v[a][6];
            r1[1] = (1-m-n)*v[a][1] + m*v[a][4] + n*v[a][7];
            r1[2] = (1-m-n)*v[a][2] + m*v[a][5] + n*v[a][8];

            for(int i=0;i<gqp_num;i++){
              //s,t, coefficients describing points chosen via Gaussian Quadrature
              s= gqpoints[i][0];t=gqpoints[i][1];

              //Set r2 = (1-s-t)*r20 + s*r21 + t*r22
              r2[0] = (1-s-t)*v[b][0] + s*v[b][3] + t*v[b][6];
              r2[1] = (1-s-t)*v[b][1] + s*v[b][4] + t*v[b][7];
              r2[2] = (1-s-t)*v[b][2] + s*v[b][5] + t*v[b][8]; 
              
              //Set r_vec = r2 - r1
              r_vec[0] = r2[0] - r1[0];
              r_vec[1] = r2[1] - r1[1];
              r_vec[2] = r2[2] - r1[2];

              wi = gqweights[i];
              f_inner = wi*factor(r_vec, normal_array[a], normal_array[b]);              
              x += f_inner;
            }
            x *= (area_array[b]/2);
            wj = gqweights[j];
            f += wj*x;
          }
          f /= - 2*M_PI;
          printf(", View Factor = %f \n", f);

          view_factors[a].push_back(VFdata(b, f));
          view_factors[b].push_back(VFdata(a, ( area_array[a] * f ) / area_array[b] )) ;
          //printf( " %d sees %d, VF= %F \n", a, b, f);
        }
        
      }
      
    }
  } 
  rtcReleaseScene(scene);
  rtcReleaseDevice(device);

printf("Code ran \n");

  double f_add;
  int zero_rows = 0;
  for(int i=0;i<num_tri;i++){
    f_add = 0.0;
    if(view_factors[i].size() > 0){
      for(int j=0;j<view_factors[i].size();j++){
        // printf("F %d %d = %f \n", i, view_factors[i][j].target, view_factors[i][j].vfactor);
        f_add += view_factors[i][j].vfactor;
      }
    }else{ zero_rows++;}
    printf( " Triangle %d radiates a total of %f \n", i , f_add );
  }
  printf( "There are %d triangles that see nothing \n", zero_rows);

  return 0;
}