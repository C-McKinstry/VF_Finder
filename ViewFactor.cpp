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

  vector<vector<double>> view_factor_array( num_tri , vector<double> (num_tri, 0));

  RTCRayHit rayhit;
  float r_vec[3];
  double f;
  double alpha_centre[3];
  double beta_centre[3];
  double f_sum;
  
  for(int a =0; a<num_tri;a++){
    for(int b=a; b<num_tri;b++){
      if(a==b || dot_prod(normal_array[a], normal_array[b]) >= 0){
        //Situation where no ray tracing is necessary
        view_factor_array[a][b]=0.0;
        view_factor_array[b][a]=0.0;
      }
      else{
        

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

        if(rayhit.hit.geomID != b){
          //Blocked
          view_factor_array[a][b]=0.0;
          view_factor_array[b][a]=0.0;
        }
        else{
          //Hit
          f_sum=0;
          for(int alpha=0;alpha<4;alpha++){
            for(int beta=0;beta<4;beta++){
              get_centre(centroid_array[a], v[a], alpha, alpha_centre);
              get_centre(centroid_array[b], v[b], beta, beta_centre);
              r_vec[0] = beta_centre[0]-alpha_centre[0];
              r_vec[1] = beta_centre[1]-alpha_centre[1];
              r_vec[2] = beta_centre[2]-alpha_centre[2];  
              f =  - dot_prod(r_vec, normal_array[a]) * dot_prod(r_vec, normal_array[b]) * area_array[b];
              f /= pow( dot_prod(r_vec, r_vec), 2 ) * M_PI * 16;
              f_sum += f;
            }
          }
          view_factor_array[a][b] = f_sum;
          view_factor_array[b][a] = ( area_array[a] * f_sum ) / area_array[b];
          printf( " %d sees %d, VF= %F \n", a, b, f);
        }
        
      }
    }
  } 
  rtcReleaseScene(scene);
  rtcReleaseDevice(device);

  double row_sum;
  for(int i=0;i<num_tri;i++){
    row_sum=0.0;
    for(int j=0;j<num_tri;j++){
      row_sum += view_factor_array[i][j];
    }
    printf("Row %d adds up to %f \n", i, row_sum);
  }


  return 0;
}