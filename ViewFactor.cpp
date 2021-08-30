// Copyright 2009-2021 Intel Corporation
// SPDX-License-Identifier: Apache-2.0

#include <embree3/rtcore.h>
#include <limits>
#include <iostream>
#include <MshIO/include/mshio/mshio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <functions.h>
#include <mpi.h>
using namespace std;

int main(int argc, char *argv[])
{
  if (argc < 2){std::cout <<"You haven't entered a file\n";return 0;}
  string pointnum = "16";
  if(argc == 3){pointnum = argv[2];}
  if(argc > 3){std::cout<< "Too many arguments \n";return 0;}

  //Import .msh file to be read
  mshio::MshSpec spec = mshio::load_msh(argv[1]);

  vector<double> nodedata = getnodedata(spec);

  
  for(int q = 0;q<spec.elements.num_entity_blocks;q++){
    printf("Element block %d is of type %d \n", q, spec.elements.entity_blocks[q].element_type);
    if(spec.elements.entity_blocks[q].element_type == 2){
      printf("Block %d has %d triangles \n",q,spec.elements.entity_blocks[q].num_elements_in_block);
    }
  }
  for(int p = 0;p<spec.nodes.num_entity_blocks;p++){
    printf("Node block %d has %d nodes, belong to %dd entity %d \n", p,
    spec.nodes.entity_blocks[p].num_nodes_in_block, spec.nodes.entity_blocks[p].entity_dim, spec.nodes.entity_blocks[p].entity_tag);
  }
  
  RTCDevice device = rtcNewDevice(NULL);
  RTCScene scene   = rtcNewScene(device);
  RTCGeometry geom;
  
  float* vb;
  unsigned* ib;
  int num_tri = 0;
  vector<vector<double>> centroid_array;
  vector<vector<double>> normal_array;
  vector<double> area_array;
  vector<vector<double>> v;
  int n1,n2,n3;
  vector<double> emptyrowd;
  
  for(int eblock_num = 0;eblock_num<spec.elements.num_entity_blocks;eblock_num++){
    if(spec.elements.entity_blocks[eblock_num].element_type == 2){
        for(int i=0;i<spec.elements.entity_blocks[eblock_num].num_elements_in_block; i++){
          auto& triangles = spec.elements.entity_blocks[eblock_num];

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
          vb[0] = nodedata[3*(n1-1)];vb[1] = nodedata[3*(n1-1) +1];vb[2] = nodedata[3*(n1-1)+2];
          vb[3] = nodedata[3*(n2-1)];vb[4] = nodedata[3*(n2-1) +1];vb[5] = nodedata[3*(n2-1)+2];
          vb[6] = nodedata[3*(n3-1)];vb[7] = nodedata[3*(n3-1) +1];vb[8] = nodedata[3*(n3-1)+2];

          v.push_back(emptyrowd);
          for(int j=0;j<9;j++){ //Float -> Double
            v[num_tri].push_back(vb[j]);
          }

          //Calculates centre of triangle for later use
          centroid_array.push_back(emptyrowd);
          centroid_array[num_tri].push_back( ( v[num_tri][0] + v[num_tri][3] + v[num_tri][6] ) / 3 );
          centroid_array[num_tri].push_back( ( v[num_tri][1] + v[num_tri][4] + v[num_tri][7] ) / 3 );
          centroid_array[num_tri].push_back( ( v[num_tri][2] + v[num_tri][5] + v[num_tri][8] ) / 3 );
          normal_array.push_back(calcuate_normal(v[num_tri]));
          area_array.push_back( area(v[num_tri]) );

          rtcCommitGeometry(geom);
          rtcAttachGeometry(scene, geom);
          rtcReleaseGeometry(geom);

          
          //printf("Triangle %d has points: \n %f %f %f \n %f %f %f \n %f %f %f \n", num_tri,v[num_tri][0], v[num_tri][1],
          //v[num_tri][2],v[num_tri][3],v[num_tri][4],v[num_tri][5],v[num_tri][6],v[num_tri][7],v[num_tri][8]);
          //printf("It has area %f, normal: %f %f %f \n", area_array[num_tri], normal_array[num_tri][0],normal_array[num_tri][1],normal_array[num_tri][2]);
          
          //printf("Triangle %d has nodes %d %d %d, and normal:  %f %f %f \n", num_tri, n1,n2,n3,
          //normal_array[num_tri][0],normal_array[num_tri][1],normal_array[num_tri][2]);
          num_tri++;
        }
    }
  }


  rtcCommitScene(scene);

  vector<vector<VFdata>> view_factors;
  vector<VFdata> empty_row;
  RTCRayHit rayhit;

  
  //Fethces the points and weights used for Gaussian quadrature integration
  vector<vector<double>> gqpoints = get_gqpoints(pointnum);
  vector<double> gqweights = get_gqweights(pointnum);
  int gqp_num = gqweights.size();
  
  //Creates the array the view factors will go into
  for(int k=0;k<num_tri;k++){
    view_factors.push_back(empty_row);
  }

  int numpairs = sumonetonminusone(num_tri);
  int myrank; int numranks;
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numranks);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

  int pairsperrank = numpairs / numranks ;
  int start = pairsperrank*myrank;
  int end = (myrank+1)*pairsperrank;
  if(myrank = numranks - 1){
    end = numpairs;
  }

  int currentpairnum = -1;
  
  for(int a =0; a<num_tri;a++){//Each emmiting traingle
    for(int b=a; b<num_tri;b++){//Each possible recieving 
    currentpairnum++;
    if(start <= currentpairnum && currentpairnum < end){
      vector<double> r_ctc;
      r_ctc.push_back( centroid_array[b][0] - centroid_array[a][0]);
      r_ctc.push_back( centroid_array[b][1] - centroid_array[a][1]);
      r_ctc.push_back( centroid_array[b][2] - centroid_array[a][2]);

      if(a !=b && dot_prod(r_ctc,normal_array[a]) >= 0 && dot_prod(r_ctc,normal_array[b]) <=0 ){
        //Checks where VF should be 0

        //Where the ray comes from
        rayhit.ray.org_x = centroid_array[a][0];
        rayhit.ray.org_y = centroid_array[a][1];
        rayhit.ray.org_z = centroid_array[a][2];
        
        //Direction the ray goes in
        rayhit.ray.dir_x = centroid_array[b][0] - centroid_array[a][0];
        rayhit.ray.dir_y = centroid_array[b][1] - centroid_array[a][1];
        rayhit.ray.dir_z = centroid_array[b][2] - centroid_array[a][2];
        
        //Ray starts at a, carries on until it hits something
        rayhit.ray.tnear  = 0.f;
        rayhit.ray.tfar   = std::numeric_limits<float>::infinity();
        rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
        RTCIntersectContext context;
        rtcInitIntersectContext(&context);
        rtcIntersect1(scene, &context, &rayhit);

        if(rayhit.hit.geomID == b){
          //If it hits something, and it's b

          //printf(" Triangle %d  sees triangle %d ", a, b);
          
          double f = calculate_viewfactor(v[a], v[b], gqp_num, gqpoints, gqweights,
           normal_array[a], normal_array[b], area_array[b]);
          //double f =  area_array[b]*sum_j;
          //double ar = area_array[a]*area_array[b]*sum_j;

          view_factors[a].push_back(VFdata(b, f));
          view_factors[b].push_back(VFdata(a, ( area_array[a] * f ) / area_array[b] )) ;
          printf("Triangle %d sees triangle %d, VF= %f \n", a,b,f);
          //printf("Area = %f\n", ar);

        }
        
      }
      
    }
    
    }
  } 

  MPI_Finalize();
  rtcReleaseScene(scene);
  rtcReleaseDevice(device);
  

  /*
  double total = 0.0;
  double area_sum;
  for(int t=0;t<4;t++){
    //printf("Triangle %d sees %d triangles,", t, view_factors[t].size());
    for(int s=0;s<view_factors[t].size();s++){
      total += (area_array[t])* view_factors[t][s].vfactor;
    }
    area_sum += area_array[t];
  }
  total /=area_sum;

  printf(" total view factor %f \n", total);
  */
  printdatatofile(num_tri, centroid_array,  view_factors);


  return 0;
}