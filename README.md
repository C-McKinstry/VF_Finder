# VF_Finder
Project done as part of PRACE 'Summer of HPC' 2021

Command line input: ./VF meshfile.msh *numberofintegrationpoints

Supports 1, 6, 16 and 64 integration points. Defaults to 16 if none is inputted. Gaussian quadrature data source https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html

To get this working:

Download the project

Download https://github.com/qnzhou/MshIO as a subdirectory of the project

Install Embree in the project's parent folder

cmake .

make 

./VF sample.msh
