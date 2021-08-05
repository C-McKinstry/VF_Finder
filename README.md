# VF_Finder
Project done as part of PRACE 'Summer of HPC' 2021

Current gaussian quadrature data sets sample 16 and 64 points respectively. Change which one in the get_gq... functions in functions.h

To get this working:

Download the project

Download https://github.com/qnzhou/MshIO as a subdirectory of the project

Install Embree in the project's parent folder

cmake .

make 

./VF sample.msh
