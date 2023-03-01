# DGG-VTP
A direct and computationally efficient method for computing discrete geodesics
To know how to compile the parallel version for DGG-VTP, please look at : https://github.com/djie-0329/PVTP to know how to compile.
Meanwhile, we provide the build solution for Single Threaded DGG-VTP, so we can use it easily. We can use Visual Studio 2013 to open the build solution. For more recent versions of Visual Studio, we need to ensure VS2013 toolset (v140) is properly installed. 

To run the single thread version, we can use :
AVTP_DGG_final.exe -m [mesh name] -s [source number] -eps [error settings] -o [txt output file name]
Before running, make sure we have the obj mesh file in AVTP_DGG_final folder first.

Meanwhile, to run the parallel version, we can use :
DGG-VTP.exe -alg 4 -m [mesh obj name] -s [source index] -eps [error setting] -tol [tolerance to measure correctness] -np [number of treads] -o 
