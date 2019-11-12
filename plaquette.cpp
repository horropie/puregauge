/** A code to calculate the plaquette at a certain point on the lattice. It should be compatible with then integrating over the whole lattice. 
The input is a random matrix U and the lattice array. The output is all the plaquette values for points on the lattice.**/

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <boost/math/constants/constants.hpp>
#include "Eigen/Eigen/Dense"

/** A problem is, that the size of the lattice has to be known for this definition to work. 
Apparently std::vector can have dynamic size, maybe there is somethings similiar for general arrays.**/
double plaquette(double &(links)[][][][3], double &(matrix)[2][2], int (&dimv)[3]){

    double xyplaquettevalue, xzplaquettevalue, yzplaquettevalue;

    for(int x = 0; x < dimv[0]; x++){
        for (int y = 0; y < dimv[1]; y++){
            for (int z = 0; z < dimv[2]; z++){
                  xyplaquettevalue = links[x][y][z][1].adj() * links[x][y+1][z][0].adj() * links[x+1][y][z][1] * links[x][y][z][0];
                  xzplaquettevalue = links[x][y][z][1].adj() * links[x][y][z+1][0].adj() * links[x+1][y][z][1] * links[x][y][z][0];
                  yzplaquettevalue = links[x][y][z][1].adj() * links[x][y][z+1][0].adj() * links[x][y+1][z][1] * links[x][y][z][0];

                  /**What output do we want? Single doubles or whole arrays? Or saving?**/
              }  
           }   
        }
    
}

