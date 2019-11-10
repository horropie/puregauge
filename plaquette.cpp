/** A code to calculate the plaquette at a certain point on the lattice. It should be compatible with then integrating over the whole lattice. 
The input is a random matrix U and the lattice array. The output is the plaquette for a certain point on the lattice.**/

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <boost/math/constants/constants.hpp>
#include "Eigen/Eigen/Dense"

/** A problem is, that the size of the lattice has to be known for this definition to work.
Can we really apply our Plaquette formula to the xyz case (It seems to be just a xy case)?  
Also, how do we safe the matrices that correspond to the links? Do we construct them on the fly?**/
void plaquette(double &(lattice)[][][][3], double &(matrix)[2][2],){
    /**Pseudocode: U(links[x][y][z][1]).adj() * U(links[x][y+1][z][0]).adj() * U(links[x+1][y][z][1]) * U(links[x][y][z][0])**/  
}

