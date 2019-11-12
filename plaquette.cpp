/** A code to calculate the plaquette at a certain point on the lattice and then sum over all points on the lattice. **/

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
/**#include <boost/math/constants/constants.hpp>**/
#include <eigen3/Eigen/Dense>
#include "lattice.cpp"
using namespace Eigen;

/** A problem is, that the size of the lattice has to be known for this definition to work. 
Apparently std::vector can have dynamic size, maybe there is somethings similiar for general arrays.**/

double plaquette(lattice l){

    Matrix<std::complex<double>,2,2> xyplaquettevalue, xzplaquettevalue, yzplaquettevalue;
    double plaquette;

    for(int x = 0; x < 1; x++){
        for (int y = 0; y < 2; y++){
            for (int z = 0; z < 4; z++){
                  xyplaquettevalue = l.links[x][y][z][1].adjoint() * l.links[x][y+1][z][0].adjoint() * l.links[x+1][y][z][1] * l.links[x][y][z][0];
                  xzplaquettevalue = l.links[x][y][z][2].adjoint() * l.links[x][y][z+1][0].adjoint() * l.links[x+1][y][z][2] * l.links[x][y][z][0];
                  yzplaquettevalue = l.links[x][y][z][2].adjoint() * l.links[x][y][z+1][1].adjoint() * l.links[x][y+1][z][2] * l.links[x][y][z][1];

                  plaquette += xyplaquettevalue.trace().real() + xzplaquettevalue.trace().real() + yzplaquettevalue.trace().real();

              }  
           }   
        }
    
}

int main()
{

    /** Initalise the lattice dimenions. **/
    int dimv[3];
    setlattice(dimv);

    lattice l;
    double x = plaquette(l);
    std::cout << x << std::endl;

    return 0;
}


