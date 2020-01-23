#ifndef LATTICE_H
#define LATTICE_H

#include <Eigen/Dense>

using namespace Eigen;
class lattice
{
    public:
        lattice(int space,int t);
         int xdim;
         int ydim;
         int zdim;
         int tdim;
        double T;
        virtual ~lattice();
        void update(int xpos, int ypos, int zpos, int direction, int time, Eigen::Matrix<std::complex<double>,2,2> newMatrix);
        double energy(void);
        double plaquette(int x, int y, int z, int dir, int time);
        Matrix<std::complex<double>,2,2> *****links;


    protected:

    private:
        Eigen::Matrix<std::complex<double>,2,2> rand( void);

};

#endif // LATTICE_H
