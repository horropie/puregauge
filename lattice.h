#ifndef LATTICE_H
#define LATTICE_H

#include <Eigen/Dense>

using namespace Eigen;

class lattice
{
    public:
        lattice();
        int xdim;
        int ydim;
        int zdim;
        double T;
        virtual ~lattice();
        void update(int xpos, int ypos, int zpos, int direction);
        double energy(void);
        double plaquette(int x, int y, int z, int dir);
        Matrix<std::complex<double>,2,2> links[2][3][5][3] ;

    protected:

    private:
        Matrix<std::complex<double>,2,2> rand(void);
};

#endif // LATTICE_H
