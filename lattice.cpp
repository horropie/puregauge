#include "lattice.h"
#include <eigen3/Eigen/Dense>
#include <random>
#include <vector>
#include <complex>
/**#include <boost/math/constants/constants.hpp>**/

using namespace std;

std::default_random_engine generator;
std::uniform_real_distribution<double> uni_real_distr(0.0,1.0);

lattice::lattice()
{
xdim= sizeof(links)/sizeof(links[0]);
ydim= sizeof(links[0])/sizeof(links[0][0]);
zdim= sizeof(links[0][0]);

     for(int i =0;i<xdim;i++){
             for(int j =0;j<ydim;j++){
                    for(int k=0;k<zdim;k++){

                        links[i][j][k][0]=rand();
                        links[i][j][k][1]=rand();
                        links[i][j][k][2]=rand();
                        }


             }

     }

}

lattice::~lattice()
{
    //dtor
}
Eigen::Matrix<std::complex<double>,2,2> lattice::rand()
{
    Eigen::Matrix<std::complex<double>,2,2> randm;
    double rand[3];
	double norm(2.0);
	while(norm>1.0){
     rand[0]=uni_real_distr(generator);
     rand[1]=uni_real_distr(generator);
     rand[2]=uni_real_distr(generator);
	 norm=sqrt(rand[0]*rand[0]+rand[1]*rand[1]+rand[2]*rand[2]);
	} // random direction on 3 sphere
	rand[0]/=norm;
	rand[1]/=norm;
	rand[2]/=norm;
	const double phase(2.0*3.1415*uni_real_distr(generator));
	const double cosn(cos(phase));
	const double sinn(sin(phase));
	randm(0,0)=std::complex<double>(cosn,0.0)+std::complex<double>(0.0,sinn*rand[2]);
	randm(0,1)=std::complex<double>(0.0,sinn*rand[0])+std::complex<double>(sinn*rand[1],0.0);
	randm(1,1)=conj(randm(0,0));
	randm(1,0)=-conj(randm(0,1));
	return randm;
}


