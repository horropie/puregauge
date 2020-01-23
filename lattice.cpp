#include "lattice.h"
#include <Eigen/Dense>
#include <random>
#include <vector>
#include <complex>
//#include <boost/math/constants/constants.hpp>

using namespace std;

std::default_random_engine generator;
std::uniform_real_distribution<double> uni_real_distr(0.0,1.0);

lattice::lattice(int space,int t)
{
xdim=space;
ydim=space;
zdim=space;
tdim=t;


links = new Eigen::Matrix<std::complex<double>,2,2>****[xdim];

        for(int i =0;i<xdim;i++){
                links[i]= new Eigen::Matrix<std::complex<double>,2,2>***[ydim];
                for(int j =0;j<ydim;j++){
                    links[i][j]= new Eigen::Matrix<std::complex<double>,2,2>**[zdim];
                    for(int k=0;k<zdim;k++){
                        links[i][j][k]=new Eigen::Matrix<std::complex<double>,2,2>*[4];
                        for(int r=0;r<4;r++){
                            links[i][j][k][r]=new Eigen::Matrix<std::complex<double>,2,2>[tdim];
                            for(int t=0;t<tdim;t++){
                                links[i][j][k][r][t]=rand();
              }
              }
              }
              }
              }



}

lattice::~lattice()
{

}
void lattice::update(int xpos, int ypos, int zpos, int direction, int time, Eigen::Matrix<std::complex<double>,2,2> newMatrix)
{
   links[xpos][ypos][zpos][direction][time]=newMatrix;
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
//cout<<"rand: "<<rand<<endl;
	const double phase(2.0*3.1415*uni_real_distr(generator));
	//cout<<"phase: "<<phase<<endl;
	const double cosn(cos(phase));
	const double sinn(sin(phase));
	//cout<<"cosn: "<<cosn<<endl;
	//cout<<"sinn: "<< sinn<<endl;
	randm(0,0)=std::complex<double>(cosn,0.0)+std::complex<double>(0.0,sinn*rand[2]);
	randm(0,1)=std::complex<double>(0.0,sinn*rand[0])+std::complex<double>(sinn*rand[1],0.0);
	randm(1,1)=conj(randm(0,0));
	randm(1,0)=-conj(randm(0,1));
	return randm;


}
