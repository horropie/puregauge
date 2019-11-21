#include <iostream>
#include "lattice.h"
#include <Eigen/Dense>
#include <random>
#include <vector>
#include <complex>
#include <boost/math/constants/constants.hpp>
#include <fstream>
using namespace std;
//definiere Parameter (Anzahl Wdh., Temperatur, Zufallsgenerator)
static int rep=10;
double beta= 0.7;
double N =1.0;

Eigen::Matrix<std::complex<double>,2,2> test;
Eigen::Matrix<std::complex<double>,2,2> temp;

fstream f;

double alt=0.0;
double neu=0.0;
double p;



std::default_random_engine Generator;
std::uniform_real_distribution<double> uni_real_dist(0.0,1.0);

// Methode zum Erzeugen von Zufallsmatrix
Eigen::Matrix<std::complex<double>,2,2> zuf()
{
    Eigen::Matrix<std::complex<double>,2,2> randm;
    double rand[3];
	double norm(2.0);
	while(norm>1.0){
     rand[0]=uni_real_dist(Generator);
     rand[1]=uni_real_dist(Generator);
     rand[2]=uni_real_dist(Generator);
	 norm=sqrt(rand[0]*rand[0]+rand[1]*rand[1]+rand[2]*rand[2]);
	} // random direction on 3 sphere
	rand[0]/=norm;
	rand[1]/=norm;
	rand[2]/=norm;
	const double phase(2.0*3.1415*uni_real_dist(Generator));
	const double cosn(cos(phase));
	const double sinn(sin(phase));
	randm(0,0)=std::complex<double>(cosn,0.0)+std::complex<double>(0.0,sinn*rand[2]);
	randm(0,1)=std::complex<double>(0.0,sinn*rand[0])+std::complex<double>(sinn*rand[1],0.0);
	randm(1,1)=conj(randm(0,0));
	randm(1,0)=-conj(randm(0,1));
	return randm;
}




int main()
{
    lattice l ;
     string datei= "Metropolis.txt";
        f.open(datei,ios::out);
       f << "O" << endl;

    // Berechne die Summe der Plaquettes des Gitters, schpeichere in plaquette

    double plaquette=0.0;
    Eigen::Matrix<std::complex<double>,2,2> xyplaquettevalue;

     for(int x = 0; x < l.xdim; x++){
        for (int y = 0; y < l.ydim; y++){
            for (int z = 0; z < l.zdim; z++){
                  xyplaquettevalue = l.links[x][y][z][1].adjoint() * l.links[x][(y+1)%l.ydim][z][0].adjoint() * l.links[(x+1)%l.xdim][y][z][1] * l.links[x][y][z][0];
                  plaquette += xyplaquettevalue.trace().real();

            }}}
            //Beginne Metroplolis
      for(int i =0; i<rep;i++){

          //Gehe linear durchs Gitter
             for(int x = 0; x < l.xdim; x++){
                for (int y = 0; y < l.ydim; y++){
                    for (int z = 0; z < l.zdim; z++){
            //Erzeuge Zufallsmatrix für x-Richtung
                        test=zuf();
                        //Berechne die alten (alt) und neuen (neu) Plaquette-Werte
                        temp=(l.links[x][y][z][0]*l.links[(x+1)%l.xdim][y][z][1]*l.links[x][(y+1)%l.ydim][z][0].adjoint()*l.links[x][y][z][1].adjoint())+(l.links[x][y][z][0]*l.links[(x+1)%l.xdim][(y-1)%l.ydim][z][1].adjoint()*l.links[x][(y-1)%l.ydim][z][0].adjoint()*l.links[x][(y-1)%l.ydim][z][1]);
                        alt=temp.trace().real();
                        temp=(test*l.links[(x+1)%l.xdim][y][z][1]*l.links[x][(y+1)%l.ydim][z][0].adjoint()*l.links[x][y][z][1].adjoint())+(test*l.links[(x+1)%l.xdim][(y-1)%l.ydim][z][1].adjoint()*l.links[x][(y-1)%l.ydim][z][0].adjoint()*l.links[x][(y-1)%l.ydim][z][1]);
                        neu=temp.trace().real();
                        //Akzeptieren oder ablehnen
                        if(exp((-beta/N)*(neu-alt))>=1){
                            l.update(x,y,z,0,test);
                            plaquette=plaquette-alt+neu;
                        }
                        else{
                           p=uni_real_dist(Generator);
                           if (exp((-beta/N)*(neu-alt))>=p){
                                l.update(x,y,z,0,test);
                                plaquette=plaquette-alt+neu;

                           }
                        }

                        // Das Gleiche für y
                         test=zuf();
                        //Berechne die alten (alt) und neuen (neu) Plaquette-Werte
                        temp=(l.links[x][y][z][1]*l.links[x][(y+1)%l.ydim][z][0]*l.links[(x+1)%l.xdim][y][z][1].adjoint()*l.links[x][y][z][0].adjoint())+(l.links[x][y][z][1]*l.links[(x-1)%l.xdim][(y+1)%l.ydim][z][0].adjoint()*l.links[(x-1)%l.xdim][y][z][1].adjoint()*l.links[(x-1)%l.xdim][y][z][0]);
                        alt=temp.trace().real();
                        temp=(test*l.links[x][(y+1)%l.ydim][z][0]*l.links[(x+1)%l.xdim][y][z][1].adjoint()*l.links[x][y][z][0].adjoint())+(test*l.links[(x-1)%l.xdim][(y+1)%l.ydim][z][0].adjoint()*l.links[(x-1)%l.xdim][y][z][1].adjoint()*l.links[(x-1)%l.xdim][y][z][0]);
                        neu=temp.trace().real();
                        //Akzeptieren oder ablehnen
                        if(exp((-beta/N)*(neu-alt))>=1){
                            l.update(x,y,z,1,test);
                            plaquette=plaquette-alt+neu;
                        }
                        else{
                           p=uni_real_dist(Generator);
                           if (exp((-beta/N)*(neu-alt))>=p){
                                l.update(x,y,z,1,test);
                                plaquette=plaquette-alt+neu;

                           }
                        }


                    }}}


//Messwert speichern
 f <<plaquette<< endl;
      }
f.close();
    return 0;
}
