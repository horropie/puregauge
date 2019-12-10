#include <iostream>
#include <iostream>
#include "lattice.cpp"

#include <eigen3/Eigen/Dense>
#include <random>
#include <vector>
#include <complex>
//#include <boost/math/constants/constants.hpp>
#include <fstream>
using namespace std;

//definiere Parameter (Anzahl Wdh., Temperatur, Zufallsgenerator)
static int rep=100;
static int therm=200;
double beta;
double N =1.0;

Eigen::Matrix<std::complex<double>,2,2> test;
Eigen::Matrix<std::complex<double>,2,2> temp;

fstream f;

double deltaS;
double p;



std::default_random_engine Generator;
std::uniform_real_distribution<double> uni_real_dist(0.0,1.0);

// Methode zum Erzeugen von Zufallsmatrix
 void zuf(Eigen::Matrix<std::complex<double>,2,2>& randm)
{

  double rand[3];
	double norm=2.0;
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

}

int modul(int x, int mod){
 x=x%mod;
 if(x<0){
    x=mod+x;
 }
 return x;
}


Eigen::Matrix<std::complex<double>,2,2> staple_x(lattice l,int x,int y, int z){
    Eigen::Matrix<std::complex<double>,2,2> s;
    //in xy-Richtung
    s=l.links[(x+1)%l.xdim][y][z][1]*l.links[x][(y+1)%l.ydim][z][0].adjoint()*l.links[x][y][z][1].adjoint();
    s+=l.links[(x+1)%l.xdim][modul(y-1,l.ydim)][z][1].adjoint()*l.links[x][modul(y-1,l.ydim)][z][0].adjoint()*l.links[x][modul(y-1,l.ydim)][z][1];
    //in xz Richtung
    //s+=l.links[(x+1)%l.xdim][y][z][2]*l.links[x][y][(z+1)%l.zdim][0].adjoint()*l.links[x][y][z][2].adjoint();
    //s+= l.links[x][y][z][0]*l.links[(x+1)%l.xdim][y][modul(z-1,l.zdim)][2].adjoint()*l.links[x][y][(z-1)%l.zdim][0].adjoint()*l.links[x][y][modul(z-1,l.zdim)][2];
     return s;
}
Eigen::Matrix<std::complex<double>,2,2> staple_y(lattice l,int x,int y, int z){
        Eigen::Matrix<std::complex<double>,2,2> s;
        //in yx-Richtung
        s=l.links[x][(y+1)%l.ydim][z][0]*l.links[(x+1)%l.xdim][y][z][1].adjoint()*l.links[x][y][z][0].adjoint();
        s+=l.links[modul(x-1,l.xdim)][(y+1)%l.ydim][z][0].adjoint()*l.links[modul(x-1,l.xdim)][y][z][1].adjoint()*l.links[modul(x-1,l.xdim)][y][z][0];
        //in yz Richtung
        //s+=l.links[x][(y+1)%l.ydim][z][2]*l.links[x][y][(z+1)%l.zdim][1].adjoint()*l.links[x][y][z][2].adjoint();
        //s+=(l.links[x][y][z][1]*l.links[x][(y+1)%l.ydim][modul(z-1,l.zdim)][2].adjoint()*l.links[x][y][modul(z-1,l.zdim)][1].adjoint()*l.links[x][y][modul(z-1,l.zdim)][2]);
    return s;
    }

Eigen::Matrix<std::complex<double>,2,2> staple_z(lattice l,int x,int y, int z){
        Eigen::Matrix<std::complex<double>,2,2> s;
        //in zx-Richtung
        //s=l.links[x][y][(z+1)%l.zdim][0]*l.links[(x+1)%l.xdim][y][z][2].adjoint()*l.links[x][y][z][0].adjoint();
        //s+=l.links[modul(x-1,l.xdim)][y][(z+1)%l.zdim][0].adjoint()*l.links[modul(x-1,l.xdim)][y][z][2].adjoint()*l.links[modul(x-1,l.xdim)][y][z][0];
        //in zy Richtung
        //s+=l.links[x][y][(z+1)%l.zdim][1]*l.links[x][(y+1)%l.ydim][z][2].adjoint()*l.links[x][y][z][1].adjoint();
        //s+=l.links[x][modul(y-1,l.ydim)][(z+1)%l.zdim][1].adjoint()*l.links[x][modul(y-1,l.ydim)][z][2].adjoint()*l.links[x][modul(y-1,l.ydim)][z][1];
    return s;
}





int main()
{

    //Initialisiere Gitter mit heißem Start
    lattice l ;
    // Berechne die Summe der Plaquettes des Gitters, speichere in plaquette
    double plaquette=0.0;
    Eigen::Matrix<std::complex<double>,2,2> xyplaquettevalue, xzplaquettevalue, yzplaquettevalue;

     for(int x = 0; x < l.xdim; x++){
        for (int y = 0; y < l.ydim; y++){
            for (int z = 0; z < 1; z++){
                  xyplaquettevalue = l.links[x][y][z][1].adjoint() * l.links[x][modul(y+1,l.ydim)][z][0].adjoint() * l.links[(x+1)%l.xdim][y][z][1] * l.links[x][y][z][0];
                  //xzplaquettevalue = l.links[x][y][z][2].adjoint() * l.links[x][y][(z+1)%l.zdim][0].adjoint() * l.links[(x+1)%l.xdim][y][z][2] * l.links[x][y][z][0];
                  //yzplaquettevalue = l.links[x][y][z][2].adjoint() * l.links[x][y][(z+1)%l.zdim][1].adjoint() * l.links[x][(y+1)%l.ydim][z][2] * l.links[x][y][z][1];
                  plaquette += (xyplaquettevalue.trace().real()); //+xzplaquettevalue.trace().real()+yzplaquettevalue.trace().real());

            }}}

    // Für verschiedene Temperaturen
    for( int t =20;t>19;t-- ){
            beta=1/(t*0.1);
            //Speicherort und Name der Ausgabedatei
            // string datei= "Metropolis"+to_string(t)+".txt";
            string datei="therm.txt";
            f.open(datei,ios::out);
            f << "O" << endl;

        //Beginne Metropolis
        for(int i =0; i<rep+therm;i++){

            //Gehe linear durchs Gitter
            for(int x = 0; x < l.xdim; x++){
                for (int y = 0; y < l.ydim; y++){
                    for (int z = 0; z < 1; z++){

                        //Erzeuge Zufallsmatrix für x-Richtung
                        zuf(test);

                        // cout<<"Spur d. Testmatrix:"<<test.trace()<<endl;

                        //Berechne die Differenz der alten (alt) und neuen (neu) Plaquette-Werte (neu-alt)
                        // Alte Plaquettes in xy Richtung
                        deltaS=((-l.links[x][y][z][0]+test)*staple_x(l,x,y,z)).trace().real();
                        //Akzeptieren oder ablehnen
                        p=uni_real_dist(Generator);
                        if (exp((-beta/N)*(deltaS))>=p){
                            l.update(x,y,z,0,test);
                            plaquette=plaquette+deltaS;
                           }


                        // Das Gleiche für y
                        zuf(test);
                        // cout<<"Spur d. Testmatrix:"<<test.trace().real()<<endl;
                        //Berechne die Differenz der alten (alt) und neuen (neu) Plaquette-Werte (neu-alt)
                        // Alte Werte in xy
                        deltaS=((test-l.links[x][y][z][1])*staple_y(l,x,y,z)).trace().real();
                        p=uni_real_dist(Generator);
                        if (exp((-beta/N)*deltaS)>=p){
                            l.update(x,y,z,1,test);
                            plaquette=plaquette+deltaS;

                           }


                        // Das Gleiche für z
                        // zuf(test);
                        // cout<<"Spur d. Testmatrix:"<<test.trace().real()<<endl;
                        //Berechne die alten (alt) und neuen (neu) Plaquette-Werte
                        //deltaS=((test-l.links[x][y][z][2])*staple_z(l,x,y,z)).trace().real();
                        //p=uni_real_dist(Generator);
                        //if (exp((-beta/N)*deltaS)>=p){
                        //    l.update(x,y,z,2,test);
                        //      plaquette=plaquette+deltaS;

                        //     }

            }}}
            
    //Messwert speichern
    if(i>=therm)
        f <<plaquette<< endl;
    }

f.close();
    }
return 0;
}
