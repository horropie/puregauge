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
double polyakovLine;
double polyakovLineSum;



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


Eigen::Matrix<std::complex<double>,2,2> staple_x(lattice l,int x,int y, int z, int t){

    Eigen::Matrix<std::complex<double>,2,2> s;
    //in xy-Richtung
    s=l.links[(x+1)%l.xdim][y][z][1][t]*l.links[x][(y+1)%l.ydim][z][0][t].adjoint()*l.links[x][y][z][1][t].adjoint();
    s+=l.links[(x+1)%l.xdim][modul(y-1,l.ydim)][z][1][t].adjoint()*l.links[x][modul(y-1,l.ydim)][z][0][t].adjoint()*l.links[x][modul(y-1,l.ydim)][z][1][t];
    //in xz Richtung
  s+=l.links[(x+1)%l.xdim][y][z][2][t]*l.links[x][y][(z+1)%l.zdim][0][t].adjoint()*l.links[x][y][z][2][t].adjoint();
      s+= l.links[x][y][z][0][t]*l.links[(x+1)%l.xdim][y][modul(z-1,l.zdim)][2][t].adjoint()*l.links[x][y][modul((z-1),l.zdim)][0][t].adjoint()*l.links[x][y][modul(z-1,l.zdim)][2][t];

     return s;
}
Eigen::Matrix<std::complex<double>,2,2> staple_y(lattice l,int x,int y, int z, int t){
        Eigen::Matrix<std::complex<double>,2,2> s;
        //in yx-Richtung
        s=l.links[x][(y+1)%l.ydim][z][0][t]*l.links[(x+1)%l.xdim][y][z][1][t].adjoint()*l.links[x][y][z][0][t].adjoint();
        s+=l.links[modul(x-1,l.xdim)][(y+1)%l.ydim][z][0][t].adjoint()*l.links[modul(x-1,l.xdim)][y][z][1][t].adjoint()*l.links[modul(x-1,l.xdim)][y][z][0][t];
        //in yz Richtung
        s+=l.links[x][(y+1)%l.ydim][z][2][t]*l.links[x][y][(z+1)%l.zdim][1][t].adjoint()*l.links[x][y][z][2][t].adjoint();
        s+=(l.links[x][y][z][1][t]*l.links[x][(y+1)%l.ydim][modul(z-1,l.zdim)][2][t].adjoint()*l.links[x][y][modul(z-1,l.zdim)][1][t].adjoint()*l.links[x][y][modul(z-1,l.zdim)][2][t]);
    return s;
    }

Eigen::Matrix<std::complex<double>,2,2> staple_z(lattice l,int x,int y, int z, int t){
        Eigen::Matrix<std::complex<double>,2,2> s;
        //in zx-Richtung
        s=l.links[x][y][(z+1)%l.zdim][0][t]*l.links[(x+1)%l.xdim][y][z][2][t].adjoint()*l.links[x][y][z][0][t].adjoint();
        s+=l.links[modul(x-1,l.xdim)][y][(z+1)%l.zdim][0][t].adjoint()*l.links[modul(x-1,l.xdim)][y][z][2][t].adjoint()*l.links[modul(x-1,l.xdim)][y][z][0][t];
        //in zy Richtung
        s+=l.links[x][y][(z+1)%l.zdim][1][t]*l.links[x][(y+1)%l.ydim][z][2][t].adjoint()*l.links[x][y][z][1][t].adjoint();
        s+=l.links[x][modul(y-1,l.ydim)][(z+1)%l.zdim][1][t].adjoint()*l.links[x][modul(y-1,l.ydim)][z][2][t].adjoint()*l.links[x][modul(y-1,l.ydim)][z][1][t];
    return s;
}

double action(lattice l){
    Eigen::Matrix<std::complex<double>,2,2> xyplaquettevalue, xzplaquettevalue, yzplaquettevalue;
    double O=0.0;
   for (int t = 0; t<l.tdim; t++){
      for(int x = 0; x < l.xdim; x++){
        for (int y = 0; y < l.ydim; y++){
          for (int z = 0; z < 1; z++){
                  xyplaquettevalue = l.links[x][y][z][1][t].adjoint() * l.links[x][modul(y+1,l.ydim)][z][0][t].adjoint() * l.links[(x+1)%l.xdim][y][z][1][t] * l.links[x][y][z][0][t];
                  xzplaquettevalue = l.links[x][y][z][2][t].adjoint() * l.links[x][y][(z+1)%l.zdim][0][t].adjoint() * l.links[(x+1)%l.xdim][y][z][2][t] * l.links[x][y][z][0][t];
                  yzplaquettevalue = l.links[x][y][z][2][t].adjoint() * l.links[x][y][(z+1)%l.zdim][1][t].adjoint() * l.links[x][(y+1)%l.ydim][z][2][t] * l.links[x][y][z][1][t];
                  O += (xyplaquettevalue.trace().real()+xzplaquettevalue.trace().real()+yzplaquettevalue.trace().real());

            }}}}
            //O=O/(2*6*6*6*4);
            return O;
}





int main()
{

    //Initialisiere Gitter mit heißem Start
    lattice l ;
    // Berechne die Summe der Plaquettes des Gitters, speichere in plaquette
    double plaquette=0.0;
     beta=1.0;
            //Speicherort und Name der Ausgabedatei
            // string datei= "Metropolis"+to_string(t)+".txt";
            string datei="test.txt";
            f.open(datei,ios::out);
            f << "O; action; Polyakov" << endl;

        //Beginne Metropolis
        for(int i =0; i<rep+therm;i++){

            //Gehe linear durchs Gitter
        for( int t = l.tdim-1;t>=0;t-- ){
            for(int x = 0; x < l.xdim; x++){
                for (int y = 0; y < l.ydim; y++){
                    for (int z = 0; z < 1; z++){


                        //Erzeuge Zufallsmatrix für x-Richtung
                        zuf(test);

                        // cout<<"Spur d. Testmatrix:"<<test.trace()<<endl;

                        //Berechne die Differenz der alten (alt) und neuen (neu) Plaquette-Werte (neu-alt)
                        // Alte Plaquettes in xy Richtung
                        deltaS=((-l.links[x][y][z][0][t]+test)*staple_x(l,x,y,z,t)).trace().real();
                        //Akzeptieren oder ablehnen
                        p=uni_real_dist(Generator);
                        if (exp((-beta/N)*(deltaS))>=p){
                            l.update(x,y,z,0,t,test);
                            plaquette=plaquette+(deltaS);
                           }


                        // Das Gleiche für y
                        zuf(test);
                        // cout<<"Spur d. Testmatrix:"<<test.trace().real()<<endl;
                        //Berechne die Differenz der alten (alt) und neuen (neu) Plaquette-Werte (neu-alt)
                        // Alte Werte in xy
                        deltaS=((test-l.links[x][y][z][1][t])*staple_y(l,x,y,z,t)).trace().real();
                        p=uni_real_dist(Generator);
                        if (exp((-beta/N)*deltaS)>=p){
                            l.update(x,y,z,1,t,test);
                             plaquette=plaquette+(deltaS);

                           }


                        // Das Gleiche für z
                         zuf(test);
                        // cout<<"Spur d. Testmatrix:"<<test.trace().real()<<endl;
                        //Berechne die alten (alt) und neuen (neu) Plaquette-Werte
                        deltaS=((test-l.links[x][y][z][2][t])*staple_z(l,x,y,z,t)).trace().real();
                        p=uni_real_dist(Generator);
                        if (exp((-beta/N)*deltaS)>=p){
                            l.update(x,y,z,2,t,test);
                             plaquette=plaquette+(deltaS);

                            }



            }}}}

    //Polyakov Line
    for(int x = 0; x < l.xdim; x++){
      for (int y = 0; y < l.ydim; y++){
        for (int z = 0; z < 1; z++){
          for (int t = l.tdim-1; t >= 0; t--){

             // cout << (l.links[x][y][z][0][t]).trace().real() << endl;
              polyakovLine *= l.links[x][y][z][0][t];
              }
              //Auf Identität setzen und Spur rausziehen

              polyakovLineSum += polyakovLine.trace().real();
              polyakovLine.setIdentity();

          }}}

    //cout <<"Polyakov:"<< polyakovLineSum << endl;

    //Messwert speichern
    if(i>=therm)
        f <<plaquette<<" ; "<< action(l)<<" ; "<<polyakovLineSum << endl;
     polyakovLineSum=0.0;
     }

f.close();
//cout<<l.links[0][0][0][0][0]<<endl;

return 0;
}
