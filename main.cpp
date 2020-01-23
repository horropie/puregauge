#include <iostream>
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
static int rep=500;
static int therm=100;
double beta;
double N =1.0;
double S_alt;
double S_neu;

double Vol;

Eigen::Matrix<std::complex<double>,2,2> test;
Eigen::Matrix<std::complex<double>,2,2> temp;



fstream f;

double deltaS;
double p;

Eigen::Matrix<std::complex<double>,2,2> polyakovLine;
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
	//cout<<"rand: "<<rand[0]<<" , "<<rand[1]<<" , "<<rand[2]<<endl;
	const double phase(2.0*3.1415*uni_real_dist(Generator));
	//cout<<"phase: "<<phase<<endl;
	const double cosn(cos(phase));
	const double sinn(sin(phase));
	//cout<<"cosn: "<<cosn<<endl;
	//cout<<"sinn: "<< sinn<<endl;
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
    s+= l.links[(x+1)%l.xdim][y][modul(z-1,l.zdim)][2][t].adjoint()*l.links[x][y][modul((z-1),l.zdim)][0][t].adjoint()*l.links[x][y][modul(z-1,l.zdim)][2][t];
// in xt-Richtung
    s+=l.links[(x+1)%l.xdim][y][z][3][t]*l.links[x][y][z][0][modul(t+1,l.tdim)].adjoint()*l.links[x][y][z][3][t].adjoint();
    s+= l.links[(x+1)%l.xdim][y][z][3][modul(t-1,l.tdim)].adjoint()*l.links[x][y][z][0][modul(t-1,l.tdim)].adjoint()*l.links[x][y][z][3][modul(t-1,l.tdim)];
    return s;
}
Eigen::Matrix<std::complex<double>,2,2> staple_y(lattice l,int x,int y, int z, int t){
        Eigen::Matrix<std::complex<double>,2,2> s;
        //in yx-Richtung
    s=l.links[x][(y+1)%l.ydim][z][0][t]*l.links[(x+1)%l.xdim][y][z][1][t].adjoint()*l.links[x][y][z][0][t].adjoint();
    s+=l.links[modul(x-1,l.xdim)][(y+1)%l.ydim][z][0][t].adjoint()*l.links[modul(x-1,l.xdim)][y][z][1][t].adjoint()*l.links[modul(x-1,l.xdim)][y][z][0][t];
        //in yz Richtung
    s+=l.links[x][(y+1)%l.ydim][z][2][t]*l.links[x][y][(z+1)%l.zdim][1][t].adjoint()*l.links[x][y][z][2][t].adjoint();
    s+=(l.links[x][(y+1)%l.ydim][modul(z-1,l.zdim)][2][t].adjoint()*l.links[x][y][modul(z-1,l.zdim)][1][t].adjoint()*l.links[x][y][modul(z-1,l.zdim)][2][t]);
 //in yt Richtung
    s+=l.links[x][modul(y+1,l.ydim)][z][3][t]*l.links[x][y][z][1][modul(t+1,l.tdim)].adjoint()*l.links[x][y][z][3][t].adjoint();
    s+= l.links[x][modul(y+1,l.ydim)][z][3][modul(t-1,l.tdim)].adjoint()*l.links[x][y][z][1][modul(t-1,l.tdim)].adjoint()*l.links[x][y][z][3][modul(t-1,l.tdim)];
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
   //in zt Richtung
    s+=l.links[x][y][modul(z+1,l.zdim)][3][t]*l.links[x][y][z][2][modul(t+1,l.tdim)].adjoint()*l.links[x][y][z][3][t].adjoint();
    s+= l.links[x][y][modul(z+1,l.zdim)][3][modul(t-1,l.tdim)].adjoint()*l.links[x][y][z][2][modul(t-1,l.tdim)].adjoint()*l.links[x][y][z][3][modul(t-1,l.tdim)];
    return s;
}

Eigen::Matrix<std::complex<double>,2,2> staple_t(lattice l,int x,int y, int z, int t){
        Eigen::Matrix<std::complex<double>,2,2> s;
        //in tx-Richtung
    s=l.links[x][y][z][0][modul(t+1,l.tdim)]*l.links[(x+1)%l.xdim][y][z][3][t].adjoint()*l.links[x][y][z][0][t].adjoint();
    s+=l.links[modul(x-1,l.xdim)][y][z][0][modul(t+1,l.tdim)].adjoint()*l.links[modul(x-1,l.xdim)][y][z][3][t].adjoint()*l.links[modul(x-1,l.xdim)][y][z][0][t];
        //in ty Richtung
    s+=l.links[x][y][z][1][modul(t+1,l.tdim)]*l.links[x][(y+1)%l.ydim][z][3][t].adjoint()*l.links[x][y][z][1][t].adjoint();
    s+=l.links[x][modul(y-1,l.ydim)][z][1][modul(t+1,l.tdim)].adjoint()*l.links[x][modul(y-1,l.ydim)][z][3][t].adjoint()*l.links[x][modul(y-1,l.ydim)][z][1][t];
   //in tz Richtung
    s+=l.links[x][y][z][2][modul(t+1,l.tdim)]*l.links[x][y][modul(z+1,l.zdim)][3][t].adjoint()*l.links[x][y][z][2][t].adjoint();
    s+= l.links[x][y][modul(z-1,l.zdim)][2][modul(t+1,l.tdim)].adjoint()*l.links[x][y][modul(z-1,l.zdim)][3][t].adjoint()*l.links[x][y][modul(z-1,l.zdim)][2][t];
    return s;
}

double plaquette_sum(lattice l){
    Eigen::Matrix<std::complex<double>,2,2> xyplaquettevalue, xzplaquettevalue, yzplaquettevalue,xtplaquettevalue,ytplaquettevalue,ztplaquettevalue;

    double O=0.0;
   for (int t = 0; t<l.tdim; t++){
      for(int x = 0; x < l.xdim; x++){
        for (int y = 0; y < l.ydim; y++){
          for (int z = 0; z < l.zdim; z++){
                xyplaquettevalue = (l.links[x][y][z][0][t]*l.links[(x+1)%l.xdim][y][z][1][t]*l.links[x][modul(y+1,l.ydim)][z][0][t].adjoint()*l.links[x][y][z][1][t].adjoint());
                xzplaquettevalue = (l.links[x][y][z][0][t]*l.links[(x+1)%l.xdim][y][z][2][t]*l.links[x][y][(z+1)%l.zdim][0][t].adjoint()*l.links[x][y][z][2][t].adjoint());
                yzplaquettevalue = (l.links[x][y][z][1][t]*l.links[x][(y+1)%l.ydim][z][2][t] * l.links[x][y][(z+1)%l.zdim][1][t].adjoint() *l.links[x][y][z][2][t].adjoint());

                xtplaquettevalue= (l.links[x][y][z][3][t]*l.links[x][y][z][0][modul(t+1,l.tdim)]*l.links[modul(x+1,l.xdim)][y][z][3][t].adjoint()*l.links[x][y][z][0][t].adjoint());
                ytplaquettevalue= (l.links[x][y][z][3][t]*l.links[x][y][z][1][modul(t+1,l.tdim)]*l.links[x][modul(y+1,l.ydim)][z][3][t].adjoint()*l.links[x][y][z][1][t].adjoint());
                ztplaquettevalue= (l.links[x][y][z][3][t]*l.links[x][y][z][2][modul(t+1,l.tdim)]*l.links[x][y][modul(z+1,l.zdim)][3][t].adjoint()*l.links[x][y][z][2][t].adjoint());

                O += (xyplaquettevalue.trace().real()+xzplaquettevalue.trace().real()+yzplaquettevalue.trace().real()+xtplaquettevalue.trace().real()+ytplaquettevalue.trace().real()+ztplaquettevalue.trace().real());

}}}}

            return O;
}
double action(lattice l,double inv_coupling){
return inv_coupling*((2*3*l.xdim*l.ydim*l.zdim*l.tdim)-plaquette_sum(l));
}




int main()
{
    lattice l(8,4) ;
    Vol=l.xdim*l.ydim*l.zdim*l.tdim;
for(int b=1;b<=100;b++){
          beta=b*0.05;
    //Initialisiere Gitter mit heißem Start

    // Berechne die Summe der Plaquettes des Gitters, speichere in plaquette
   double plaquette=beta*plaquette_sum(l);


   //
   // cout<<"Gesamt:"<<plaquette<<endl;

            //Speicherort und Name der Ausgabedatei
             string datei= "Metropolis"+to_string(beta)+".txt";
            //string datei="test.txt";
            f.open(datei,ios::out);
            f << "O ; Polaykov"<<endl;

        //Beginne Metropolis
        for(int i =0; i<rep+therm;i++){

            //Gehe linear durchs Gitter
        for( int t = l.tdim-1;t>=0;t-- ){
            for(int x = 0; x < l.xdim; x++){
                for (int y = 0; y < l.ydim; y++){
                    for (int z = 0; z < l.zdim; z++){
                        //Erzeuge Zufallsmatrix für x-Richtung
                        zuf(test);
                        //Berechne die Differenz der alten (alt) und neuen (neu) Plaquette-Werte (neu-alt)
                        // Alte Plaquettes in xy Richtung
                        deltaS=beta*((-l.links[x][y][z][0][t]+test)*staple_x(l,x,y,z,t)).trace().real();
                        //Akzeptieren oder ablehnen
                        p=uni_real_dist(Generator);
                        if (exp(-deltaS)>=p){
                            l.update(x,y,z,0,t,test);
                            plaquette=plaquette+(deltaS);
                        }

                        // Das Gleiche für y
                        zuf(test);
                        //Berechne die Differenz der alten (alt) und neuen (neu) Plaquette-Werte (neu-alt)
                        // Alte Werte in xy
                        deltaS=beta*((test-l.links[x][y][z][1][t])*staple_y(l,x,y,z,t)).trace().real();
                        p=uni_real_dist(Generator);
                        if (exp(-deltaS)>=p){
                            l.update(x,y,z,1,t,test);
                            plaquette=plaquette+(deltaS);
                        }


                        // Das Gleiche für z
                         zuf(test);
                        //Berechne die alten (alt) und neuen (neu) Plaquette-Werte
                        deltaS=beta*((test-l.links[x][y][z][2][t])*staple_z(l,x,y,z,t)).trace().real();
                        p=uni_real_dist(Generator);
                        if (exp(-deltaS)>=p){
                            l.update(x,y,z,2,t,test);
                            plaquette=plaquette+(deltaS);
                        }

                        // Das Gleiche für t
                         zuf(test);
                        //Berechne die alten (alt) und neuen (neu) Plaquette-Werte
                        deltaS=beta*((test-l.links[x][y][z][3][t])*staple_t(l,x,y,z,t)).trace().real();
                        p=uni_real_dist(Generator);
                        if (exp(-deltaS)>=p){
                            l.update(x,y,z,3,t,test);
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
       f << plaquette/(24*Vol)<<" ; "<<polyakovLineSum/(2*Vol) << endl;
     polyakovLineSum=0.0;
     }

f.close();

}
//Speicher freigeben
        for(int i =0;i<l.xdim;i++){
            for(int j =0;j<l.ydim;j++){
                for(int k=0;k<l.zdim;k++){
                    for(int r=0;r<4;r++){
                       delete[] l.links[i][j][k][r];
                    }


                delete[] l.links[i][j][k];
                }
                delete[] l.links[i][j];
                }
                delete[] l.links[i];
                }
                delete[] l.links;


return 0;
}
