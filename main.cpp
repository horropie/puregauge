#include <iostream>
#include <random>
#include <math.h>
#include <fstream>

using namespace std;
double f;
fstream file;

std::default_random_engine generator;
std::uniform_real_distribution<double> uni_real_distr(0.0,1.0);



double function(double x)
{
    f= 3*x*x;
    return  f;
}

double rectint(double start,double finish, int step)
{
    double area=0.0;
    double interval= finish-start;
    double addstep=interval/step;
    //cout<< addstep<<endl;
    for(int i=0;i<step; i=i+1 ){
            //cout << start+(i*addstep)<<endl;
        area += addstep*(function(start+(i*addstep)));
    }
    return area;
}
double simpson(double start,double finish, int step){
    double area=0.0;
    double interval= finish-start;
    double addstep=interval/step;

    for(int i=1;i<=step/2;i++){
        area+=(addstep/3.0)*(function(start+addstep*(2*i-2))+(4*function(start+addstep*(2*i-1)))+function(start+addstep*(2*i)));
    }
    return area;
}

double montecarlo(int rep, int dimension){

    double in=0.0;
    double area=0.0;
    //double rand[dimension];
    double sum=0.0;


    for(int i =0; i<rep;i++){
            for(int j=0;j<dimension;j++){
              //rand[j]= ;
               sum=sum+pow(uni_real_distr(generator),2);
            }

        if(sum<=1)
            in=in+1.0;
        sum=0.0;

    }
    area=pow(2.0,dimension)*(in/rep);
    return area;

}
int faculty(int number){
    if(number ==1)
        return number;
    else
        return number*faculty(number-1);

}
 int double_fac(int number){
    if(number ==1)
        return number;
    else
        return number*double_fac(number-2);

 }

int main()
{
  double calculate_even;
  double calculate_odd;
  int dim_odd;
  int dim_even;

    //cout << montecarlo(1000,4)<< endl;
    string datei= "montecarlo.txt";
        file.open(datei,ios::out);
       file << "dimension;calculated;compare" << endl;
       for(int l=1;l<=10;l=l+1){
            dim_even=2*l;
            dim_odd=(2*l)+1;
         calculate_even=pow(3.1415,l)/faculty(l);
         calculate_odd=(pow(2,l+1)*pow(3.1415,l))/double_fac(dim_odd);
         //cout << calculate_even<<endl;
         //cout<< calculate_odd<<endl;

        file<<to_string(dim_even)+";"+to_string(montecarlo(10000,dim_even))+";"+to_string(calculate_even)<<endl;
        file<<to_string(dim_odd)+";"+to_string(montecarlo(10000,dim_odd))+";"+to_string(calculate_odd)<<endl;
       }
       file.close();
        //cout<<double_fac(5);


    return 0;
}
