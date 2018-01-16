#include<iostream>
#include<cmath>

#define G 1

using namespace std;

double phi(double x,double y,double z,double m1,double m2,double a)
{
   double m=m1+m2;
   double phi=-G*m1/sqrt(pow(x+m2/m*a,2)+y*y+z*z)-G*m2/sqrt(pow(x-m1/m*a,2)+y*y+z*z)
             -G*m/(2.0*a*a*a)*(x*x+y*y);
   return phi;
}

double phi_x(double x,double y,double z,double m1,double m2,double a)
{
   double m=m1+m2;
   double phi_x=G*m1*(x+m2/m*a)/pow(pow(x+m2/m*a,2)+y*y+z*z,3.0/2.0)+G*m2*(x-m1/m*a)/pow(pow(x-m1/m*a,2)+y*y+z*z,3.0/2.0)-G*m*x/(a*a*a);
   return phi_x;
}

double phi_xx(double x,double y,double z,double m1,double m2,double a)
{
   double m=m1+m2;
   double phi_xx=-3.0*G*m1*pow(x+m2/m*a,2.0)/pow(pow(x+m2/m*a,2)+y*y+z*z,5.0/2.0)+G*m1/pow(pow(x+m2/m*a,2)+y*y+z*z,3.0/2.0)-3.0*G*m2*pow(x-m1/m*a,2.0)/pow(pow(x-m1/m*a,2)+y*y+z*z,5.0/2.0)+G*m2/pow(pow(x-m1/m*a,2)+y*y+z*z,3.0/2.0)-G*m/(a*a*a);
   return phi_xx;
}

double phi_y(double x,double y,double z,double m1,double m2,double a)
{
   double m=m1+m2;
   double phi_y=G*m1*y/pow(pow(x+m2/m*a,2)+y*y+z*z,3.0/2.0)+G*m2*y/pow(pow(x-m1/m*a,2)+y*y+z*z,3.0/2.0)-G*m*y/(a*a*a);
   return phi_y;
}

double phi_z(double x,double y,double z,double m1,double m2,double a)
{
   double m=m1+m2;
   double phi_z=G*m1*z/pow(pow(x+m2/m*a,2)+y*y+z*z,3.0/2.0)+G*m2*z/pow(pow(x-m1/m*a,2)+y*y+z*z,3.0/2.0);
   return phi_z;
}

double newton_raphson(double f(double,double,double,double,double,double),double fx(double,double,double,double,double,double),double x0,double x,double y,double z,double m1,double m2,double a)
{
  double x1;
  for(;;){
     x1 = x0 - f(x0,y,z,m1,m2,a)/fx(x0,y,z,m1,m2,a);
     if (abs(x1-x0)<1.e-10) break;
     x0=x1;
  }
  return x0;
}


double Phi(double m1,double m2,double a){
   double phi1,root;
   double x0=0.0;
   root = newton_raphson(&phi_x,&phi_xx,x0,0,0,0,m1,m2,a);
   phi1 = phi(root,0,0,m1,m2,a);
   return phi1;
}

extern "C" 
{
   extern void phi_(double *phi1,double *xL1,double *m1,double *m2, double *a)
      {
           *phi1 = Phi(*m1,*m2,*a);
           double x0=0.0;
           *xL1 = newton_raphson(&phi_x,&phi_xx,x0,0,0,0,*m1,*m2,*a);
      }
}
