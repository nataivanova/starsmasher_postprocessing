#include<iostream>
#include<cmath>

const double mol=1.0,g=1.0,erg=1.0,K=1.0,eV=1.0,s=1.0;
const double NA=6.02214179e23/mol;
const double me=9.1093897e-28*g;
const double kcgs=1.380658e-16*erg/K;
const double chiH0=13.598*eV;
const double chiHe0=24.587*eV;
const double chiHe1=54.416*eV;
const double keV=8.617343e-5*eV/K;
const double hp=6.6260755e-27*erg/s;

using namespace std;

double nion(double rho,double X, double Y){return rho*NA*(X+Y/4.0);}

double fH0(double T){
    return 2.0*pow(2.0*M_PI*me*kcgs*T,3.0/2.0)*exp(-chiH0/(keV*T)) / (hp*hp*hp);
}

double fHe0(double T){
    return 2.0*pow(2.0*M_PI*me*kcgs*T,3.0/2.0)*exp(-chiHe0/(keV*T)) / (hp*hp*hp);
}

double fHe1(double T){
        return 2.0*pow(2.0*M_PI*me*kcgs*T,3.0/2.0)*exp(-chiHe1/(keV*T)) / (hp*hp*hp);
}

double XH1(double rho,double T,double X,double Y,double xe){
    return 0.5*fH0(T)/nion(rho,X,Y)/(xe+0.5*fH0(T)/nion(rho,X,Y));
}

double XHe1(double rho,double T,double X,double Y,double xe){
    return 2.0*fHe0(T)/nion(rho,X,Y)/( xe + 2.0*fHe0(T)/nion(rho,X,Y) + 
               fHe0(T)*fHe1(T)/(nion(rho,X,Y)*nion(rho,X,Y)*xe) );
}

double XHe2(double rho,double T,double X,double Y,double xe){
    return 0.5 * XHe1(rho,T,X,Y,xe) / xe * fHe1(T) / nion(rho,X,Y);
}

double Xe(double rho,double T,double X,double Y,double xe){
    return (X*XH1(rho,T,X,Y,xe) + Y/4.0 * (XHe1(rho,T,X,Y,xe) + 2.0 * XHe2(rho,T,X,Y,xe)))/(X + Y/4.0);
}

// nata: note that this routine is not written well.
// I placed templorary cut off, but it does not work well at neat neutral ionization, didnt check whether
// it works for a full ionization at rho>1 g/sm^3.
void Solve(double rho,double T,double X,double Y,
             double &xh1n,double &xhe1n,double &xhe2n,double &xen,int &iter,double xe0=0.01){
    int i;
    for(i=1;i<10000;i++){
      if(T > 1000) {
	xh1n = XH1(rho,T,X,Y,xe0);
	xhe1n = XHe1(rho,T,X,Y,xe0);
	xhe2n = XHe2(rho,T,X,Y,xe0);
        xen = Xe(rho,T,X,Y,xe0);}      
      else {
	xh1n = 0.; xhe1n = 0.;  xhe2n = 0.; xen=0.; 
      }
        if(abs(xe0-xen) <= 1e-10 ) {iter=i;return;}
        else xe0 =xen;
    }
    iter=i;
}

extern "C"
{
    extern void deg_of_ion_(const double *rho,const double *T,const double *X,const double *Y,
                            double *xe,double *xh1,double *xhe1,double *xhe2,int *iter)
    {
        //cout<<"Input in C++: "<<*rho<<" "<<*T<<" "<<*X<<" "<<*Y<<endl;
        Solve(*rho,*T,*X,*Y,*xh1,*xhe1,*xhe2,*xe,*iter);
        //cout<<"Output in C++: "<<*xe<<endl;
    }
}

