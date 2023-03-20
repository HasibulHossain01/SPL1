#include<bits/stdc++.h>
 using namespace std;
 //dx/dt = f(t,x,y)=x+2+y;
 //dy/dt=  g(t,x,y)=3+x+2+y;

 double f(const double&t,const double &x,const double &y);
 double g(const double&t,const double &x,const double &y);

 void Runge_Kutta(vector<double>&x,vector<double>&y,vector<double>&t,const double &h,const int& n){
    double F1,F2,F3,F4,G1,G2,G3,G4;
    for(int i=0;i<n;++i){
        F1=f(t[i],x[i],y[i]);
        G1=g(t[i],x[i],y[i]);
        F2=f(t[i]+(h/2.),x[i]+F1*(h/2.),y[i]+G1*(h/2.));
        G2=g(t[i]+(h/2.),x[i]+F1*(h/2.),y[i]+G1*(h/2.));
        F3=f(t[i]+(h/2.),x[i]+F2*(h/2.),y[i]+G2*(h/2.));
        G3=g(t[i]+(h/2.),x[i]+F2*(h/2.),y[i]+G2*(h/2.));
        F4=f(t[i]+h,x[i]+F3*h,y[i]+G3*h);
        G4=g(t[i]+h,x[i]+F3*h,y[i]+G3*h);

        x[i+1]=x[i]+(h/6)*(F1+2*F2+2*F3+F4);
        y[i+1]=x[i]+(h/6)*(G1+2*G2+2*G3+G4);
        t[i+1]=t[i]+h;
        
    }
 }
 
 int main (){
   int n ; //no of iteration::specified bt the user.
   std::cout<<"Please inter the number of iteration"<<std::endl;
   std::cin>>n;
   vector<double>t[n];
   vector<double>x[n];
   vector<double>y[n];
   double h=0.02;//step size

   t[0]=0;
   x[0]=6;
   y[0]=4;
   
}

 double f(const double&t,const double &x,const double &y){

    return x+2*y;
 }
 double g(const double&t,const double &x,const double &y){
    return 3*x+y;
 }
