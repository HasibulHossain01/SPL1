#include "bits/stdc++.h"
#include "matplotlibcpp.h"

#define MOD (1e9 + 7);
#define FIXED_FLOAT(x) std::fixed<<std::setprecision(2)<<(x)
using ll = int64_t;
using ull = uint64_t;
#define ll long long

namespace plt=matplotlibcpp;
using namespace std;

int inputFunctionK();
int inputFuntionU(int enoughData,double V,double t,double theta);
int printValues(double tx,double Vx,double Vy,double Ax,double Ay,double x,double y,double V);
int withFriction(double x,double y,double V,double t,double theta,double D);
int withoutFriction(double x,double y,double V,double t,double theta);
double func(double x);
int derivFunc(double x);
int distanceFromOrigin(double x,double y);
int preservationofEnergy(double V,double y);
void newtonRaphson(double x);
double lenthofthepath(double V,double theta,double t);
string diffofX(double V,double theta);
string diffofY(double V,double theta);
void EnergyGraph();
double rootFunc( double t);
double calculate(double lower_limit, double upper_limit,int interval_limit );
string diffTerm(string pTerm);
string diffstr(string& poly);

double D,V,p,C,A,r,distanceFromO;
double t,T,Tf,H,R,TotalEnergy,PotentialEnergy,KineticEnergy;
double x,y,theta,m;
const double g=9.81,pi=3.1416;
int enoughData=0;
char chF;

int main(void){

    cout<<" Is it a projectile motion with air resistance?(y/n): ";
    cin>>chF;
   
    inputFunctionK();
   
    if((enoughData==3)&&(chF='y')) withFriction(x,y,V,t,theta,D);
    
    else withoutFriction(x,y,V,t,theta);

     
     return 0;

}

int inputFunctionK(){
    cout<<"Enter initial Co-ordinate of x: ";
    cin>>x;
    if(cin.fail()) {
      cout<<" taking 0 as a ideal position";
      x=0;
    }
    cout<<endl;

    cout<<"Enter initial Co-ordinate of y(mandatory): ";
    cin>>y;
    if(cin.fail()) {
    cout<<" taking 0 as a ideal position";
      y=0;
    }
    cout<<endl;

    cout<<"Enter Initial Velocity V in m/s: ";
    cin>>V;
    if(V>0) enoughData++;
    cout<<endl;

    cout<<"Enter Time t in second: ";
    cin>>t;
    if(t>0) enoughData++;
    cout<<endl;

     cout<<"Enter Initial Angle in Degree theta: ";
     cin>>theta;
     if(theta>=0) enoughData++; //value of theta can be 0
     cout<<endl;
    
     cout<<"Enter mass of the object in kg m(mandatory): " ;
     cin>>m;
     cout<<endl;

     if(chF=='y'){
     cout<<"Enter Friction Coefficient of air D(not mandatory,can be extracted): " ;
     cin>>D;
     if(D==0){    
         cout<<"     Enter Density p(mandatory): ";
         cin>>p;
         cout<<"     Enter Drag-Coeeficient(mandatory): ";
         cin>>C;
         cout<<"     Enter radius(mandatory): ";
         cin>>r;
         A=pi*(r*r);
         D=(p*C*A)/2;
         if(p==0||C==0||r==0) {
          cout<<"Not enough Data...";
          EXIT_FAILURE;
          }
     }
    }
     cout<<endl;
     
     if(enoughData<3) {
        inputFuntionU(enoughData,V,t,theta);
     }

     theta=theta*(3.1416/180);
     return 0;
}

int inputFuntionU(int enoughData,double V,double t,double theta){
          double Vx=V*cos(theta);
          double  Vy=V*sin(theta);
        int enoughData_1=enoughData;
        cout<<"Enter any of the three of Time of flight(T), MAximum Vertical Range(H) or Maximum horizontal Range(R)[if unknown enter 0]: "<<endl;
        cout<<"T : ";
        cin>>Tf;
        if(Tf!=0) enoughData_1++;
        cout<<"H: ";
        cin>>H;
        if(H!=0) enoughData_1++;
        cout<<"R : ";
        cin>>R;
        if(R!=0) enoughData_1++;

        if(enoughData+enoughData_1>=3){
          if(enoughData==2&&enoughData_1>=1){
            
            if(V!=0&&t!=0&&T!=0){
              theta=asin((g*T)/(2*V));
            }
            if(V!=0&&theta!=0&&T!=0){
              t=T;
              enoughData++;
            }
            else if(theta!=0&&t!=0&&T!=0){
               V=(g*T)/(2*sin(theta));
               enoughData++;
            }

            else if(V!=0&&t!=0&&H!=0){
              theta=asin((2*g*sqrt(H))/V);
              enoughData++;
            }
            else if(V!=0&&theta!=0&&H!=0){
              t=(2*V*sin(theta)/g);
              enoughData++;
            }
            else if(theta!=0&&t!=0&&H!=0){
               V=((2*g*sqrt(H))/sin(theta));
               enoughData++;
            }

            else if(V!=0&&t!=0&&R!=0){
              theta=0.5*((asin((g*R)/(V*V))));
              enoughData++;
            }
            else if(V!=0&&theta!=0&&R!=0){
              t=(2*V*sin(theta)/g);
              enoughData++;
            }
            else if(theta!=0&&t!=0&&R!=0){
               V=sqrt((g*R)/sin(2*theta));
               enoughData++;
            }                      
           }
          else if(enoughData==1&&enoughData_1==1){
             
              if(V!=0&&Tf!=0){
              theta=asin((g*Tf)/(2*V));
              t=Tf;
              enoughData++;
            }
             else if(theta!=0&&Tf!=0){
              theta=asin((g*Tf)/(2*V));
              t=Tf;
              enoughData++;
            }
             else if(t!=0&&Tf!=0){
               
            }
            
            
           else if(V!=0&&H!=0){
              theta=asin((2*g*sqrt(H))/V);
              t=(2*V*sin(theta))/g;
              enoughData++;
            }
            else if(theta!=0&&H!=0){
              V=sqrt(2*g*H)/sin(theta);
              t=(2*V*sin(theta))/g;
              enoughData++;
            }
            else if(t!=0&&H!=0){
               
            }


            else if(V!=0&&R!=0){
              theta=0.5*((asin((g*R)/(V*V))));
              t=(2*V*sin(theta))/g;
              enoughData++;
            }
            else if(theta!=0&&R!=0){
              theta=asin((g*T)/(2*V));
              V=sqrt((g*R)/sin(2*theta));
              enoughData++;
            }
            else if(t!=0&&R!=0){
                
            }
              
           }
          else if(enoughData==0&&enoughData_1>1){
            if(Tf!=0&&H!=0){
              t=2*sqrt(H);
              Vy=(g*Tf)/2;
              cout<<"Is the value of Horizontal component of V (Vx) known?(if not enter 0)";
              cin>>Vx;
              if(Vx!=0){
               theta=atan(Vy/Vx);                                                         //value of Vx couldnt be zero
               V=Vx/cos(theta);
              }

            }
            else if(T!=0&&R!=0){
              Vx=(2*R/Tf);
              cout<<"Is the value of Vertical component of V (Vy) known?(if not enter 0)";
              cin>>Vy;
              theta=atan(Vy/Vx);                                       //value of Vy could be zero too
              V=Vx/cos(theta);
            }
            else if(H!=0&&R!=0){
              theta=atan((4*H)/R);
              V=sqrt((g*R)/sin(2*theta));
              t=(2*V*sin(theta)/g);
            }

          }
           
        }

        return 0;
   
}

int withFriction(double x,double y,double V,double t,double theta,double D){
    double tempX=x,tempY=y,tempV=V,Ax=0,Ay=0;
    double Tm=t/100000,tx=0;    
    double tmax=t*100000;

    double Vx=V*cos(theta);
    double  Vy=V*sin(theta);

   
    std::vector<double> positionX,positionY,time,pE,kE,tE,HA,VA,VC,HC;
    
    cout<<"                                     OUTPUT DISPLAY                             "<<endl<<endl<<endl;
    for(int i=0;i<tmax;i++) {    
          
    if(tx>t) break;   
    if(tempY<0) break; 

    Ax=-(D/m)*tempV*Vx;
    Ay=-g-(D/m)*tempV*Vy;
    Vx=Vx+Ax*Tm;
    Vy=Vy+Ay*Tm;
    tempV=sqrt((Vx*Vx)+(Vy*Vy));
    tempX=tempX+Vx*Tm+0.5*Ax*(Tm*Tm);
    tempY=tempY+Vy*Tm+0.5*Ay*(Tm*Tm);
    tx=tx+Tm;
    
    }
    cout<<endl;
  
    cout<<"Horizontal Component of Velocity at "<<tx<< " second is "<<Vx<<endl;
    cout<<"Vertical Component of Velocity at "<<tx<< " second is "<<Vy<<endl;
    cout<<"Horizontal Component of Accelaration at "<<tx<< " second is "<<Ax<<endl;
    cout<<"Vertical Component of Accelaration at "<<tx<< " second is "<<Ay<<endl;
    cout<<"Distance from origin: "<<distanceFromOrigin(tempX,tempY)<<endl;
    cout<<endl;
  
  
    tx=0,tempX=x,tempY=y,tempV=V,Ax=0,Ay=0;
    Vx=V*cos(theta);
    Vy=V*sin(theta);

    cout<<"Time"<<"    x "<<"      y "<<"       V "<<"       Dis0 "<<"     kE "<<"      pE "<<"      tE"<<endl;
    for(int i=0;i<((2*V*sin(theta)/g)*100000);i++) {  

        if(tempY<0) break; 

        Ax=-(D/m)*tempV*Vx;
        Ay=-g-(D/m)*tempV*Vy;
      
      if(i%25000==0){
        cout<<FIXED_FLOAT(tx)<<"    ";
        cout<<FIXED_FLOAT(tempX)<<"    ";
        cout<<FIXED_FLOAT(tempY)<<"    ";
        cout<<FIXED_FLOAT(tempV)<<"    ";
        distanceFromOrigin(tempX,tempY);
        preservationofEnergy(tempV,tempY);
        cout<<endl;
      }
       positionX.push_back(tempX);
       positionY.push_back(tempY);
       time.push_back(tx);
       kE.push_back(0.5*m*(tempV*tempV));
       pE.push_back(m*g*tempY);
       tE.push_back((0.5*m*(tempV*tempV))+(m*g*tempY));
       HA.push_back(Ax);
       VA.push_back(Ay);
       VC.push_back(Vx);
       HC.push_back(Vy);
    
       Vx=Vx+Ax*Tm;
       Vy=Vy+Ay*Tm;
       tempV=sqrt((Vx*Vx)+(Vy*Vy));
       tempX=tempX+Vx*Tm+0.5*Ax*(Tm*Tm);
       tempY=tempY+Vy*Tm+0.5*Ay*(Tm*Tm);
       tx=tx+Tm;
       if(i%5000==0){
       plt::clf();
       plt::named_plot("y=xtan(theta)-(gx^2)/2u^2cos(theta)",positionX,positionY);
       plt::legend();
       plt::pause(0.001);
       }
    }
  
    
    printValues(tx,Vx,Vy,Ax,Ay,tempX,tempY,tempV);   

    EnergyGraph();
    plt::clf(); 
    plt::named_plot("y=xtan(theta)-(gx^2)/2u^2cos(theta)",positionX,positionY);
    plt::legend();
    plt::save("Trajectory.png");
    plt::show();
    plt::clf(); 
    plt::named_plot("Total Energy",time,tE);
    plt::named_plot("Potential Energy",time,pE);
    plt::named_plot("Kinetic Energy",time,kE);
    plt::legend();     
    plt::save("PreservationOfEnergy.png");
    plt::show();
    plt::clf(); 
    plt::named_plot("Horizontal Component of Velocity",time,HC);
    plt::named_plot("Vertical Component of Velocity",time,VC);
    plt::legend();
    plt::save("ComponentsOfVelocity.png");    
    plt::show();
    plt::clf();
    plt::named_plot("Horizontal Component of Accelaration",time,HA);
    plt::named_plot("Vertical Component of Accelaration",time,VA);
    plt::legend();
    plt::save("ComponentsOfAccelaration.png");
    plt::show();

    return 0;
}
    
int withoutFriction(double x,double y,double V,double t,double theta){
    double tempX=x,tempY=y,tempV=V,Ax=0,Ay=0;
    double Tm=t/100000,tx=0;    
    double tmax=t*100000;
    D=0;
    double Vx=V*cos(theta);
    double  Vy=V*sin(theta);

   
    std::vector<double> positionX,positionY,time,pE,kE,tE,HC,VC,HA,VA;
    cout<<"                                     OUTPUT DISPLAY                             "<<endl<<endl<<endl;
    
    for(int i=0;i<tmax;i++) {    
          
    if(tx>t) break;   
    if(tempY<0) break; 

    Ax=-(D/m)*tempV*Vx;
    Ay=-g-(D/m)*tempV*Vy;
    Vx=Vx+Ax*Tm;
    Vy=Vy+Ay*Tm;
    tempV=sqrt((Vx*Vx)+(Vy*Vy));
    tempX=tempX+Vx*Tm+0.5*Ax*(Tm*Tm);
    tempY=tempY+Vy*Tm+0.5*Ay*(Tm*Tm);
    tx=tx+Tm;
    
    }
    cout<<endl;
    cout<<"Horizontal Component of Velocity at "<<tx<< " second is "<<Vx<<endl;
    cout<<"Vertical Component of Velocity at "<<tx<< " second is "<<Vy<<endl;
    cout<<"Horizontal Component of Accelaration at "<<tx<< " second is "<<Ax<<endl;
    cout<<"Vertical Component of Accelaration at "<<tx<< " second is "<<Ay<<endl;
    cout<<"Distance from origin: "<<distanceFromOrigin(tempX,tempY)<<endl;
    cout<<endl;
  
  
    tx=0,tempX=x,tempY=y,tempV=V,Ax=0,Ay=0;
    Vx=V*cos(theta);
    Vy=V*sin(theta);
    
    double x0 = -20; // Initial values assumed
    newtonRaphson(x0);
  cout<<"Time"<<"    x "<<"      y "<<"       V "<<"       Dis0 "<<"     kE "<<"      pE "<<"      tE"<<endl;
  for(int i=0;i<(((2*V*sin(theta)/g)*100000));i++) {  

        if(tempY<0) break; 

        Ax=-(D/m)*tempV*Vx;
        Ay=-g-(D/m)*tempV*Vy;
      
       if(i%25000==0){
        cout<<tx<<"    ";
        cout<<tempX<<"    ";
        cout<<tempY<<"    ";
        cout<<tempV<<"    ";
        distanceFromOrigin(tempX,tempY);
        preservationofEnergy(tempV,tempY);
        cout<<endl;
      }
       positionX.push_back(tempX);
       positionY.push_back(tempY);
       time.push_back(tx);
       kE.push_back(0.5*m*(tempV*tempV));
       pE.push_back(m*g*tempY);
       tE.push_back((0.5*m*(tempV*tempV))+(m*g*tempY));
    
       Vx=Vx+Ax*Tm;
       Vy=Vy+Ay*Tm;
       tempV=sqrt((Vx*Vx)+(Vy*Vy));
       tempX=tempX+Vx*Tm+0.5*Ax*(Tm*Tm);
       tempY=tempY+Vy*Tm+0.5*Ay*(Tm*Tm);
       tx=tx+Tm;
       if(i%500==0){
       plt::clf();
       plt::named_plot("y=xtan(theta)-(gx^2)/2u^2cos(theta)",positionX,positionY);
       plt::named_plot("Total Energy",time,tE);
       plt::named_plot("Potential Energy",time,pE);
       plt::named_plot("Kinetic Energy",time,kE);
       plt::legend();
       plt::pause(0.01);
       }
    }
  
    printValues(tx,Vx,Vy,Ax,Ay,tempX,tempY,tempV);
    
    EnergyGraph();
    
    plt::clf(); 
    plt::named_plot("y=xtan(theta)-(gx^2)/2u^2cos(theta)",positionX,positionY);
    plt::legend();
    plt::save("Trajectory.png");
    plt::show();
    plt::clf(); 
    plt::named_plot("Total Energy",time,tE);
    plt::named_plot("Potential Energy",time,pE);
    plt::named_plot("Kinetic Energy",time,kE);
    plt::legend();     
    plt::save("PreservationOfEnergyWithoutFriction.png");
    plt::show();
    plt::clf(); 
    plt::named_plot("Horizontal Component of Velocity",time,HC);
    plt::named_plot("Vertical Component of Velocity",time,VC);
    plt::legend();
    plt::save("ComponentsOfVelocity.png");
    plt::show();
    plt::clf();
    plt::named_plot("Horizontal Component of Accelaration",time,HA);
    plt::named_plot("Vertical Component of Accelaration",time,VA);
    plt::legend();
    plt::save("ComponentsOfVelocity.png");
    plt::show();

    return 0;
}

int printValues(double txdash,double Vx,double Vy,double Ax,double Ay,double xdash,double ydash,double Vdash){
      cout<<endl;
      cout<<"Time of flight(T): "<<txdash<<endl;          
      cout<<"Vertical Range(Maximum Distance Form Origin): "<<xdash<<endl;
      cout<<"Velocity when the object hits the ground: "<<Vdash <<endl;
      lenthofthepath( V,theta,t);
      cout<<endl;
      return 0;
  }

int distanceFromOrigin(double xdash,double ydash){
       distanceFromO=sqrt((xdash*xdash)+(ydash*ydash));
      
       cout<<FIXED_FLOAT(distanceFromO)<<"    ";
}

int preservationofEnergy(double V,double y){
     KineticEnergy=0.5*m*(V*V);
     PotentialEnergy=m*g*y;
     TotalEnergy=KineticEnergy+PotentialEnergy;
     cout<<FIXED_FLOAT(KineticEnergy)<<"    ";
     cout<<FIXED_FLOAT(PotentialEnergy)<<"    ";
     cout<<FIXED_FLOAT(TotalEnergy)<<"    ";
     return 0;
}

void EnergyGraph(){
    double xdash=x;
    double ydash=y;
    double vdash=V;
    double Ax=0,Ay=0;
    double Tm=t/100000,tx=0;    
    double tmax=t*100000;

    double Vx=V*cos(theta);
    double Vy=V*sin(theta);

   
    std::vector<double> time,pE,kE,tE;
    
    KineticEnergy=0;
    PotentialEnergy=0;
    TotalEnergy=0;

    for(int i=0;i<tmax;i++) {    
          
  
    if(ydash<0) break; 

    Ax=-(D/m)*vdash*Vx;
    Ay=-g-(D/m)*vdash*Vy;
    Vx=Vx+Ax*Tm;
    Vy=Vy+Ay*Tm;
    vdash=sqrt((Vx*Vx)+(Vy*Vy));
    xdash=xdash+Vx*Tm+0.5*Ax*(Tm*Tm);
    ydash=ydash+Vy*Tm+0.5*Ay*(Tm*Tm);
    tx=tx+Tm;
    time.push_back(tx);
    kE.push_back(0.5*m*(vdash*vdash));
    pE.push_back(m*g*ydash);
    tE.push_back(0.5*m*(vdash*vdash)+m*g*ydash);
    
    if(i%500==0)
            {
       plt::clf();
       plt::named_plot("Total Energy",time,tE);
       plt::named_plot("Potential Energy",time,pE);
       plt::named_plot("Kinetic Energy",time,kE);
       plt::legend();
       plt::pause(0.02);
       }
    
  }
}

void newtonRaphson(double x)
{
    double h = func(x) / derivFunc(x);
    while (abs(h) >= EPSILON)
    {
        h = func(x) / derivFunc(x);

        x = x - h;
    }

    cout << "The value of the root is : " << x;
}

double func(double x)
{
	return x*tan(theta)-(g*x*x)/(2*u*u*cos(theta)*cos(theta));
}

int derivFunc(double x){
  return (tan(theta)-(2*g*x)/(2*u*u*cos(theta)));
}

double rootFunc( double tdash)
{
   return (sqrt((-g*tdash+V*sin(theta))*(-g*tdash+V*sin(theta))+(V*cos(theta))*(V*cos(theta))));
}

double calculate(double lower_limit, double upper_limit,int interval_limit )
{
    double value;
    double interval_size = (upper_limit - lower_limit)
                          / interval_limit;
    double sum = rootFunc(lower_limit) + rootFunc(upper_limit);
 

    for (int i = 1 ; i < interval_limit ; i++)
    {
        if (i % 3 == 0)
            sum = sum + 2 * rootFunc(lower_limit + i * interval_size);
        else
            sum = sum + 3 * rootFunc(lower_limit + i * interval_size);
    }
    return ( 3 * interval_size / 8 ) * sum ;
}

string diffTerm(string pTerm)
{

	string coeffStr = "", S = "";
	int i;

	for (i = 0; pTerm[i] != 't'; i++)
		coeffStr.push_back(pTerm[i]);

	long long coeff
		= atol(coeffStr.c_str());


	string powStr = "";
	for (i = i + 2; i != pTerm.size(); i++)
		powStr.push_back(pTerm[i]);

	long long power
		= atol(powStr.c_str());
	string a, b;

	Converting the value
	to the string
	ostringstream str1, str2;

	For ax^n, we find (n)*a*x^(n-1)
	coeff = coeff * power;
	str1 << coeff;
	a = str1.str();
	power--;
	str2 << power;
	b = str2.str();
	 S += a + "t^" + b;

	return S;
  return 0;
}

string diffstr(string& poly)
{


	istringstream is(poly);

	string pTerm, S = "";

	while (is >> pTerm) {


		if (pTerm == "+") {
			S += " + ";
			continue;
		}

		if (pTerm == "-") {
			S += " - ";
			continue;
		}

		else
			S += diffTerm(pTerm);
	}
	return S;

  return 0;
}

string diffofY(double V,double theta){
    double Vx = V * cos(theta);
    double Vy=V*sin(theta);
     
     string x=to_string(Vx);
     string y=to_string(Vy);
     string z=to_string(0.5*g);

     string str=z+"t^2+"+y+"t";

     return diffstr(str);

}

string diffofX(double V,double theta){

     double Vx = V * cos(theta);
     double Vy=V*sin(theta);

     string x=to_string(Vx);
     string y=to_string(Vy);
     string z=to_string(0.5*g);

     string str=x+"t";
     return diffstr(str);
   }

double lenthofthepath(double V,double theta,double t){
     
    string DfuncofX= diffofX(V,theta);
    string DfuncofY= diffofY(V,theta);

    int interval_limit = 10;
    double lower_limit = 1;
    double upper_limit = 10;
    double integral_res = calculate(lower_limit, upper_limit,
                                   interval_limit);
    cout << "length of the Arc of the projectile motion is:"<<integral_res;
    return 0;
}
