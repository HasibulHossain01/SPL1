#include<bits/stdc++.h>
using namespace std;

int inputFunctionK();
int inputFuntionU(int enoughData,float V,float t,float theta);
int printValues(float tx,float Vx,float Vy,float Ax,float Ay,float x,float y,float V);
int withFriction(float x,float y,float V,float Vx,float Vy,float m,float t,float theta,float D);
//int withoutFriction(float x,float y,float V,float t,float theta);
int distanceFromOrigin(float x,float y);
int preservationofEnergy(float V,float y);

float D,Ax,Ay,V,Vx,Vy,p,C,A,r,distanceO;
float t,T,tx,tmax,Tf,H,R;
float x,y,theta,m;
const float g=9.81,pi=3.1416;
int enoughData=0;
char chF;

int main(void){

    cout<<" Is it a projectile motion with air resistance?(y/n): ";
    cin>>chF;
    inputFunctionK();
   
    if((enoughData==3)&&(chF='y')) withFriction(x,y,V,Vx,Vy,m,t,theta,D);
    
    //else withoutFriction(x,y,V,t,theta);   

     
     return 0;

}

int inputFunctionK(){
    
    cout<<"Enter initial Co-ordinate of x(mandatory): ";
    cin>>x;
    cout<<endl;

    cout<<"Enter initial Co-ordinate of y(mandatory): ";
    cin>>y;
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
     if(theta>0) enoughData++; //value of theta can be 0
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
         A=pi*(r*r);
         D=(p*C*A)/2;
     }
    }
     cout<<endl;
     
     if(enoughData<3) {
        inputFuntionU(enoughData,V,t,theta);
     }

     theta=theta*3.1416/180;
     Vx=V*cos(theta);
     Vy=V*sin(theta);

     return 0;
}

int inputFuntionU(int enoughData,float V,float t,float theta){
        int enoughData_1=enoughData;
        cout<<"Enter any of the three of Time of flight(T), MAximum Vertical Renge(H) or Maximum horizontal Renge(R)[if unknown enter 0]: "<<endl;
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
                //need more work

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
                //need more work
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
                //need more work
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

int withFriction(float x,float y,float V,float Vx,float Vy,float m,float t,float theta,float D){
   
    T=t/100000,tx=0;    
    tmax=t*100000;
   
    for(int i=0;i<tmax;i++) {    
          
    if(tx>t) break;   
    if(y<0) break; 

    Ax=-(D/m)*V*Vx;
    Ay=-g-(D/m)*V*Vy;
    Vx=Vx+Ax*T;
    Vy=Vy+Ay*T;
    V=sqrt((Vx*Vx)+(Vy*Vy));
    x=x+Vx*T+0.5*Ax*(T*T);
    y=y+Vy*T+0.5*Ay*(T*T);
    tx=tx+T;
    }
    cout<<endl;
    cout<<"Horizontal Component of Velocity at "<<tx<< " second is "<<Vx<<endl;
    cout<<"Vertical Component of Velocity at "<<tx<< " second is "<<Vy<<endl;
    cout<<"Horizontal Component of Accelaration at "<<tx<< " second is "<<Ax<<endl;
    cout<<"Vertical Component of Accelaration at "<<tx<< " second is "<<Ay<<endl;
    cout<<"Distance from origin: "<<distanceO<<endl;
    cout<<endl;

    
  for(int i=0;i<(20*100000);i++) {    
    if(y<0) break; 

    Ax=-(D/m)*V*Vx;
    Ay=-g-(D/m)*V*Vy;
    
    if(i==0||i%10000==0) {
 
    cout<<tx<<"    ";
    cout<<x<<"    ";
    cout<<y<<"    ";
    cout<<V<<"    ";
    distanceO=distanceFromOrigin(x,y);
    preservationofEnergy(V,y);
    cout<<endl;
    
    }
    Vx=Vx+Ax*T;
    Vy=Vy+Ay*T;
    V=sqrt((Vx*Vx)+(Vy*Vy));
    x=x+Vx*T+0.5*Ax*(T*T);
    y=y+Vy*T+0.5*Ay*(T*T);
    tx=tx+T;

    }

    printValues(tx,Vx,Vy,Ax,Ay,x,y,V);
    
    return 0;
}

//int withoutFriction(float x,float y,float V,float t,float theta){
    // return 0;   
//}

int printValues(float tx,float Vx,float Vy,float Ax,float Ay,float x,float y,float V){
      cout<<"Time of flight(T): "<<tx<<endl;          
      cout<<"Vertical Renge(Maximum Distance Form Origin): "<<x<<endl;
      cout<<"Velocity when the object hits the ground: "<<V <<endl;
      cout<<endl;
      return 0;
  }

int distanceFromOrigin(float x,float y){
      float distance=sqrt((x*x)+(y*y));
      cout<<distance<<"    ";
      return distance;
}

int preservationofEnergy(float V,float y){
     float KineticEnergy=0.5*m*(V*V);
     float PotentialEnergy=m*g*y;
     float TotalEnergy=KineticEnergy+PotentialEnergy;
     cout<<KineticEnergy<<"    ";
     cout<<PotentialEnergy<<"    ";
     cout<<TotalEnergy<<"    ";
      return 0;
}