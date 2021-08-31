#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <iomanip>

using namespace std;

#define n 100     ///the number of section
#define nghost 2  ///extra section
#define small 0.0000001

void orderChoice(double order,double u[],double ul[],double ur[],int m)

{
    if (m==0)
    {
        if(order==1) cout<<"1st order"<<endl;
        else if (order>2 & order<3) 
        {
            cout<<"2nd order"<<endl;
            if (order==2.1) cout<<" Beam Warning method"<<endl;
            else if (order==2.2) cout<<" Lax Wendrof method"<<endl;
            else if (order==2.3) cout<<" Foumn method"<<endl;
        }
    }
    
    for(int i=nghost; i<=n+nghost+1; i++)
    {
        if (order==1)
        {
            ul[i]=u[i-1];
            ur[i]=u[i];
        }

        /////Beam Warning method/////
        else if (order==2.1)
        {
            ul[i]=u[i-1]+0.5*(u[i-1]-u[i-2]);
            ur[i]=u[i]-0.5*(u[i+1]-u[i]);

            ///time pass
            //ul[i]=ul[i]-c*dt*0.5*(u[i]-u[i-1])/dx;
            //ur[i]=ur[i]-c*dt*0.5*(u[i+2]-u[i-1])/dx;
        }

        /////Lax Wendrof method/////
        else if (order==2.2)
        {
            ul[i]=u[i-1]+0.5*(u[i]-u[i-1]);
            ur[i]=u[i]-0.5*(u[i]-u[i-1]);
        }

        /////Foumn method/////
        else if (order==2.3)
        {
            ul[i]=u[i-1]+0.25*(u[i]-u[i-2]);
            ur[i]=u[i]+0.25*(u[i+1]-u[i-1]);
        }
    }
}

void SoundofSpeedCal(double a_L[],double a_R[],double Pre_L[],double Pre_R[],double Rho_L[],double Rho_R[],double gamma)
{
    for(int i=nghost; i<=n+nghost+1; i++)
    {
        a_L[i]=sqrt(gamma*Pre_L[i]/Rho_L[i]);
        a_R[i]=sqrt(gamma*Pre_R[i]/Rho_R[i]);
    }
}

void calcflux(double Rho_L[],double Rho_R[],double Vel_L[],double Vel_R[],double Pre_L[],double Pre_R[],double Rho_E_L[],double Rho_E_R[],double a_L[],double a_R[],double F1[],double F2[],double F3[],int m,double dx,double dt)
{
    double *F1_L=new double[n+nghost*2+1];
    double *F2_L=new double[n+nghost*2+1];
    double *F3_L=new double[n+nghost*2+1];
    double *F1_R=new double[n+nghost*2+1];
    double *F2_R=new double[n+nghost*2+1];
    double *F3_R=new double[n+nghost*2+1];

    for(int i=nghost; i<=n+nghost+1; i++)
    {
        F1_L[i]=Rho_L[i]*Vel_L[i];                  //[F-1]//
        F2_L[i]=Rho_L[i]*pow(Vel_L[i],2)+Pre_L[i];  //[F-2]//
        F3_L[i]=Vel_L[i]*(Rho_E_L[i]+Pre_L[i]);     //[F-3]//

        F1_R[i]=Rho_R[i]*Vel_R[i];                  //[F-1]//
        F2_R[i]=Rho_R[i]*pow(Vel_R[i],2)+Pre_R[i];  //[F-2]//
        F3_R[i]=Vel_R[i]*(Rho_E_R[i]+Pre_R[i]);     

        /*if ( Riemann_type == 0) 
        {
        F1[i]=0.5*(F1_L[i]+F1_R[i]) - 0.5*(max(a_L[i],a_R[i]) + max(abs(Vel_L[i]), abs(Vel_R[i])))*(Rho_R[i] - Rho_L[i]); //Lax-friendrichs
        F2[i]=0.5*(F2_L[i]+F2_R[i]) - 0.5*(max(a_L[i],a_R[i]) + max(abs(Vel_L[i]), abs(Vel_R[i])))*(Rho_R[i]*Vel_R[i] - Rho_L[i]*Vel_L[i]);
        F3[i]=0.5*(F3_L[i]+F3_R[i]) - 0.5*(max(a_L[i],a_R[i]) + max(abs(Vel_L[i]), abs(Vel_R[i])))*(Rho_E_R[i] - Rho_E_L[i]);
        //std::cout <<"F1["<<i<<"]="<< F1[i] << std::endl;
        }
        else if (Riemann_type == 1)
        {*/
        F1[i]=0.5*(F1_L[i]+F1_R[i]) - dx/dt/2.0*(Rho_R[i] - Rho_L[i]); //Lax-friendrichs
        F2[i]=0.5*(F2_L[i]+F2_R[i]) - dx/dt/2.0*(Rho_R[i]*Vel_R[i] - Rho_L[i]*Vel_L[i]);
        F3[i]=0.5*(F3_L[i]+F3_R[i]) - dx/dt/2.0*(Rho_E_R[i] - Rho_E_L[i]);
        //std::cout <<"F1["<<i<<"]="<< F1[i] << std::endl;
        //}
    }
}

void Equation(double u[],double f[],double dx,double dt)
{
    for(int i=nghost; i<=n+nghost; i++)
    {
        u[i]=u[i]-dt/dx*(f[i+1]-f[i]); 
    }
}

void EdgeBoundary(double u[])
{
    u[ nghost - nghost ] = u[ nghost - nghost + 1];
    u[ nghost - nghost + 1] = u[ nghost - nghost + 2 ];
    u[ n + nghost*2 - 1] = u[ n + nghost*2 - 2 ];
    u[ n + nghost*2] = u[ n + nghost*2 - 1 ];
}



main()
{
    int i;
    int m=0;
    double gamma=1.4; //heat ratio(比熱比)

    /////time
    double tSta=0.0;
    double tEnd=0.03;

    //spatial domain
    double xmax=1.0;
    double xmin=0.0;

    /* 
    1st order
        1.0 : 
     2nd order 
        2.1 : Beam Warning method
        2.2 : Lax Wendrof method
        2.3 : Foumn method
    */
    double order=1.0;       

    double dx=(xmax-xmin)/n;
    double dt=0.00001;

    double *x=new double[n+nghost*2];   //grid

    //Rho(質量)
    double *Rho=new double[n+nghost*2];
    double *Rho_L=new double[n+nghost*2+1];
    double *Rho_R=new double[n+nghost*2+1];

    //Pressure(圧力)
    double *Pre=new double[n+nghost*2];
    double *Pre_L=new double[n+nghost*2+1];
    double *Pre_R=new double[n+nghost*2+1];

    //Velocity(速度)
    double *Vel=new double[n+nghost*2];
    double *Vel_L=new double[n+nghost*2+1];
    double *Vel_R=new double[n+nghost*2+1];

    //speed of sound (音速) =sqrt(gamma * p / rho)
    double *a_L=new double[n+nghost*2+1];
    double *a_R=new double[n+nghost*2+1];

    double *Rho_Vel=new double[n+nghost*2];      //Rho x Velocity
    double *Rho_Vel2_Pre=new double[n+nghost*2]; //Rho x Velocity^2+Pressure
    
    double *Rho_E=new double[n+nghost*2];        //Rho x Internal Energy
    double *Rho_E_L=new double[n+nghost*2+1];
    double *Rho_E_R=new double[n+nghost*2+1];

    double *F3=new double[n+nghost*2];

    /////initial condition(t=0)////
    x[0]=0.0 - dx * nghost;
    for(i=0; i<=n+nghost*2;i++)
    {
        x[i+1]=x[i]+dx;
    }

    for(i=0; i<=n+nghost*2; i++)
    {
        ///high position
        if(x[i]<=0.5)
        {
            Pre[i]=1.0;
            Vel[i]=0.0;
            Rho[i]=1.0;                                           //[U-1]//
        }
        ///low position
        else
        {
            Pre[i]=0.0;
            Vel[i]=0.0;
            Rho[i]=0.125;                                         //[U-1]//
        }

        Rho_Vel[i]=Rho[i]*Vel[i];                                 //[U-2]//
        Rho_E[i]=0.5*Rho[i]*pow(Vel[i],2)+Pre[i]/(gamma-1);       //[U-3]//

    }

    /////calculation (t>0)/////
    while(tSta<tEnd)
    {
        orderChoice(order,Rho,Rho_L,Rho_R,m);
        orderChoice(order,Vel,Vel_L,Vel_R,m);
        orderChoice(order,Pre,Pre_L,Pre_R,m);
        orderChoice(order,Rho_E,Rho_E_L,Rho_E_R,m);
        
        SoundofSpeedCal(a_L,a_R,Pre_L,Pre_R,Rho_L,Rho_R,gamma);
        calcflux(Rho_L,Rho_R,Vel_L,Vel_R,Pre_L,Pre_R,Rho_E_L,Rho_E_R,a_L,a_R,Rho_Vel,Rho_Vel2_Pre,F3,m,dx,dt);
        
        Equation(Rho,Rho_Vel,dx,dt);
        Equation(Rho_Vel,Rho_Vel2_Pre,dx,dt);
        Equation(Rho_E,F3,dx,dt);

        for(int i=nghost; i<=n+nghost; i++)
        {
            Vel[i]=Rho_Vel[i]/max(Rho[i], small);
            Pre[i]=(gamma-1)*(Rho_E[i] - 0.5*Rho[i]*pow(Vel[i],2));
        }

        EdgeBoundary(Rho);
        EdgeBoundary(Vel);
        EdgeBoundary(Pre);
        EdgeBoundary(Rho_E);
        EdgeBoundary(Rho_Vel);

        tSta=tSta+dt;
        m++;
        for (int i = 0; i < n; i++)
        {

        
        std::cout <<"Rho["<<i<<"]="<< Rho[i] << std::endl;
        }
    }
    
}

//double *e=new double[n+nghost*2];   //internal(内部エネルギー)