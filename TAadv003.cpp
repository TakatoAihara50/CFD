#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <iomanip>

using namespace std;

//num of grid point
#define n 100
#define nghost 2


void orderChoice(double order,double c,double dx,double dt,double u[],double ul[],double ur[],int m)
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

            ///time pass
            //ul[i]=ul[i]-c*dt*0.5*(u[i]-u[i-1])/dx;
            //ur[i]=ur[i]-c*dt*0.5*(u[i+2]-u[i-1])/dx;
        }

        /////Foumn method/////
        else if (order==2.3)
        {
            ul[i]=u[i-1]+0.25*(u[i]-u[i-2]);
            ur[i]=u[i]+0.25*(u[i+1]-u[i-1]);

            ///time pass
            //ul[i]=ul[i]-c*dt*0.5*(u[i]-u[i-1])/dx;
            //ur[i]=ur[i]-c*dt*0.5*(u[i+2]-u[i-1])/dx;
        }
    }
}



void calcflux(int riemann_type,double c,double dt,double dx,double f[],double ul[],double ur[],double fluxl[],double fluxr[],int m)
{
    if (m==0)
    {
        if(riemann_type==1) cout<<"rho riemann solver"<<endl;
        else if (riemann_type==2) cout<<"Lax-Friedrich Riemann solver"<<endl;
    }

    for(int i=nghost; i<=n+nghost+1; i++)
    {
        /////rho riemann solver/////
        if(riemann_type==1)
        {
            if(c>=0) f[i]=c*ul[i];
            else f[i]=c*ur[i];
        }

        /////Lax-Friedrich Riemann solver/////
        else if(riemann_type==2)
        {
             for(int i=nghost; i<=n+nghost+1; i++)
            {
                fluxl[i]=c*ul[i];
                fluxr[i]=c*ur[i];
            }
            f[i]=0.5*(fluxl[i]+fluxr[i]) - 0.5*dx/dt*(ur[i] - ul[i]);
        }
    }
}



int main()
{
    int i;
    int m=0;

    /////choices
    double order;       
    /*
     1st order
        1.0 : 
     2nd order 
        2.1 : Beam Warning method
        2.2 : Lax Wendrof method
        2.3 : Foumn method
    */
    int riemann_type=1;   // 1 : rho      , 2 : lax-friedrich

    /////time
    double tSta=0.0;
    double tEnd=1.0;

    /////wave speed
    double c = 1;

    /////spatial domain
    double xmax=1;
    double xmin=0;

    double dx=(xmax-xmin)/n;
    double dt=0.1*dx/abs(c); //dt=0.001

    double *x=new double[n+nghost*2];
    double *u=new double[n+nghost*2];
    double *f=new double[n+nghost*2+1];
    double *ul=new double[n+nghost*2+1];
    double *ur=new double[n+nghost*2+1];
    double *fluxl=new double[n+nghost*2+1];
    double *fluxr=new double[n+nghost*2+1];

    ofstream fk;
    fk.open("test003.txt");


    /////initial condition (t=0)/////
    x[0]=0.0 - dx * nghost;
    for(i=0; i<=n+nghost*2;i++)
    {
        x[i+1]=x[i]+dx;
    }
    for(i=0; i<=n+nghost*2; i++)
    {
        u[i]=exp(-0.5*pow(((x[i]-0.5)/0.08),2));
        // sin curve
        //u[i]= (1.1+sin(4.0*3.14*(x[i] -x[2])));
    }


    /////clculation (t>0) /////
    cout<<"Chose clculation order."<<endl;
    cout<<"1st order\n 1.0 :\n2nd order\n 2.1 : Beam Warning method\n 2.2 : Lax Wendrof method\n 2.3 : Foumn method"<<endl;
    cin>>order;

    while(tSta<tEnd)
    {
        orderChoice(order, c, dx, dt, u, ul, ur, m);
        calcflux(riemann_type, c, dt, dx, f, ul, ur, fluxl, fluxr, m);
        
        for(int i=nghost; i<=n+nghost; i++)
        {
            u[i]=u[i]-dt/dx*(f[i+1]-f[i]);
            fk<<u[i]<<" ";
        }
        fk << endl;

        /////respawn cycle
        u[ nghost - nghost ] = u[ n + nghost - 1];
        u[ nghost - nghost + 1] = u[ n + nghost ];
        u[ n + nghost*2 - 1] = u[ nghost ];
        u[ n + nghost*2] = u[ nghost + 1];

        tSta=tSta+dt;
        m++;
    }
    fk.close();
    
    return 0;
}





//cout<<"x["<<i<<"]="<<x[i]<<endl;
//cout<<u[i]<<endl;

/////time end output/////
//    for(int i=nghost; i<=n+nghost; i++)
//    {
//        std::cout << u[i] << std::endl;
//    }