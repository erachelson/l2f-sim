#include <L2Fsim/flight_zone/thermal/std_thermal.hpp>

using namespace std;

namespace L2Fsim {
    
/*---------------------------------------------------------------------------------------
 ---------------------------       CONSTRUCTOR / DESTRUCTOR     -------------------------
 --------------------------------------------------------------------------------------*/


std_thermal::std_thermal(int mod,double tB,double XC0,double YC0,double ZC0,
                         double Zi,double wstar, int Lifetime, double Ksi,bool read):
                         tBirth(tB),xc0(XC0),yc0(YC0),zc0(ZC0)
{
    if (!read)
    {
        //All time parameters should be used in same units of time. [s]
        //||--------------------------------------------------||
        //|| var           ||      min      ||     max        ||
        //||--------------------------------------------------||
        //|| Tlife [s]     ||      600      ||     1200       ||
        //|| wind [m/s]    ||     -0.2      ||     0.2        ||
        //|| ksi           ||      0.1      ||     0.35       ||
        //||--------------------------------------------------||
        
        // DATA Peak velocity and Mixing layer (cf 'version3')
        int size_DATA=13;
        double DATA[size_DATA][3];
        DATA[0][0]=1.14;
        DATA[1][0]=1.48;
        DATA[2][0]=1.64;
        DATA[3][0]=1.97;
        DATA[4][0]=2.53;
        DATA[5][0]=2.38;
        DATA[6][0]=2.56;
        DATA[7][0]=2.69;
        DATA[8][0]=2.44;
        DATA[9][0]=2.25;
        DATA[10][0]=1.79;
        DATA[11][0]=1.31;
        DATA[12][0]=1.26;
        
        // The lifetime inversely proportional to the Peak velocity
        int minLifeTime=600; int maxLifeTime=1200;
        double wmin=DATA[0][0];double wmax=DATA[0][0];
        for (int i=0;i<size_DATA;i++)
        {
            if (DATA[i][0]<wmin) wmin = DATA[i][0];
            if (DATA[i][0]>wmax) wmax = DATA[i][0];
        }
        double a=minLifeTime;double b=maxLifeTime;
        double c=1./wmax;double d=1./wmin;
        for (int i=0;i<size_DATA;i++)
        {
            DATA[i][1]=floor(a+(b-a)*((1./DATA[i][0]-c)/(d-c)));
        }
        
        //Thermal initialized with all default values
        double coef = 1.7; // boost w_peak
        int numberconfig=rand_int(0,13);
        model=mod;
        w_star   = coef * DATA[numberconfig][0];
        zi       = Zi; //+ 30*normalLaw();
        lifeTime = DATA[numberconfig][1];
        ksi = rand_double(0.1,0.35);
    }
    else
    {
        model    = mod;
        w_star   = wstar;
        zi       = Zi;
        lifeTime = Lifetime;
        ksi      = Ksi;
    }
}
    
// destructor
    // no pointer in attributes
    // nothing to delete
    

/*---------------------------------------------------------------------------------------
 ---------------------------              GETTEURS              -------------------------
 --------------------------------------------------------------------------------------*/

// getCenter
vector<double> std_thermal::getCenter()
{
    vector<double> w;
    w.push_back(xc0);
    w.push_back(yc0);
    w.push_back(zc0);
    return w;
}
    
/*---------------------------------------------------------------------------------------
 ---------------------------            METHODS USEFUL         --------------------------
 --------------------------------------------------------------------------------------*/

/* The function calculates the distance between the given point(x,y)
 and the thermal updraft center at a given altitude z.
 Effect of ambient winds and thermal drifting are considered.   */
double std_thermal::distToUpdraftCenter(double x, double y, double z)
{
    double r,xcz,ycz;
    xcz=xc0+windx*z; // xcz is the drifted center at alt z
    ycz=yc0+windy*z; // ycz is the drifted center at alt z
    r=sqrt((xcz-x)*(xcz-x) + (ycz-y)*(ycz-y));
    return(r);
}
    
// isAlive or not
bool std_thermal::isAlive(double currentTime)
{
    return currentTime-tBirth>lifeTime?(timeCoeff(currentTime)>0?1:0):1;
}
    
// This function returns the value of the thermal life cycle coefficient.
double std_thermal::timeCoeff(double currentTime)
{
    double tIndex,tau,D,T;
    tIndex = currentTime - tBirth;
    tau = tIndex-(lifeTime/2.0);
    T = (1+ksi)/lifeTime;
    D = (1-ksi) / (2.0*T);
    
    if (fabs(tau) <= D)
    {
        return(1);
    }
    else if (( D < fabs(tau) ) && ( fabs(tau) <= (1+ksi)/(2.0*T) ))
    {
        return(0.5 * (1+cos( (PI*T)/ksi * (fabs(tau) - D))));
    }
    else
    {
        return(0);
    }
}


/*---------------------------------------------------------------------------------------
 ---------------------------                 WIND              --------------------------
 --------------------------------------------------------------------------------------*/

// THE function to calculate the updraft
std_thermal& std_thermal::wind(double x,double y,double z,double t,vector<double> &w)
{
    //	z is above CML or below the ground
    if(z>zi || z<zc0){w[2]=0.0;}
    
    // Calculate distance based on Ambient Wind
    double r;
    r=distToUpdraftCenter(x,y,z);
    
    // Calulate the updraftValue according to 'thermal_model'
    switch(model)
    {
        case 1:
            // updraft calculation based on Allen model
            w[2] += Allen(r,z)*timeCoeff(t);
            break;
            
        case 2:
            // Calculation using Childress model
            w[2] += Childress(r, z)*timeCoeff(t);
            break;
            
        case 3:
            // updraft Calculation based on Lenschow with Gaussian distribution
            w[2]+=Lenschow(r,z,1)*timeCoeff(t);
            break;
            
        case 4:
            //updraft Calculation based on Lenschow with Geodon model
            w[2]+=Lenschow(r,z,0)*timeCoeff(t);
            break;
            
        case 5:
            // Calculating using Lawrance model
            Lawrance(w,x,y,z,t,windx,windy);
            w[0]*=timeCoeff(t);
            w[1]*=timeCoeff(t);
            w[2]*=timeCoeff(t);
            break;
            
        default:
            cout << "Caution! 'thermal_model' chosen isn't available!!" << endl;
            cout << "Method 'wind' in std_thermal didn't calculate vector w" << endl;
    }
    
    return *this;
}

/*---------------------------------------------------------------------------------------
 ---------------------------               MODELS               -------------------------
 --------------------------------------------------------------------------------------*/

/*----------------------------------
           Allen's model
 ----------------------------------*/
double std_thermal::Allen(double r, double z)
{
    double w_total;
    double w_peak,wd, w_ ,wL,s_wd;
    double r1,r2,r1_r2;
    double k1,k2,k3,k4;
    
    r1_r2 = 0.36;    //r1 / r2
    k1 = 1.4866;     //ki values valid for r1_r2 = r1 / r2 = 0.36
    k2 = 4.8354;
    k3 = -0.0320;
    k4 = 0.0001;
    
    r2 = 0.102 * pow((z / zi),1.0/3.0) * (1.0 - 0.25*z/zi) * zi;
    r2 = max(10.0, r2);
    r1 = r1_r2*r2;
    
    //Calculating w_ and W_peak using
    w_ = w_star * pow((z / zi),1.0/3.0) * (1.0 - 1.1*z/zi);
    w_peak = 3.0 * w_ * (r2*r2*r2 - r1*r2*r2) / (r2*r2*r2 - r1*r1*r1);
    
    //Checking for downdraft conditions
    wL = ((r2 < r) && (r < (2.0 * r2))) ? PI/6.0 * sin(PI * r / r2) : 0.;
    s_wd = ((0.5 < (z / zi)) && ((z / zi) < 0.9)) ? 2.5 * (z/zi - 0.5) : 0.;
    wd = s_wd * wL;
    
    // Total updraft
    w_total = w_peak * ( 1/(1+pow((fabs(k1*r/r2 + k3)),k2)) + k4 * r / r2 + wd );
    
    return (r>2.0*r2)?0.:w_total+w_total*normalLaw()/100.;
}
    
double std_thermal::integralWzAllen(double h)
{
    double output;
    if(h>1100)
        h=1100;
    output = 1/(2.64 * pow((h/1400.0),1.0/3.0) * (1.0 - 1.1*h/1400.0));
    return(output);
}
    
/*----------------------------------
         Childress model
----------------------------------*/
double std_thermal::Childress(double r, double z)
{
    /*
     Reference:
     An Empirical Model of thermal Updrafts Using Data Obtained From a Manned Glider
     Christopher E. Childress
     */
    
    double w_total, w_peak,wd,w_dec, w_;
    double d_T;
    double r1,r2;
    double z_star;
    z_star = z/zi;
    
    //Calculation of radius of the thermal
    d_T = zi*(0.4 * pow(z_star,1.0/3.0) * (1 - 0.5*z_star)) + (z*z_star*z_star - 0.6*z*z_star)/PI;
    r2 = 0.5*d_T;
    
    //Core downdraft radius
    r1 = 0.5*(0.17*d_T + 0.5*(z_star - 0.6)*d_T);
    
    // Calculating w_ and w_peak based on Allen
    w_ = w_star * pow((z / zi),1.0/3.0) * (1 - 1.1*z/zi);
    w_peak = 3 * w_ * (r2*r2*r2 - r1*r2*r2) / (r2*r2*r2 - r1*r1*r1);
    
    // Calculation of downdraft terms
    w_dec = (-zi/(1.275*w_star*w_star))*(12.192/z);
    wd = w_dec*(z_star + 0.45) - 0.5*w_peak;
    
    if(z_star>1.0)
    //The flight level is higher than the CBL so no effect of thermal updraft.
    {
        return(0.0);
    }
    
    //Calculation of Updraft based on equations 14,15,16 from Childress
    
    if(z_star<0.5 && r<=r2)
    {
        w_total = w_peak*cos((r/r2)*PI*0.5);
    }
    else if(z_star<0.9)
    {
        if(r<r1)
        {
            w_total = wd*cos((r/r1)*PI*0.5);
        }
        else if(r1<=r && r<=(r1+r2))
        {
            w_total = w_peak*sin((r-r1)/r2 * 1.212 *PI);
            
            if(w_total<0.0)
                w_total=0.0;
        } 
        else
        {
            w_total = 0.0;
        }
    }
    else
    {
        if(r<r1)
        {
            w_total = 0.5*wd*cos((r/r1)*PI*0.5);
        }
        else if(r1<=r && r<=(r1+r2))
        {
            w_total = (1-z_star)*w_peak*sin((r-r1)/r2 * 1.212 *PI);
            
            if(w_total<0.0)
                w_total=0.0;
        }
        else
        {
            w_total=0.0;
        }				
    }
    return(w_total);
}
    
/*----------------------------------
          Lenschow model
 ----------------------------------*/
double std_thermal::Lenschow(double r, double z, bool choice)
{
    if(z>zi)
    {
        //	The flight level is higher than the CBL so no effect of thermal updraft.
        return(0.0);
    }
    
    double w_total,w_,d;
    double w_peak;
    
    //The diameter of the thermal is given by
    d = 0.16 * pow((z / zi),1.0/3.0) * (1 - 0.25*z/zi) * zi;
    
    //normalized updraft velocity
    w_= w_star*pow((z / zi),1.0/3.0) * (1 - 1.1*z/zi);
    
    //variance of the updraft velocity
    //double var = 1.8*pow((z / zi),2.0/3.0) * (1-0.8*z/zi)*(1 - 0.8*z/zi);
    
    //w_peak assuming a gaussian distribution is
    w_peak=w_;//*(2.0/pow(PI,0.5));
    
    if (choice==1)
    {
        w_total=w_peak*exp(-(4*r*r/(d*d)));
    }
    else
    {
        w_total=w_peak*exp(-(4*r*r/(d*d)))*(1-(4*r*r/(d*d)));
    }
    return(w_total);
}
    
/*----------------------------------
         Lawrance's model
----------------------------------*/
void std_thermal::Lawrance(vector<double> &w,double x, double y, double z, double t,double Wx,double Wy)
{
    double x0=0.0,y0=0.0,z0=0.0,xt,yt,zt;
    
    double rT,r1,r1_r2=0.36, k=3.0;
    rT = 0.102 * pow((z / zi),1.0/3.0) * (1.0 - 0.25*z/zi) * zi;
    rT = max(10.0, rT);
    //cout<<r2<<endl;
    r1 = r1_r2*rT;
    double w_,w_core ;
    
    //Calculating w_ and W_peak using
    w_ = w_star * pow((z / zi),1.0/3.0) * (1.0 - 1.1*z/zi);
    w_core = 3.0 * w_ * (rT*rT*rT - r1*rT*rT) / (rT*rT*rT - r1*r1*r1);
    
    z0 = 800.0;
    
    if(z0<k*rT)
    //The bubble has not detached from the ground yet
    {
        x0=xc0;
        y0=yc0;
    }
    else
    //The bubble is completely formed and it can detach itself from the ground and move along with the wind
    {
        x0=xc0; //+ simpsons(integralWzAllen(z),100.0,z,1000)*Wx;
        y0=yc0; //+ simpsons(integralWzAllen(z),100.0,z,1000)*Wy;
    }
    
    zt=z-z0;
    xt=x-x0;
    yt=y-y0;
    
    double dH;
    dH=sqrt(xt*xt + yt*yt);
    
    //Calculation of Wz
    if(dH==0)
    {
        w[2]+= w_core * 0.5*(cos(PI*zt/(k*rT))+1.0);
    }
    else if(dH<=2*rT)
    {
        w[2]+=0.5*(w_core*rT)/(PI*dH) * sin(PI*dH/rT) * (cos(PI*zt/(k*rT))+1.0);
    }
    else
    {
        w[2]+=0.0;
    }
    
    if(fabs(zt)>(k*rT))
    {
        w[2]+=0.0;
    }
    
    //calculation of Wx Wy
    if(dH!=0.0||dH<rT)
    {
        w[0]+=-(w[2]*zt*xt)/(dH*(dH-rT)*k*k);
        w[1]+=-(w[2]*zt*yt)/(dH*(dH-rT)*k*k);
    }
    else if(dH==rT)
    {
        w[0]+=-w_core/(2.0*k*rT) * (1.0 + cos(PI*z/(k*rT)));
        w[1]+=-w_core/(2.0*k*rT) * (1.0 + cos(PI*z/(k*rT)));
    }
    else
    {
        w[0]+=0;
        w[1]+=0;
    }
}
    
double std_thermal::simpsons( double (*f)(double x), double a, double b, int n)
{
    double h = (b - a) / n;
    double x;
    double r;
    char m = 0;
    double s = 0.0;
    
    for (x = a; x <= b; x+=h)
    {
        r = f(x);
        if (x == a || x == b)
        {
            s += r;
        }
        else
        {
            m = !m;
            s += r * (m+1) * 2.0;
        }
    }
    return s * (h/3.0);
}
    
}
