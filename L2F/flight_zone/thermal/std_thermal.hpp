#ifndef L2F_STD_THERMAL_HPP_
#define L2F_STD_THERMAL_HPP_

#include <L2F/flight_zone/thermal/thermal.hpp>
#include <L2F/utils/utils.hpp>
#include <iostream>
#include <vector>
#include <math.h>
#include <ctime>
#include <vector>

namespace L2F {

class std_thermal : public thermal {
    
    /*--------------------------------------
     ------------- Attributes --------------
     -------------------------------------*/
    
protected:
    
    int model;             // chosen model for the thermals
    const double tBirth;   // TimeBirth of the thermal
    const int xc0,yc0,zc0; // Position of the thermal center at tBirth
    double w_star;         // Convective velocity scaling parameter [m/s]
    double zi;             // Convective mixing layer thickness [m]
    double lifeTime;       // lifetime of a thermal
    double ksi;            // ksi shape factor linked to thermal life cycle
    
    /*--------------------------------------
     --------------- Methods ---------------
     -------------------------------------*/
    
public:
    
    // constructor
    std_thermal(int model,double tB=0.,double XC0=0.,double YC0=0.,double ZC0=0.,
                double Zi=1200.,double wstar=2.5, int lifetime=1200, double ksi=0.2,bool read=false);
    
    // destructor
    ~std_thermal();
    
    // getteurs
    double getw_star() {return w_star;}
    double gettBirth() {return tBirth;}
    double getlifeTime() {return lifeTime;}
    double getzi() {return zi;}
    int getModel() {return model;}
    double getksi() {return ksi;}
    std::vector<double> getCenter();

    // methods useful
    double distToUpdraftCenter(double x, double y, double z);
    bool isAlive(double currentTime);
    double timeCoeff(double currentTime);
    double integralWzAllen(double h);
    double simpsons( double (*f)(double x), double a, double b, int n);
    
    // wind
    virtual std_thermal& wind(double x,double y,double z,double t,std::vector<double> &w);
    
    
    /*--------------------------------------
     ------------ Thermal Models -----------
     -------------------------------------*/
    
    double Allen(double r,double z);
    double Childress(double r,double z);
    double Lenschow(double r, double z, bool choice);
    void Lawrance(std::vector<double> &w,double x, double y, double z, double t,double Wx,double Wy);
};

}

#endif
