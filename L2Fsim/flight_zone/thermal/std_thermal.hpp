#ifndef L2FSIM_STD_THERMAL_HPP_
#define L2FSIM_STD_THERMAL_HPP_

#include <L2Fsim/flight_zone/thermal/thermal.hpp>
#include <L2Fsim/utils/utils.hpp>
#include <iostream>
#include <vector>
#include <math.h>
#include <ctime>
#include <vector>

namespace L2Fsim {

/// The abstract class std_thermal is a subclass of thermal. It is a specialization of thermal.

class std_thermal : public thermal {

    /*--------------------------------------
     ------------- Attributes --------------
     -------------------------------------*/

protected:
    /// The chosen model for the thermals.
    int model;
    /// The birth time of the thermal.
    const double tBirth;
    /// Position of the thermal center at tBirth.
    const int xc0,yc0,zc0;
    /// Convective velocity scaling parameter [m/s].
    double w_star;
    /// Convective mixing layer thickness [m].
    double zi;
    /// The life time of the thermal.
    double lifeTime;
    /// The shape factor linked to thermal life cycle.
    double ksi;

    /*--------------------------------------
     --------------- Methods ---------------
     -------------------------------------*/

public:

    /// Constructor of a thermal.
    /**
        \param model the chosen model for the thermal.
        \param tB the birth time of the thermal.
		\param XC0 the x-axis coordinate of the thermal center at tBirth.
		\param YC0 the x-axis coordinate of the thermal center at tBirth.
		\param ZC0 the x-axis coordinate of the thermal center at tBirth.
		\param Zi the convective mixing layer thickness.
		\param wstar the convective velocity scaling parameter or peak velocity. Random draw between 13 values.
		\param lifetime the life time of the thermal. Random draw between 600 and 1200.
		\param ksi the shape factor linked to thermal life cycle. Random draw between 0.1 and 0.35.
		\param read a boolean that specifies if you create or read the thermal.
	*/
    std_thermal(int model,double tB=0.,double XC0=0.,double YC0=0.,double ZC0=0.,
                double Zi=1200.,double wstar=2.5, int lifetime=1200, double ksi=0.2,bool read=false);

    /// Destructor.
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
    ///Compute the Allen's model.
    /**
		\param r the radius of the thermal.
		\param z the height of the thermal.
	*/
    double Allen(double r,double z);

    ///Compute the Childress's model.
    /**
		\param r the radius of the thermal.
		\param z the height of the thermal.
	*/
    double Childress(double r,double z);

    ///Compute the Lenschow's model.
    /**
		\param r the radius of the thermal.
		\param z the height of the thermal.
	*/
    double Lenschow(double r, double z, bool choice);

    ///Compute the Lawrance's model.
    /**
		\param r the radius of the thermal.
		\param z the height of the thermal.
	*/
    void Lawrance(std::vector<double> &w,double x, double y, double z, double t,double Wx,double Wy);
};

}

#endif
