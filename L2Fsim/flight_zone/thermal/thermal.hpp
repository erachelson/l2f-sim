#ifndef L2FSIM_THERMAL_HPP_
#define L2FSIM_THERMAL_HPP_

#include <vector>

namespace L2Fsim {

/// The abstract class thermal calculates the wind vector w linked with a thermal, at time t.

class thermal {

    /*--------------------------------------
     ------------- Attributes --------------
     -------------------------------------*/

protected:
    /// The horizontal wind in the flight zone, implying the drift of the thermal [m/s].
    double windx,windy;

    /*--------------------------------------
     --------------- Methods ---------------
     -------------------------------------*/

public:

    // get
    ///Get the peak velocity of the thermal.
    virtual double getw_star() =0;
    ///Get the birth time of the thermal.
    virtual double gettBirth() =0;
    ///Get the life time of the thermal.
    virtual double getlifeTime() =0;
    ///Get the height of the thermal.
    virtual double getzi() =0;
    ///Get the model of the thermal.
    virtual int getModel() =0;
    ///Get the shape factor linked to thermal life cycle.
    virtual double getksi() =0;
    ///Get the vector of thermal centers.
    virtual std::vector<double> getCenter() =0;

    // method
    ///Get the state of the thermal : alive or not.
    /**
		\param t the simulated time.
	*/
    virtual bool isAlive(double t) =0;

    ///Get the age of the thermal. Affect the peak velocity.
    /**
		\param t the simulated time.
	*/
    virtual double timeCoeff(double t) =0;

    ///Get the Euclidean distance between a point (x,y,z) and the center of the thermal. If z =! 0 then it considers the drift of the thermal center at the height z.
    /**
		\param x the x-axis coordinate of the point.
		\param y the y-axis coordinate of the point.
		\param z the z-axis coordinate of the point.
	*/
    virtual double distToUpdraftCenter(double x, double y, double z) =0;

    // set
    ///Set the wind speed onto the x-axis and y-axis.
    /**
		\param wx the wind speed onto the x-axis.
		\param wy the wind speed onto the y-axis.
	*/
    void setwind(double wx,double wy)
    {
        windx=wx;
        windy=wy;
    }

    // wind
    /// Computes the wind vector w, at point (x,y,z), at time t.
	/**
		\param x a horizontal coordinate of the system (x,y,z).
		\param y a horizontal coordinate of the system (x,y,z).
		\param z the vertical coordinate of the system (x,y,z).
		\param t the time.
		\param w the wind vector: windx, windy, windz.
	*/
	virtual thermal& wind(double x,double y,double z,double t,std::vector<double> &w) =0;
};

}

#endif
