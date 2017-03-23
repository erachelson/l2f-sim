#ifndef L2FSIM_FLAT_ZONE_HPP_
#define L2FSIM_FLAT_ZONE_HPP_

#include <vector>

namespace L2Fsim{

/// The abstract class flat_zone is a subclass of flight_zone. It implements a flat ground at sea level (z=0) with a horizontal wind.


class flat_zone : public flight_zone
{

    /*--------------------------------------
     ------------- Attributes --------------
     -------------------------------------*/

protected:
    /// The wind speed onto the x-axis.
	double windx;
	/// The wind speed onto the y-axis.
	double windy;


    /*--------------------------------------
     --------------- Methods ---------------
     -------------------------------------*/

public:

    /// Computes the wind vector w, at point (x,y,z), at time t.
	/**
		\param x a horizontal coordinate of the system (x,y,z).
		\param y a horizontal coordinate of the system (x,y,z).
		\param z the vertical coordinate of the system (x,y,z).
		\param t the time.
		\param w the wind vector: windx, windy, windz.
	*/

	virtual flat_zone& wind(double x, double y, double z, double t, std::vector<double> &w)
    {
		std::vector<double>(3,0.).swap(w);
		w.at(0) = windx;
		w.at(1) = windy;
		return *this;
	}

};

}

#endif
