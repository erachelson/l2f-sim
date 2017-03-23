#ifndef L2FSIM_FLIGHT_ZONE_HPP_
#define L2FSIM_FLIGHT_ZONE_HPP_

#include <vector>

namespace L2Fsim{

/// The abstract class flight_zone implementing the zone in which the aircraft moves.
/**
 * A flight zone holds two important concepts:
 * - it has a characterization of the wind w in the flight zone at a given time;
 * - it has an altitude z of the ground surface at all points in the flight zone;
 */

class flight_zone
{

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
		\param w the wind vector: windx, windy, windz. These three components are the projections of the wind speed onto the coordinate system (x,y,z).
	*/

	virtual flight_zone& wind(double x, double y, double z, double t, std::vector<double> &w) =0;

    /// Computes the z coordinate of the ground at (x,y). Set at 0 as default value.
	/**
        \param x a horizontal coordinate of the system (x,y,z).
		\param y a horizontal coordinate of the system (x,y,z).
		\param z the vertical coordinate of the system (x,y,z).
	*/

    flight_zone& ground(double x, double y, double &z)
    {
        z=0.;
        return *this;
    }
};

}

#endif
