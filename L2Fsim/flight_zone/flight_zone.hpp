#ifndef L2FSIM_FLIGHT_ZONE_HPP_
#define L2FSIM_FLIGHT_ZONE_HPP_

#include <vector>

namespace L2Fsim{

class flight_zone
{
    
    /*--------------------------------------
     --------------- Methods ---------------
     -------------------------------------*/
    
public:
    
    //Computes the wind vector w, at point (x,y,z), at time t
	virtual flight_zone& wind(double x, double y, double z, double t, std::vector<double> &w) =0;
	
    //Computes the z coordinate of the surface at (x,y)
    flight_zone& ground(double x, double y, double &z)
    {
        z=0.;
        return *this;
    }
};

}

#endif
