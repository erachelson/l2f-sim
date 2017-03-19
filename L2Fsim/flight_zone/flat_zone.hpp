#ifndef L2FSIM_FLAT_ZONE_HPP_
#define L2FSIM_FLAT_ZONE_HPP_

#include <vector>

namespace L2Fsim{

class flat_zone : public flight_zone
{
    
    /*--------------------------------------
     ------------- Attributes --------------
     -------------------------------------*/
    
protected:
    
	double windx;
	double windy;
    
    
    /*--------------------------------------
     --------------- Methods ---------------
     -------------------------------------*/

public:
    
    // wind
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
