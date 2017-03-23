#ifndef L2FSIM_PILOT_HPP_
#define L2FSIM_PILOT_HPP_

#include <vector>

namespace L2Fsim {

class pilot {
	/* methods */
public:
    /// this method take the observation of the pilot (what he sees) (obs) and return his command (command) when he is in the thermal zone
	virtual pilot& operator()(std::vector<double> &observation,
							  std::vector<double> &command) =0;
    /// this method take the observation of the pilot (what he sees) (obs) and return his command (command) when he's out of the range of the thermal zone
    virtual pilot& out_of_range(std::vector<double> &observation,
                              std::vector<double> &command) =0;
    
};
}

#endif
