#ifndef L2FSIM_PASSIVE_PILOT_HPP_
#define L2FSIM_PASSIVE_PILOT_HPP_

#include <L2Fsim/pilot/pilot.hpp>
#include <vector>

namespace L2Fsim {

class passive_pilot : public pilot {
public:
	/* attributes */
	unsigned int output_dim=3;
	/* methods */
    /// A passive pilot implementation. A pilot that send a constant response to any given aircraft state.
    /** The pilot take the observation vector (what he sees in the aircraft) and return a null command vector :
     * [alpha beta sigma] = [0 0 0]
     */
	passive_pilot& operator()(std::vector<double> &obs,
									  std::vector<double> &command) {
		std::vector<double>(output_dim, 0.).swap(command);
        command.at(2) = 0.;
		return *this;
	}
    
    /// In the outrange of the thermal the pilot just go for the command : [alpha beta sigma]= [0 0 0.4] trying to be in the thermal zone again
    passive_pilot& out_of_range(std::vector<double> &obs,
                                  std::vector<double> &command)
    {
        std::vector<double>(3, 0.).swap(command);
        command.at(2) = 0.4;
        
    }
};

}

#endif
