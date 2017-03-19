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
	passive_pilot& operator()(std::vector<double> &obs,
									  std::vector<double> &command) {
		std::vector<double>(output_dim, 0.).swap(command);
        command.at(2) = 0.3;
		return *this;
	}
};

}

#endif
