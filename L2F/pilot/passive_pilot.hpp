#ifndef L2F_PASSIVE_PILOT_HPP_
#define L2F_PASSIVE_PILOT_HPP_

#include <L2F/pilot/pilot.hpp>
#include <vector>

namespace L2F {

class passive_pilot : public pilot {
public:
	/* attributes */
	unsigned int output_dim;
	/* methods */
	passive_pilot& operator()(std::vector<double> &obs,
									  std::vector<double> &command) {
		std::vector<double>(output_dim, 0.).swap(command);
        command.at(2) = 0.1;
		return *this;
	}
};

}

#endif
