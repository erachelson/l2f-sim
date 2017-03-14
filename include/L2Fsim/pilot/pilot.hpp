#ifndef L2F_PILOT_HPP_
#define L2F_PILOT_HPP_

#include <vector>

namespace L2F {

class pilot {
	/* methods */
public:
	virtual pilot& operator()(std::vector<double> &observation,
							  std::vector<double> &command) =0;
};
}

#endif
