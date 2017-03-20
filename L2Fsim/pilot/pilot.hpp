#ifndef L2FSIM_PILOT_HPP_
#define L2FSIM_PILOT_HPP_

#include <vector>

namespace L2Fsim {

class pilot {
	/* methods */
public:
	virtual pilot& operator()(std::vector<double> &observation,
							  std::vector<double> &command) =0;
    
    virtual pilot& out_of_range(std::vector<double> &observation,
                              std::vector<double> &command) =0;
    
};
}

#endif
