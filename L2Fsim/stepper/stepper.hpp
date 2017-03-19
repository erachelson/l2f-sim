#ifndef L2FSIM_STEPPER_HPP_
#define L2FSIM_STEPPER_HPP_

#include <L2Fsim/flight_zone/flight_zone.hpp>
#include <L2Fsim/aircraft/aircraft.hpp>
#include <L2Fsim/pilot/pilot.hpp>

namespace L2Fsim {

struct stepper {
	/* methods */
public:
	virtual void operator()(flight_zone &fz,
							aircraft &ac,
							pilot &pl,
							double dt) =0;
};

}

#endif
