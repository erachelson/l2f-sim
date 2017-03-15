#ifndef L2F_STEPPER_HPP_
#define L2F_STEPPER_HPP_

#include <L2F/flight_zone/flight_zone.hpp>
#include <L2F/aircraft/aircraft.hpp>
#include <L2F/pilot/pilot.hpp>

namespace L2F {

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
