#ifndef L2FSIM_SIMULATION_HPP_
#define L2FSIM_SIMULATION_HPP_

#include <L2Fsim/flight_zone/flight_zone.hpp>
#include <L2Fsim/aircraft/aircraft.hpp>
#include <L2Fsim/stepper/stepper.hpp>
#include <L2Fsim/pilot/pilot.hpp>

namespace L2Fsim {

class simulation{
public:
	/* Attributes */
	flight_zone *fz;
	aircraft *ac;
	stepper *st;
	pilot *pl;
	/* Methods */
public:
	/**
	 * Constructor
	 */
	
	/**
	 * Stepping function
	 */
	void step(double dt) {
		(*st)(*fz, *ac, *pl, dt);
	}
};

}

#endif
