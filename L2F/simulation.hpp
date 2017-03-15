#ifndef L2F_SIMULATION_HPP_
#define L2F_SIMULATION_HPP_

#include <L2F/flight_zone/flight_zone.hpp>
#include <L2F/aircraft/aircraft.hpp>
#include <L2F/stepper/stepper.hpp>
#include <L2F/pilot/pilot.hpp>

namespace L2F {

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
