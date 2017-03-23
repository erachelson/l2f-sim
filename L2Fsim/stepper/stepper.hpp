#ifndef L2FSIM_STEPPER_HPP_
#define L2FSIM_STEPPER_HPP_

#include <L2Fsim/flight_zone/flight_zone.hpp>
#include <L2Fsim/aircraft/aircraft.hpp>
#include <L2Fsim/pilot/pilot.hpp>

namespace L2Fsim {

struct stepper {
	/* methods */
public:
    /// The stepper is the temporal integrator of the model.
    /** The stepper take the thermal model, the aircraft model and the command law (pilot)
     * and make an intégration on the different componant.
     */
	virtual void operator()(flight_zone &fz,
							aircraft &ac,
							pilot &pl,
							double dt) =0;
};

}

#endif
