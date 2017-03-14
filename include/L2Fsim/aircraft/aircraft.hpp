#ifndef L2F_AIRCRAFT_HPP_
#define L2F_AIRCRAFT_HPP_

#include <vector>
#include <L2F/flight_zone/flight_zone.hpp>

namespace L2F {

/// The abstract class Aircraft implementing the aircraft models.
/**
 * An aircraft holds three important concepts:
 * - it has an internal state x which characterizes uniquely the aircraft's configuration at a given time;
 * - it has a dynamics function f, such that, given a command u, x'=f(x,u);
 * - it has an observation function on the current state.
 */
class aircraft {
	/* methods */
public:
	/// Updates the aircraft's state vector.
	/**
		\param st the aircraft's state vector.
		\param cd the aircraft's command vector: alpha, beta, sigma.
		\param fz the flight zone.
		\param t the time.
		\param statedot the derivatives vector: x_dot, y_dot, z_dot, V_dot, khi_dot, gamma_dot.
	*/
	virtual const aircraft& state_dynamics(const std::vector<double> &st,
		 							 const std::vector<double> &cd,
									 flight_zone& fz,
									 double t,
									 std::vector<double> &statedot) const =0;
	/// Gets the observation vector.
	/**
		\param observation the observation vector.
	*/
	virtual aircraft& observation(std::vector<double> &observation)  =0;
	/// Gets the aircraft's state vector.
	/**
		\param state the vector to copy the aircraft's state in.
	*/
	virtual aircraft& get_state(std::vector<double> &state) =0;
	/// Sets the aircraft's state vector.
	/**
		\param state the vector to be copied to the aircraft's state.
	*/
	virtual aircraft& set_state(const std::vector<double> &state) =0;
    /**
        \param The model decide if the value of the state are out of the range of the model
     */
    virtual aircraft& is_in_model(std::vector<double> &state,std::vector<double> &command) = 0;

};

}

#endif
