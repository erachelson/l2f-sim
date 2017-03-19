#ifndef L2FSIM_RKF45_HPP_
#define L2FSIM_RKF45_HPP_

#include <L2Fsim/stepper/stepper.hpp>

namespace L2Fsim {

class RKF45 : public stepper {
	/* methods */
public:
	int nb_steps;
	/** State dimension */
	int xdim;
	/** Control dimension */
	int udim;
	double time;
public:
	virtual void operator()(flight_zone &fz,
							aircraft &ac,
							pilot &pl,
							double dt) {

// Initialisation des variables
					std::vector<double> xdot(xdim);
				  std::vector<double> x;
					std::vector<double> obs;
					std::vector<double> u;
					ac.get_state(x);
					ac.observation(obs);
// Le pilot fait son action
					pl(obs, u);

// Initialisation temporel
					double delta_t;
					delta_t = dt/((double) nb_steps);
					double time;


// Initialisation des paramètres pour Runge-Kutta-Fehlberg
					std::vector<double>  k1(xdim);
					std::vector<double>  k2(xdim);
					std::vector<double>  k3(xdim);
					std::vector<double>  k4(xdim);
					std::vector<double>  k5(xdim);
					std::vector<double>  k6(xdim);

					std::vector<double> temp(xdim);

					for(int n=0; n<nb_steps; ++n) {

						ac.get_state(x);
						ac.state_dynamics(x,u, fz,time, k1);

						for(int i=0; i<x.size(); ++i) {

							k1[i] = delta_t * k1[i];
							temp[i] = x[i] + 1/4 * k1[i];

						}

						ac.state_dynamics(temp,u, fz,time + 1/4*delta_t, k2);

						for(int i=0; i<x.size(); ++i) {

							k2[i] = delta_t * k2[i];
							temp[i] = x[i] + 3/32 * k1[i] + 9/32 * k2[i];

						}

						ac.state_dynamics(temp,u, fz,time + 3/8*delta_t, k3);

						for(int i=0; i<x.size(); ++i) {

							k3[i] = delta_t * k3[i];
							temp[i] = x[i] + 1932/2197 * k1[i] - 7200/2191 * k2[i] + 7296/2197 * k3[i];

						}

						ac.state_dynamics(temp,u, fz,time + 12/13*delta_t, k4);

						for(int i=0; i<x.size(); ++i) {

							k4[i] = delta_t * k4[i];
							temp[i] = x[i] + 439/216 * k1[i] - 8 * k2[i] + 3680/513 * k3[i] - 845/4104 * k4[i];

						}

						ac.state_dynamics(temp,u, fz,time + 1*delta_t, k5);

						for(int i=0; i<x.size(); ++i) {

							k5[i] = delta_t * k5[i];
							temp[i] = x[i] - 8/27 * k1[i] + 2 * k2[i] - 3544/2565 * k3[i] + 1859/4104 * k4[i] - 11/40*k5[i];

						}

						ac.state_dynamics(temp,u, fz,time + 1/2*delta_t, k6);

						for(int i=0; i<x.size(); ++i) {

							k6[i] = delta_t * k6[i];

						}

						// Then we reassemble everything to calculate the true value :
						for(int i=0; i<x.size(); ++i) {

							x[i] = x[i]+25/216 * k1[i] + 1408/2565 * k3[i] + 2197/ 4107 * k4[i] - 1/5 * k5[i];

						}

						ac.set_state(x);
						time += delta_t;

					}

	}
};

}

#endif
