#ifndef L2FSIM_EULER_INTEGRATOR_HPP_
#define L2FSIM_EULER_INTEGRATOR_HPP_

#include <L2Fsim/stepper/stepper.hpp>
#include <cstdio>
#include <L2Fsim/utils/utils.hpp>
#include <stdio.h>
#include <stdlib.h>

namespace L2Fsim {

class euler_integrator : public stepper {
	/* attribute */
    public:
	/** Integration dt */
	double delta_t;
    /** Current time of the simulation */
	double time;
    /** File name for the log*/
    std::string name_file = "data/data_plane.txt";
    /** File name for the log*/
    std::string name_file1 = "data/energy.txt";

    /** Precedant observation */
    std::vector<double> obs_old;

	/* methods */
public:
	/** Constructor */
	euler_integrator(double dt_t=0.0001) : delta_t(dt_t) {
        time = 0;
    }
	/** Stepping operator 
     * We use the the euler integrator over the time dt
     * We use the aircraft class to calculte d_state and we obtain the new state
     */
	 void operator()(flight_zone &fz,
							aircraft &ac,
							pilot &pl,
							double dt) {

		// X dot -> parameter of aircraft ? (-> no recreate each time)
		std::vector<double> xdot;
		std::vector<double> x;
		std::vector<double> obs;
		std::vector<double> u;

		ac.get_state(x);
		ac.observation(obs);
		if (time==0.) {obs_old=obs;}

		pl(obs, u);

        std::vector<double> concat;
        concat.reserve( x.size() + u.size() ); // preallocate memory
        concat.insert( concat.end(), x.begin(),x.end() );
        concat.insert( concat.end(), u.begin(), u.end() );

    	//printf("saving data\n");
    	save_state_into_file(name_file,concat,time);
    	save_energy_into_file(name_file1,time,dt,obs_old,obs);

		for(int n=0; n<dt/delta_t; ++n) {

			//printf("Aircraft simulation\n" );
			ac.get_state(x);

			ac.state_dynamics(x,u, fz,time, xdot);
			//printf("State dynamic\n" );

			for(int i=0; i<x.size(); ++i) {
				x[i] += xdot[i] * delta_t;
			}

			ac.set_state(x);

			time += delta_t;

		}
		//printf("Finishing of the first step\n" );
         
		obs_old=obs;
         
         // We verify if the model is correct :
         ac.is_in_model(x,u);
         
     }
};

}

#endif
