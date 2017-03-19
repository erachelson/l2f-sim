#include <L2Fsim/simulation.hpp>
#include <L2Fsim/aircraft/BeelerGlider.hpp>
#include <L2Fsim/pilot/passive_pilot.hpp>
#include <L2Fsim/flight_zone/flat_thermal_soaring_zone.hpp>
#include <L2Fsim/flight_zone/flat_zone.hpp>
#include <L2Fsim/stepper/euler_integrator.hpp>
#include <L2Fsim/stepper/RKF45.hpp>
#include <L2Fsim/pilot/pilot_genQlearn.hpp>
#include <L2Fsim/pilot/pilot_etienne.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace L2Fsim;
using namespace std;

int main()
{
    srand (time(NULL));
    double time_limit = 1000;
    
    /*--------------------------------------
     ------------- FLIGHT ZONE -------------
     -------------------------------------*/
    
    // read config1.txt
    flat_zone my_zone;
    //flat_thermal_soaring_zone my_zone("DATA/config1.txt");
    
    /*--------------------------------------
     --------------- AIRCRAFT --------------
     -------------------------------------*/
	double m = 1.35;
	double ws = 2.;
	double ar=16.;
	
    // Type of aircraft
	BeelerGlider my_glider(m,ws,ar);

    /*--------------------------------------
     ---------------- PILOT ----------------
     -------------------------------------*/
    //passive_pilot my_pilot;
    pilot_genQlearn my_pilot;

    /*--------------------------------------
     --------------- STEPPER ---------------
     -------------------------------------*/
	euler_integrator my_stepper;

    /*--------------------------------------
     ------------- SIMULATION --------------
     -------------------------------------*/
    // Initialization of the main class
    simulation mysim;

    // We attach those class the main class simulation
	mysim.fz = &my_zone;
	mysim.ac = &my_glider;
	mysim.st = &my_stepper;
	mysim.pl = &my_pilot;

    // Initialization of the main parameter
	// Here we choose the time limit for our simulation (second)
	double dt = 0.01; // 1 ms as dt
	double time_current = 0;

	// Initilation of the position // TODO with a config file
	double x0 = -500.;// m
	double y0 = 0.;// m
	double z0 = 1000.;// m
	double V0 = 20.; // m/s
	double gamma0 = 0.; // rad
	double khi0 = 0.; // rad

	std::vector<double> state_init = {x0,y0,z0,V0,gamma0,khi0};
	mysim.ac->set_state(state_init);

	// for loop on the stepper
    
	while(time_current < time_limit)
    {

        printf("time stamp : %f\n", time_current);
		// We increment the system -> evolution of the aircraft
		mysim.step(dt);
		time_current += dt;

	}
     



}
