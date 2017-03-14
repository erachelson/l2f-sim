#include <L2F/simulation.hpp>
#include <L2F/aircraft/BeelerGlider.hpp>
#include <L2F/pilot/passive_pilot.hpp>
#include <L2F/flight_zone/flat_thermal_soaring_zone.hpp>
#include <L2F/stepper/euler_integrator.hpp>
#include <L2F/stepper/RKF45.hpp>
#include <L2F/pilot/pilot_genQlearn.hpp>
#include <L2F/pilot/pilot_etienne.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace L2F;
using namespace std;

int main()
{
    srand (time(NULL));
    /*--------------------------------------
     ------------- FLIGHT ZONE -------------
     -------------------------------------*/
    double time_limit = 1000;
    double deltaT     = 100;  // actualization of thermals
    int minX  = -100;        // definition of box
    int maxX  = 100;         // definition of box
    int minY  = -100;        // definition of box
    int maxY  = 100;         // definition of box
    int minZ  = 0;            // definition of box
    int maxZ  = 2000;         // definition of box
    double wx = 0.;           // definition of wind
    double wy = 0.;           // definition of wind
    int zi    = 1400.;        // definition of height of all thermals

    /* Initialization of a windfield */
    // Constructor 1 : flat_zone(Tend,minX,maxX,minY,maxY,minZ,maxZ,wx,wy,zi):tstart=0.
    //flat_thermal_soaring_zone my_zone(time_limit,minX,maxX,minY,maxY,minZ,maxZ,wx,wy,zi);

    // Constructor 2 : flat_zone(filename)
    flat_thermal_soaring_zone my_zone("DATA/config1.txt");
    // flat_thermal_soaring_zone my_zone("DATA/configUneSeuleTherm.txt");

    /* Simulate a scenario of thermals */
    // L2F::createScenario(deltaT,model);
    // 1 : Allen model
    // 2 : Childress model
    // 3 : Lenschow model with Gaussian distribution
    // 4 : Lenschow with Geodon model
    // 5 : Lawrance model
    // my_zone.createScenario(deltaT,1);

    /* Write minimum info of scenario to rebuilt it */
    // my_zone.saveConfig("DATA/config2.txt");
    // my_zone.saveConfigToCSV("DATA/config2CSV.txt");

    /* Write the whole wind DATA info in a file*/
    // writeScenario(deltaT,deltax,deltay,zslice,filename)
    my_zone.writeScenario(deltaT,5.,5.,650.,"data/wind.txt");

    /*--------------------------------------
     --------------- AIRCRAFT --------------
     -------------------------------------*/
	double m = 1.35;
	double ws = 2.;
	double ar=16.;
	
    // Type of aircraft
	BeelerGlider my_glider( m,ws,ar);

    /*--------------------------------------
     ---------------- PILOT ----------------
     -------------------------------------*/
    //passive_pilot my_pilot;

    pilot_genQlearn my_pilot;
    my_pilot.output_dim = 3;

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
	double dt = 0.001; // 1 ms as dt
	double time_current = 0;

	// Initilation of the position // TODO with a config file
	double x0 = -800.;// m
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
