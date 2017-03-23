#ifndef L2FSIM_FLAT_THERMAL_SOARING_ZONE_HPP_
#define L2FSIM_FLAT_THERMAL_SOARING_ZONE_HPP_

#include <L2Fsim/flight_zone/flight_zone.hpp>
#include <L2Fsim/flight_zone/flat_zone.hpp>
#include <L2Fsim/flight_zone/thermal/std_thermal.hpp>
#include <string>
#include <fstream>
#include <sstream>

namespace L2Fsim {

/// The abstract class flat_thermal_soaring_zone is a subclass of flat_zone.
/**
 * This class implements the z-components of the wind vector by introducing the thermals in the flat zone.
 * All over the simulation, thermals appeared and disappeared randomly in the flat zone.
 */

class flat_thermal_soaring_zone : public flat_zone
{


    /*--------------------------------------
     ------------- Attributes --------------
     -------------------------------------*/

protected:
    /// The start time of the simulation. Set at 0.
    double tstart;
    /// The end time of the simulation.
    double tend;

    ///The boundaries of the flat_thermal_soaring_zone onto the x-axis.
    int minX,maxX;
    ///The boundaries of the flat_thermal_soaring_zone onto the y-axis.
    int minY,maxY;
    ///The boundaries of the flat_thermal_soaring_zone onto the z-axis.
    int minZ,maxZ;

    ///The radius of a thermal.
    double dmin=200.;

    ///The height of a thermal.
    double zi;

    ///The list of the thermals created in the simulation.
    std::vector<thermal *> thermals;

    /*--------------------------------------
     --------------- Methods ---------------
     -------------------------------------*/

    //private
    int nbMaxThermals();
    bool createThermalCenter(std::vector<double>& center, double t);
    int numberAliveAtTime(double t);
    double environSink(double z, double t);

public:

    ///Constructor. Create an empty zone defined by its dimension, wind, end time and height of thermals.
    flat_thermal_soaring_zone(double Tend,int miX,int maX,int miY,int maY,int miZ,int maZ,double wx,double wy,double zi);

    ///Constructor. Read a saved thermal scenario in order to play an identical simulation of the flight zone.
    /**
		\param filename a file including a thermal scenario, saved thanks to the saveConfig function.
	*/
    flat_thermal_soaring_zone(std::string filename);

     /// Computes the wind vector w, at point (x,y,z), at time t.
	/**
		\param x a horizontal coordinate of the system (x,y,z).
		\param y a horizontal coordinate of the system (x,y,z).
		\param z the vertical coordinate of the system (x,y,z).
		\param t the time.
		\param w the wind vector: windx, windy, windz.
	*/
    flat_thermal_soaring_zone& wind(double x, double y, double z, double t, std::vector<double> &w);

    ///Simulate a scenario of thermals.
	/**
		\param deltaT the interval of thermal actualization.
		\param model the chosen model for the simulation.
	*/
    void createScenario(double deltaT, int model);

    ///Write the whole wind data for the visualization of a zslice in a file.
	/**
		\param deltaT the interval of thermal actualization.
		\param deltax the definition of the mesh precision onto x-axis.
		\param deltay the definition of the mesh precision onto y-axis.
		\param zslice the height of the windfield you want to write.
		\param filename the name of the file you want to write data.
	*/
    void writeScenario(double deltaT, double deltax, double deltay, double zslice,std::string filename);

    ///Save a thermal scenario in order to play it again in an other simulation.
	/**
		\param filename the name of the .txt file you want to save the scenario.
	*/
    void saveConfig(std::string filename);

    ///Save a thermal scenario in order to play it again in an other simulation.
	/**
		\param filename the name of the .csv file you want to save the scenario.
	*/
    void saveConfigToCSV(std::string filename);

};

}

#endif
