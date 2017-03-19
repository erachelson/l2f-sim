#ifndef L2FSIM_FLAT_THERMAL_SOARING_ZONE_HPP_
#define L2FSIM_FLAT_THERMAL_SOARING_ZONE_HPP_

#include <L2Fsim/flight_zone/flight_zone.hpp>
#include <L2Fsim/flight_zone/flat_zone.hpp>
#include <L2Fsim/flight_zone/thermal/std_thermal.hpp>
#include <string>
#include <fstream>
#include <sstream>

namespace L2Fsim {

class flat_thermal_soaring_zone : public flat_zone
{
    
    
    /*--------------------------------------
     ------------- Attributes --------------
     -------------------------------------*/
    
protected:
    
    double tstart;
    double tend;
    
    int minX,maxX;
    int minY,maxY;
    int minZ,maxZ;
    
    double dmin=200.; // dmin represent the radius of a thermal
    
    double zi;           // Height of layer thickness
    
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
    
    // constructor
    flat_thermal_soaring_zone(double Tend,int miX,int maX,int miY,int maY,int miZ,int maZ,double wx,double wy,double zi);
    flat_thermal_soaring_zone(std::string filename);

    // wind
    flat_thermal_soaring_zone& wind(double x, double y, double z, double t, std::vector<double> &w);
    
    // scenario
    void createScenario(double deltaT, int model);
    void writeScenario(double deltaT, double deltax, double deltay, double deltaz,std::string filename);

    // save config
    void saveConfig(std::string filename);
    void saveConfigToCSV(std::string filename);
    
};
    
}

#endif
