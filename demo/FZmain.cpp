#include <L2F/flight_zone/flat_thermal_soaring_zone.hpp>
#include <iostream>


using namespace std;

int main(int argc, const char * argv[])
{
    /*--------------------------------------
     ------------- FLIGHT ZONE -------------
     -------------------------------------*/
    double time_limit = 1000;
    double deltaT     = 100;  // actualization of thermals
    int minX  = -100;         // definition of box
    int maxX  = 100;          // definition of box
    int minY  = -100;         // definition of box
    int maxY  = 100;          // definition of box
    int minZ  = 0;            // definition of box
    int maxZ  = 2000;         // definition of box
    double wx = 0.;           // definition of wind
    double wy = 0.;           // definition of wind
    int zi    = 1400.;        // definition of height of all thermals
    
    // Initialization of a windfield
    // Constructor 1 : flat_zone(Tend,minX,maxX,minY,maxY,minZ,maxZ,wx,wy,zi):tstart=0.
    flat_thermal_soaring_zone my_zone(time_limit,minX,maxX,minY,maxY,minZ,maxZ,wx,wy,zi);
    
    // Constructor 2 : flat_zone(filename)
    // flat_thermal_soaring_zone my_zone("DATA/config1.txt");
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
    
    //flat_zone(Tend,minX,maxX,minY,maxY,minZ,maxZ,dmin,wx,wy,zi):tstart=0.
    L2F::flat_thermal_soaring_zone FZ(1000,-1000,1000,-1000,1000,0,2000,0.,0.,1200);
    
    //flat_zone(filename)
    //L2F::flat_thermal_soaring_zone FZ("../data/config.txt");
    
    // saveConfig2
    //FZ.saveConfig("../data/config2.txt");
    
    // createScenario(deltaT,model);
    // 1 : Allen model
    // 2 : Childress model
    // 3 : Lenschow model with Gaussian distribution
    // 4 : Lenschow with Geodon model
    // 5 : Lawrance model
    FZ.createScenario(deltaT,1);
    
    // writeScenario(deltaT,deltax,deltay,deltaz,filename)
    FZ.writeScenario(deltaT,5.,5.,400.,"/Users/itienne/Desktop/3A/learningtofly/Code/2016-C++_L2F_library/FZ/DATA/wind.txt");
    
    // saveConfig(filename)
    FZ.saveConfig("DATA/config.txt");
}
