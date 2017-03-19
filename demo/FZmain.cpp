#include <L2Fsim/flight_zone/flat_thermal_soaring_zone.hpp>
#include <iostream>


using namespace std;
using namespace L2Fsim;

int main(int argc, const char * argv[])
{
    /*--------------------------------------
     ------------- FLIGHT ZONE -------------
     -------------------------------------*/
    double time_limit = 1000;
    double deltaT     = 100;  // actualization of thermals
    int minX  = -1000;        // definition of box
    int maxX  = 1000;         // definition of box
    int minY  = -1000;        // definition of box
    int maxY  = 1000;         // definition of box
    int minZ  = 0;            // definition of box
    int maxZ  = 2000;         // definition of box
    double wx = 0.;           // definition of wind
    double wy = 0.;           // definition of wind
    int zi    = 1400.;        // definition of height of all thermals
    
    int model = 1;            // definition of the model of thermals
    // 1 : Allen     model
    // 2 : Childress model
    // 3 : Lenschow  model
    // 4 : Geodon    model
    // 5 : Lawrance  model
    
    //------------------------//
    // 1 : create an empty box
    //------------------------//
    flat_thermal_soaring_zone FZ(time_limit,minX,maxX,minY,maxY,minZ,maxZ,wx,wy,zi);
    
    //-----------------------------------//
    // 2 : Simulate a scenario of thermals
    //-----------------------------------//
    FZ.createScenario(deltaT,model);
    
    //---------------------------------------//
    // 3 : write DATA for visualization zslice
    //---------------------------------------//
    // Write the whole wind DATA of a zslice in a file
    double deltax = 5.;    // definition of the mesh precision in x direction
    double deltay = 5.;    // definition of the mesh precision in y direction
    double zslice = 650.;  // height of the windfield you want to write
    FZ.writeScenario(deltaT,deltax,deltay,zslice,"DATA/wind.txt");
    
    //---------------------------//
    // 4 : save this configuration
    //---------------------------//
    FZ.saveConfig("DATA/config2.txt");
    
    //---------------------------------//
    // 5 : call a previous configuration
    //---------------------------------//
    // Constructor 2
    flat_thermal_soaring_zone FZ1("DATA/config1.txt");
    flat_thermal_soaring_zone FZ_one_therm("DATA/configUneSeuleTherm.txt");
    
}
