#include <L2F/flight_zone/flat_thermal_soaring_zone.hpp>
#include <iostream>


using namespace std;

int main(int argc, const char * argv[])
{
    double deltaT =100.;
    
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
