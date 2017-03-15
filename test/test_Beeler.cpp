//
// Created by adrien on 1/18/17.
//
#include <L2F/aircraft/BeelerGlider.hpp>
#include <L2F/flight_zone/flat_thermal_soaring_zone.hpp>

using namespace L2F;

int main(int argc, const char * argv[])
{

    BeelerGlider my_glider;
    int xdim = 6;
    int udim = 3;

    std::vector<double> xdot(xdim);
    std::vector<double> x(xdim);
    std::vector<double> u(udim);

    std::vector<double>(3, 0.).swap(u);

    double time_limit = 1000;
    double deltaT = 100;
    L2F::flat_thermal_soaring_zone my_zone("../data/config.txt");

    // Initilation of the position
    double x0 = 0.;// m
    double y0 = 0.;// m
    double z0 = 100.;// m
    double V0 = 10.; // m/s
    double gamma0 = 0.; // rad
    double khi0 = 0.; // rad

    printf("ok\n");
    printf("Initialisation des variables\n" );



    double time = 0;
    std::vector<double> w(3);

    //my_zone.wind(0., 0, 10., 0.0, w);
    // test sur un element simple
    std::vector<double> state_init = {0,0,100,20,3.14/6,0};
    my_glider.state_dynamics(state_init,u, my_zone,time, xdot);

    // Ici on essaie dans une autre configuration
    state_init = {0,0,100,30,0,0};
    my_glider.state_dynamics(state_init,u, my_zone,time, xdot);

    state_init = {0,0,100,50,3.14/6,0};
    my_glider.state_dynamics(state_init,u, my_zone,time, xdot);

    state_init = {0,0,100,20,3.14/6,0};
    my_glider.state_dynamics(state_init,u, my_zone,time, xdot);

    state_init = {0,0,100,10,3.14/6,0};
    my_glider.state_dynamics(state_init,u, my_zone,time, xdot);

    state_init = {0,0,100,5,3.14/6,0};
    my_glider.state_dynamics(state_init,u, my_zone,time, xdot);

    state_init = {0,0,100,10,0,0};
    my_glider.state_dynamics(state_init,u, my_zone,time, xdot);

    state_init = {0,0,100,10,0.2,0};
    my_glider.state_dynamics(state_init,u, my_zone,time, xdot);

    state_init = {0,0,100,20,0.1,0};
    my_glider.state_dynamics(state_init,u, my_zone,time, xdot);

    state_init = {0,0,100,30,0.1,0};
    my_glider.state_dynamics(state_init,u, my_zone,time, xdot);

    state_init = {0,0,100,40,0.1,0};
    my_glider.state_dynamics(state_init,u, my_zone,time, xdot);

    return(1);

}



