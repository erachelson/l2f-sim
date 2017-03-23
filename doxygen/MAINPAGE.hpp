/*!
  \mainpage l2f-sim


Stand-alone simulator for the "Learning to fly" project

This project was developed with a team of ISAE-SUPAERO students, tutored by Pr. E. Rachelson. The goal is to provide a simple, efficient, stand-alone library, for the simulation of the flight dynamics of an autonomous glider in convective soaring conditions and the development of Reinforcement Learning control algorithms. The code hosted here contains the library and a few demonstration tools allowing to simulate different controls strategies.

The project website can be found at [Learning to Fly](http://websites.isae.fr/learning-to-fly/).

# Launching the simulation

Launch a simulation :
- make all

When you want to explore the results after the simulation :
- make energy (it shows the variation of energy during the simulation)
- make traj (it shows the trajectory of the plane with the thermic centers)
- make angle (its shows the evolution of the aircraft angle during the fligth)

You want to have some help :
- make help (it shows you the help section)

You want to modify the initial parameters :
- Go in the configInit.txt file and modify the parameter you want
(x_init,y_init,z_init,V_init,gamma_init,khi_init) // TODO

You can modify directly in the code those parameter in the /L2F/test/simulation.cpp file.
In this /L2F/test/simulation.cpp file you can also choose what kind on pilot, thermic, stepper, aircraft model you want ...

# Principle (design patterns) of the L2F simulator

## L2F simulations
A L2F simulation is primarily defined via the L2F::simulation class.
This class is defined in the root L2F directory.
It comprises a L2F::flight_zone, a L2F::aircraft, a L2F::pilot and a
L2F::stepper.
Each of these classes is defined in the homonym respective directory.
These classes are interface classes: they define the base attributes
and methods that are called by the simulation. In practice, they need to be
refined by inheritance, in order to define instantiable classes, such as, for
example :
- L2F::BeelerGlider for L2F::aircraft
- L2F::flat_thermal_soaring_zone for L2F::flight_zone
- L2F::euler_integrator for L2F::stepper
- L2F::passive_pilot for L2F::pilot

## Using the library
Example provided in the main.cpp file

## Extending the library
The library can be extended by defining new aircrafts, new steppers, new pilots
or new flight zones. To do so, just create your own class that inherits from the
L2F::aircraft, L2F::stepper, L2F::pilot or L2F::flight_zone base classes (or any
of the already defined classes in the respective directories). Then you can just
place the corresponding files in the respective directories and start using
them.

*/
