#ifndef L2F_PILOT_ETIENNE_HPP_
#define L2F_PILOT_ETIENNE_HPP_

#include <L2F/pilot/pilot.hpp>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <L2F/utils/utils.hpp>

namespace L2F {

class pilot_etienne : public pilot
{

protected:
    double masse = 1.35;

	int input_dim;       // All the information we have
	int output_dim;      // All the action we can take

  	std::vector<double> obs_prec;   // previous state
  	std::vector<double> command;    // current command 

  double time=0;

public:
    // Constructor
    pilot_etienne(int dimin = 4,int dimout=3)
    {
        input_dim = dimin;
        output_dim = dimout;
        obs_prec = {1000,20,0,0}; 
        command={0,0,0};
    }

    pilot_etienne& operator()( std::vector<double> &obs,std::vector<double> &newcommand) 
    {
		// Choosing and applying command
		choose_action(obs,obs_prec,newcommand);
		command = newcommand;

		// Save the previous obs
    	obs_prec = obs;
    	double deltaT=0.001;
    	time+=deltaT;

		return *this;
	}

private:
	void choose_action( std::vector<double> &obs, std::vector<double> &obs_prev,std::vector<double> &newcommand){

		double sigma_action =0.5;
		if(this->inThermal(obs,time,"DATA/energy.txt"))
		{
			// on tourne d'un coté
			newcommand[2] = sigma_action;

			std::cout << "--------------"<< std::endl;
			std::cout << "--------------"<< std::endl;
			std::cout << "--------------"<< std::endl;
			std::cout << "--------------"<< std::endl;

			// si vario tombe, fait un 270
				// implementer une loi permettant de savoir si la thermique vaut le coup
				// sinon ==> partir 
		}

		// else == > mode transition / recherche (louvoyage)
	}


	inline bool inThermal(std::vector<double> &obs,double time,std::string filename)
	{
		const static double g = 9.81;
		double h0 = obs_prec[0];
		double h1 = obs[0];
		double V0 = obs_prec[1];
		double V1 = obs[1];
		double dt = 0.001;

		/*
		double energie = masse * (g*h0 + V0*V0/2);
		double vario =  masse * (g*(h1-h0) + V0*(V1-V0));
		double varioLouis = ((h1-h0) + 0.5*(V1*V1-V0*V0)/g) /dt;
		*/

		double dEc = masse * (V1*V1-V0*V0) /2.;
		double dEp = masse * g * (h1-h0);
		double dEtot = dEc + dEp;

		std::ofstream outfile;
        if (time==0.)
        {
        	outfile.open(filename);
        	std::string header="time V0 V1 dEc dEp dEtot";
        	outfile << header << std::endl;
        }
        else{outfile.open(filename,std::ios::app);}

        std::string info=std::to_string(time)+" "
        				+std::to_string(V0)+" "
        				+std::to_string(V1)+" "
        				+std::to_string(dEc)+" "
        				+std::to_string(dEp)+" "
        				+std::to_string(dEtot);

        outfile << info << std::endl;

        outfile.close();

		bool intherm = (dEtot>0)?1:0;

		std::cout << dEtot << std::endl;

		return intherm;
	}

};
}

#endif