#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef L2F_UTILS_HPP_
#define L2F_UTILS_HPP_

#include <random>
#include <fstream>
#include <iostream>

namespace L2F {

    inline double sgn(double T) {
        if (T > 0){ return 1.0;}
        if (T < 0){ return -1.0;}
        if (T == 0){ return 1.0;}
        return 0;
    }
    /*--------------------------------------
     -------------   RANDOM   --------------
     -------------------------------------*/

    // return a double between double "a" and double "b"
    inline double rand_double(double a, double b)
    {
        return ( rand()/(double)RAND_MAX ) * (b-a) + a;
    }

    // return an int between int "a" and int "b"
    // implicitly a<b
    inline int rand_int(int a, int b)
    {
        return rand()%(b-a) +a;
    }

    // produce a number following a normal distribution
    inline double normalLaw()
    {
        double W,V1,V2;
        do
        {
            double U1=rand_double(0., 1.);
            double U2=rand_double(0., 1.);
            V1=2*U1-1;
            V2=2*U2-1;
            W=V1*V1+V2*V2;
        }while(W>1);
        double W_function=sqrt(-2*log(W)/W);
        return V1*W_function;
    }

    /*
     * Dans cette partie on sauvegarde les états (states) de l'aircraft dans un fichier pour plotter plus tard
     */
    inline void save_state_into_file(std::string filename,std::vector<double> &state,double time)
    {
        std::ofstream outfile;
        if (time==0.){outfile.open(filename);}
        else{outfile.open(filename,std::ios::app);}
        // La taille d'un vecteur
        int size = state.size();
        
        std::string info="";

        for(int i = 0;i<size;i++)
        {
            info = info + std::to_string(state.at(i)) + " ";
        }

        info = info + std::to_string(time);

        outfile << info << std::endl;

        outfile.close();
    }

    inline void save_energy_into_file(std::string filename,double time,double dt,
                                       std::vector<double> &obs_prec,
                                       std::vector<double> &obs)
    {
        const static double g = 9.81;
        double h0 = obs_prec[0];
        double h1 = obs[0];
        double V0 = obs_prec[1];
        double V1 = obs[1];
        /*
        double energie = masse * (g*h0 + V0*V0/2);
        double varioEtienne =  masse * (g*(h1-h0) + V0*(V1-V0));
        double varioAdrien = g * masse *(h1 - h0) + *masse/2*(V1*V1 - V0*V0);
        double varioLouis = ((h1-h0) + 0.5*(V1*V1-V0*V0)/g) /dt;
        */

        double dEc = ((V1*V1-V0*V0) / (2.*g)) /dt;
        double dEp =  (h1-h0)/dt ;
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
    }


    inline void save_Qtable(std::string filename,std::vector<double> &state,std::vector<double> & Qtable)
    {
        std::ofstream outfile;
        outfile.open(filename);
        
        // La taille d'un vecteur
        int size = state.size();
        
        std::string info = "";
        
        for(int i = 0;i<size;i++)
        {
            info = info + std::to_string(Qtable.at(i)) + " ";
            if (i==size-1)
            {
                info = info + std::to_string(Qtable.at(i)) + "  ";
            }
        }
        
        outfile << info << std::endl;
        
        outfile.close();
    }

    inline void load_Qtable(std::string filename,std::vector<double> & Qtable){
        std::ifstream file(filename, std::ios::in);
        std::string lign;
        getline(file,lign);
        while(lign!=" ")
        {
            unsigned long cut = lign.find(" ");
            unsigned long rest =lign.size()-cut;
            double num = stod(lign.substr(0,cut));
            lign = lign.substr(cut+1,rest);
            Qtable.push_back(num);
        }
    }

    inline double sigmoid(double x,double a,double c)
    {
        double exp_value;
        float return_value;

        /*** Exponential calculation ***/
        exp_value = exp( -a*(x-c));

        /*** Final sigmoid value ***/
        return_value = 1 / (1 + exp_value);

        return return_value;
    }

}

#endif
