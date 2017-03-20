#ifndef L2FSIM_PILOT_GENQLEARN_HPP_
#define L2FSIM_PILOT_GENQLEARN_HPP_

#include <L2Fsim/pilot/pilot.hpp>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <L2Fsim/utils/utils.hpp>

namespace L2Fsim {

class pilot_genQlearn : public pilot
{
    public:
    unsigned int output_dim;
    // La dimension de l'observable
    int dim_in;
    // La dimension de la commande
    int dim_out;
    // La dimension de la Qtable
    int dim_Q;
    int dim_Theta;
    // L'incrementation de la Qtable
    double alpha=1;
    // Propagation dans la Qtable
    double gamma=0.99;
    // Le degree d'exploration
    double epsilon=0.01;
    // La Qtable
    std::vector<double> Qtable;
    // Les observables precedants
    std::vector<double> obs_prev;
    std::vector<double> obs_prevprev;
    std::vector<double> action_prev;
    std::vector<double> command_actual;
    // L'incrementation des commandes
    double incrementation_command_sigma=0.003;
    double incrementation_command_beta=0.000;
    double incrementation_command_alpha=0.003;
    double AngleMax = 0.2;
    double limit = 0.5;
    double coeff = -1. / AngleMax * std::log(2./(limit+1.) - 1.);
    double dt = 0.001;

    // Constructor
    pilot_genQlearn(int dimin = 4,int dimout=3,std::vector<double> obsinit={800.,20.,0.,0.}){

        dim_in = dimin;
        dim_out = dimout;
        dim_Q = 4;
        dim_Theta = dim_Q*(dim_Q+1)/2 + dim_Q;
        // Initialisation de la Qtable

        obs_prev = obsinit;
        obs_prevprev = obsinit;
        std::vector<double>(dimout, 0.).swap(action_prev);
        std::vector<double>(dimout, 0.).swap(command_actual);

        std::vector<double>(dim_Theta, 0.).swap(Qtable);
        //Qtable = {-42.136900, 56.290775, 10.637509, -0.017635, 0.004555, 13.542521, 0.265058, 0.891439, 0.018497, 0.013512, -0.003882};


    }
    // Operator pour avoir la commande
    pilot_genQlearn& operator()(std::vector<double> &obs,std::vector<double> &command){
        // Ici on actualise la Q table en fonction des resultats
        
        // We print the Qtable
        print_Qtable();
        
        actualize_Q_table(obs);

        // On choisit l'action que l'on va prendre
        choose_action(obs,command);
        
        // Save the previous obs
        obs_prevprev = obs_prev;
        obs_prev = obs;
        command_actual = command;

        return *this;
    }

    // Generate the true vector
    void get_vQ(std::vector<double> &obs,std::vector<double> &obs_prec,std::vector<double> &command,std::vector<double> &vQ){
        // intput obs : h, V, gamma, Khi
        // input command : alpha,beta, sigma
        // True input : dh, dgamma,
        // ouput : produit - couplage + terme simple
        int i,j;

        double dh = (obs.at(0) - obs_prec.at(0))/dt;
        double dV = (obs.at(1) - obs_prec.at(1))/dt;
        double dgamma = (obs.at(2) - obs_prec.at(2))/dt;
        double beta_ = command.at(1);
        double alpha_ = command.at(0);
        double sigma_ = command.at(2);
        double gamma_angle = obs.at(2);

        std::vector<double> concat = {sigmoid( dgamma, 0.05, 0) - 0.5,sigmoid( dh, 2, 0) - 0.5,sigmoid( sigma_, 0.3, 0) - 0.5,sigmoid( alpha_, 0.1, 0) - 0.5};

        if (vQ.size() != 0) {
            printf("size vQ = %d : error\n",vQ.size());
            exit(-1);
        }
        
        //vQ = {1.};
        
        // Boucle sur les elements simples
        for(i = 0;i<concat.size();i++){
            vQ.push_back(concat.at(i));
        }

        for(i = 0;i<concat.size();i++){
            for(j=0;j<i+1;j++){
                
                    vQ.push_back(concat.at(i)*concat.at(j));
                
            }
        }
        
    }

    // Get Q value
    void get_Q_value(std::vector<double> &obs,std::vector<double> &obs_prec,std::vector<double> &command,double &Q){

        Q = 0;
        std::vector<double> vQ;
        get_vQ(obs,obs_prec,command,vQ);

        int i;
        for(i=0;i<vQ.size();i++){
            
                Q += Qtable.at(i)*vQ.at(i);
            
            
        }
    }
    // Get dQ value
    void get_dQ_value(std::vector<double> &obs,std::vector<double> &obs_prec,std::vector<double> &command,int &var,double &dQ){

        std::vector<double> vQ;
        get_vQ(obs,obs_prec,command,vQ);
        dQ = vQ.at(var);
    }

    // Actualize Qtable
    void actualize_Q_table(std::vector<double> &obs){

        double reward;
        get_reward(obs,obs_prev,reward);
        double maxQ;
        double Q;
        double dQ;
        int n = Qtable.size();
        int i;
        
        get_max_action_obs(obs,obs_prev,maxQ);
        get_Q_value(obs_prev,obs_prevprev,command_actual,Q);
        
        std::vector<double> vQ;
        get_vQ(obs_prev,obs_prevprev,command_actual,vQ);
        
        printf("begin print parameters\n");
        printf("maxQ : %lf\n",maxQ);
        printf("reward : %lf\n",reward);
        printf("Q : %lf\n",Q);
        /*
        printf("dgamma : %lf\n",vQ.at(0));
        printf("dh : %lf\n",vQ.at(1));
        printf("sigma : %lf\n",vQ.at(2));
        printf("obs, h : %lf\n",obs.at(0));
        printf("obs, V : %lf\n",obs.at(1));
        */
        std::vector<double> newQtable;
        std::vector<double>(dim_Theta, 0.).swap(newQtable);
        
        for(i=0;i<n;i++){

            dQ = vQ.at(i);
            newQtable.at(i) = Qtable.at(i) + alpha*(reward + gamma*maxQ - Q)*dQ;
            /*

            printf("dQ : %lf\n",dQ);
            */
            if(std::abs(Qtable.at(i)) > 10000000){
                print_Qtable();
                printf("The Qtable have too high value\n");
                exit(-1);
            }
        }
        Qtable = newQtable;

    }

    void choose_action(std::vector<double> &obs,std::vector<double> &command){

        std::vector<double> action_tot_alpha;
        std::vector<double> action_tot_beta;
        std::vector<double> action_tot_sigma;
        
        get_possible_action(obs,action_tot_alpha,action_tot_beta,action_tot_sigma);
        
        double rand_x = rand_double(0,1);

        if(rand_x < epsilon){

            std::cout << "Exploration." << std::endl;
            // Exploration
            int action_alpha =  rand_int(0,action_tot_alpha.size());
            action_prev[0] = action_tot_alpha.at(action_alpha);

            int action_beta =  rand_int(0,action_tot_beta.size());
            action_prev[1] = action_tot_beta.at(action_beta);

            int action_sigma =  rand_int(0,action_tot_sigma.size());
            action_prev[2] = action_tot_sigma.at(action_sigma);

        }else{

            std::cout << "Exploitation" << std::endl;
            // Exploitation
            get_max_action_action(obs,obs_prev,action_prev);

        }
        
        command = {action_prev[0] + command_actual.at(0),action_prev[1] + command_actual.at(1),action_prev[2] + command_actual.at(2)};
        
        //printf("choice : alpha = %lf, beta = %lf, sigma = %lf\n",action_prev[0],action_prev[1],action_prev[2]);

    }

    void get_max_action_obs(std::vector<double> &obs,std::vector<double> &obs_prec,double &maxQ){

        // On enumére toutes les actions possibles
        std::vector<double> action_tot_alpha;
        std::vector<double> action_tot_beta;
        std::vector<double> action_tot_sigma;
        
        get_possible_action(obs,action_tot_alpha,action_tot_beta,action_tot_sigma);
        
        int i,j,p;
        maxQ = -10000000000;
        double Q;
        std::vector<double> command_a;

        /*
        printf("Size alpha : %d\n",action_tot_alpha.size());
        printf("Size beta : %d\n",action_tot_beta.size());
        printf("Size sigma : %d\n",action_tot_beta.size());
        */
        
        if(action_tot_alpha.size() != 3 &&action_tot_beta.size() != 3&&action_tot_beta.size() != 3){
            printf("Size not ok\n");
            exit(-1);
        }
        
        //printf("choice of maxQ\n");
        for(i=0;i<action_tot_alpha.size();i++){
            for(j=0;j<action_tot_beta.size();j++){
                for(p=0;p<action_tot_sigma.size();p++) {
                    //printf("choice\n");
                    command_a ={command_actual.at(0) + action_tot_alpha.at(i),command_actual.at(1) + action_tot_beta.at(j),command_actual.at(2) + action_tot_sigma.at(p)};
                    get_Q_value(obs,obs_prec, command_a, Q);
                    if (Q > maxQ) {
                        maxQ = Q;
                    }
                }
            }
        }


    }


    void get_max_action_action(std::vector<double> &obs,std::vector<double> &obs_prec,std::vector<double> &action){

        // On enumére toutes les actions possibles
        // On enumére toutes les actions possibles
        std::vector<double> action_tot_alpha;
        std::vector<double> action_tot_beta;
        std::vector<double> action_tot_sigma;
        
        get_possible_action(obs,action_tot_alpha,action_tot_beta,action_tot_sigma);
        
        int i,j,p;
        double max_Q = -10000000000;
        
        if(action_tot_alpha.size() != 3 &&action_tot_beta.size() != 3&&action_tot_beta.size() != 3){
            printf("Size not ok\n");
            exit(-1);
        }
        
        
        
        double Q;
        std::vector<double> command_a;
        //printf("choice of action\n");
        for(i=0;i<action_tot_alpha.size();i++){
            for(j=0;j<action_tot_beta.size();j++){
                for(p=0;p<action_tot_sigma.size();p++) {
                    //printf("choice\n");
                    command_a ={command_actual.at(0) + action_tot_alpha.at(i),command_actual.at(1) + action_tot_beta.at(j),command_actual.at(2) + action_tot_sigma.at(p)};
                    get_Q_value(obs,obs_prec, command_a, Q);
                    if (Q > max_Q) {
                        max_Q = Q;
                        action = {action_tot_alpha.at(i), action_tot_beta.at(j), action_tot_sigma.at(p)};
                    }
                    
                    
                }
            }
        }
        

    }

    // Calculation of the reward
    void get_reward(std::vector<double> &obs,std::vector<double> &obs_prec,double &reward){

        // Obs[0] -> la hauteur h
        // Obs[1] -> la vitesse V
        
        reward = ((obs[0] - obs_prec[0]) + 1/(2*9.81)*(obs[1]*obs[1] - obs_prec[1]*obs_prec[1]))/dt;

    }
    
    void get_possible_action(std::vector<double> &obs,std::vector<double> &action_tot_alpha,std::vector<double> &action_tot_beta,std::vector<double> &action_tot_sigma){
    
        
        action_tot_beta = {-incrementation_command_beta,0.,incrementation_command_beta};
        
        // sigma condition
        action_tot_sigma = {-incrementation_command_sigma,0.,incrementation_command_sigma};
        if (command_actual.at(2) >= 0.6) {
            action_tot_sigma = {-2*incrementation_command_sigma,-incrementation_command_sigma,0.};
        }
        
        if (command_actual.at(2) < -0.6) {
            action_tot_sigma = {2*incrementation_command_sigma,incrementation_command_sigma,0.};
        }
        
        
        // If everything is ok
        if(std::abs(obs.at(2) + command_actual.at(0)) < 0.3){
            
        //printf("no problem for command\n");
        action_tot_alpha = {-incrementation_command_alpha,0.,incrementation_command_alpha};
        
        
        
        } else {
            
            //printf("problem for command\n");
            // If not
            
            if (obs.at(2) + command_actual.at(0) >= 0.3) {
                action_tot_alpha = {-4*incrementation_command_alpha,-2*incrementation_command_alpha};
            }
            
            if (obs.at(2) + command_actual.at(0) <= -0.3) {
                action_tot_alpha = {4*incrementation_command_alpha,2*incrementation_command_alpha};
            }
            

        }
        
        
    }
    
    // Function use to print the Qtable
    void print_Qtable(){
        
        printf("Qtable : ");
        int i;
        for (i=0; i<Qtable.size(); i++) {
            printf(" %lf",Qtable[i]);
        }
        printf("\n");
        
    }
    
    pilot_genQlearn& out_of_range(std::vector<double> &obs,
                          std::vector<double> &command)
    {
        std::vector<double>(3, 0.).swap(command);
        command.at(2) = 0.4;
        
        obs_prevprev = obs_prev;
        obs_prev = obs;
        command_actual = command;
    }
    

};
}
#endif
